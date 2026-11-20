# Scientific definitions: SurvPT point predictions never sample the posterior.
suppressPackageStartupMessages({library(SurvPT); library(interval); library(spBayesSurv)})
stopifnot(as.character(packageVersion('SurvPT')) == '0.4.2',
          as.character(packageVersion('spBayesSurv')) == '1.1.9')
for (file in c('Fiducial.R', 'LinInterpolation.R', 'Simulation_MPT.R', 'Simulation_TBP.R', 'Winkler.R'))
    source(file.path(.simu_dir, file), local = TRUE)

all_settings <- function() {
    rows <- list()
    for (n in c(50L, 75L, 100L, 200L, 500L, 800L, 1000L)) for (type in 1:4)
        rows[[length(rows) + 1L]] <- data.frame(suite = 'fiducial', n = n,
            type = as.character(type), K = NA_real_, Lambda = NA_real_, q = NA_real_,
            label = sprintf('fiducial_s%s_n%s', type, n))
    # K=1 is exponential; retain the original setting ordering and data seeds.
    for (n in c(200L, 500L)) for (type in c('right', 'current', 'interval'))
        for (q in c(.3, .5)) for (K in c(1, 2)) {
            suite <- if (K == 1) 'exp' else 'weibull'
            rows[[length(rows) + 1L]] <- data.frame(suite = suite, n = n, type = type,
                K = K, Lambda = .12, q = q,
                label = paste0(suite, '_', type, '_n', n, '_lambda0.12_q', q))
        }
    # Append new settings: historical setting IDs and data seeds remain unchanged.
    for (n in c(200L, 500L)) for (type in c('right', 'current', 'interval'))
        for (q in c(.3, .5)) rows[[length(rows) + 1L]] <- data.frame(
            suite = 'mixture', n = n, type = type, K = NA_real_, Lambda = NA_real_, q = q,
            label = paste0('mixture_', type, '_n', n, '_q', q))
    ans <- do.call(rbind, rows)
    ans$setting_id <- seq_len(nrow(ans))
    ans
}

pt_methods <- function(sensitivity = FALSE, gamma_sensitivity = FALSE) {
    rows <- list()
    for (prior in c(1L, 2L, 3L)) for (flat in c(FALSE, TRUE)) {
        name <- paste0(c('1' = 'Unif', '2' = 'Exp', '3' = 'Weibull_PT')[[as.character(prior)]],
                       if (flat) '_Flat' else '')
        rows[[length(rows) + 1L]] <- data.frame(method = name, prior = prior,
            flat = flat, kappa = 0.1, sensitivity = FALSE, gamma = 0.5, gamma_sensitivity = FALSE)
    }
    main <- do.call(rbind, rows)
    if (sensitivity) for (kappa in c(0.05, 0.1, 0.2, 0.3)) {
        extra <- main
        extra$method <- paste0(main$method, '_kappa', kappa)
        extra$kappa <- kappa
        extra$sensitivity <- TRUE
        rows[[length(rows) + 1L]] <- extra
    }
    if (gamma_sensitivity) for (gamma in c(0.5, 1, 2)) {
        extra <- main[!main$flat, , drop = FALSE]
        extra$method <- paste0(extra$method, '_gamma', gamma)
        extra$gamma <- gamma
        extra$gamma_sensitivity <- TRUE
        rows[[length(rows) + 1L]] <- extra
    }
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    out
}

canonical_method <- function(method) sub('_(kappa0[.]1|gamma0[.]5)$', '', method)

scientific_config <- function() {
    defaults <- function(f) lapply(formals(f)[-c(1L, 2L, 3L)], eval)
    list(max_M = 8L, max_EM_iter = 2000L, EM_tol = 1e-6, level = .95,
        survpt_defaults = lapply(formals(SurvPT::survPT)[
            c('allocated', 'gamma', 'kappa', 's0', 'min_count', 'split_epsilon')], eval),
        sensitivity_kappa = c(0.05, 0.1, 0.2, 0.3),
        sensitivity_gamma = c(0.5, 1, 2),
        mixture = list(weights = c(.5, .5), shape = c(3, 4), scale = c(3, 10),
            generator = 'SurvPT::data_mixture(seed, n, type, q)', quantile_tol = 1e-12),
        weibull_prior = 'SurvPT::estimate(L, R); c(inverse_scale, shape), fixed during PT sampling',
        survpt_sampling = lapply(formals(getS3method('predict', 'survPT'))[
            c('n_samples', 'burn_in', 'thin')], eval),
        MPT = defaults(spBayesSurv::survregbayes2), TBP = defaults(spBayesSurv::survregbayes),
        fiducial = lapply(formals(fiducial)[c('mfid', 'mburn')], eval),
        fiducial_evaluation = list(point = 'original linear-interpolation sample median',
            posterior_mean = 'mean of linear-interpolation survival samples',
            interval = 'equal-tailed fiducial interval; original order-statistic convention'),
        metrics = c('SSE', 'RMSE', 'RMSE_nolast'),
        probability_bounds = list(tolerance = 1e-8,
            scope = 'posterior estimate, mean, lower, upper',
            action = 'clip finite boundary roundoff only; retain point-grid values'),
        dataset_summary = 'package tau/proportion attributes; true population S(mean(tau)); all datasets',
        posterior_estimate = 'PT: predict estimate (MAP-EM); Fiducial: original median; MPT/TBP: Shat',
        posterior_bias = 'separate mean(abs(estimate-p)) and mean(abs(posterior_mean-p)); independent valid counts',
        survival_probabilities = c(.75, .5, .25), grid_points = 201L,
        seed_rule = 'data=replicate; method=100000000+1000000*method_id+10000*setting_id+replicate',
        bayes_adapter = 'Existing baseline adapters: zero MPT design; deterministic singular initializer recovery')
}

mixture_survival <- function(x) .5 * pweibull(x, 3, 3, lower.tail = FALSE) +
    .5 * pweibull(x, 4, 10, lower.tail = FALSE)
true_survival <- function(x, s) {
    if (s$suite == 'mixture') return(mixture_survival(x))
    if (s$suite == 'fiducial') return(pgamma(x,
        shape = switch(s$type, '1' = 1, '2' = 3, '3' = 2, '4' = 1),
        rate = 1, lower.tail = FALSE))
    pweibull(x, shape = s$K, scale = 1 / s$Lambda, lower.tail = FALSE)
}
mixture_quantiles <- function(probabilities = c(.75, .5, .25)) {
    vapply(probabilities, function(p) {
        stopifnot(is.finite(p), p > 0, p < 1)
        upper <- 10
        while (mixture_survival(upper) > p) upper <- upper * 2
        uniroot(function(x) mixture_survival(x) - p, c(0, upper), tol = 1e-12)$root
    }, numeric(1))
}

generate_input <- function(t, s) {
    if (s$suite == 'fiducial') {
        d <- getExportedValue('SurvPT', paste0('fiducial_', s$type))(t, s$n)
        shape <- switch(s$type, '1' = 1, '2' = 3, '3' = 2, '4' = 1)
        survival <- function(x) pgamma(x, shape = shape, rate = 1, lower.tail = FALSE)
        quantiles <- qgamma(c(.25, .5, .75), shape = shape, rate = 1)
    } else if (s$suite == 'mixture') {
        d <- SurvPT::data_mixture(seed = t, n = s$n, type = s$type, q = s$q)
        survival <- mixture_survival
        quantiles <- mixture_quantiles()
    } else {
        generator <- getExportedValue('SurvPT', switch(s$type,
            right = 'data_rc', current = 'data_cs', interval = 'data_ic'))
        d <- generator(t, s$n, s$K, s$Lambda, s$q)
        survival <- function(x) pweibull(x, s$K, 1 / s$Lambda, lower.tail = FALSE)
        quantiles <- qweibull(c(.25, .5, .75), s$K, 1 / s$Lambda)
    }
    endpoints <- as.matrix(d[, c('L', 'R')])
    tau <- attr(d, 'tau', exact = TRUE)
    proportions <- attr(d, 'censoring_proportions', exact = TRUE)
    stopifnot(length(tau) == 1L, is.finite(tau), tau > 0,
        identical(tau, max(endpoints[is.finite(endpoints)])),
        identical(names(proportions), c('Exact', 'Left', 'Interval', 'Right')),
        all(is.finite(proportions)), all(proportions >= 0 & proportions <= 1),
        abs(sum(proportions) - 1) < 1e-12,
        isTRUE(all.equal(unname(proportions), c(
            mean(is.finite(d$R) & d$L == d$R),
            mean(is.finite(d$R) & d$L == 0 & d$L < d$R),
            mean(is.finite(d$R) & d$L > 0 & d$L < d$R),
            mean(is.infinite(d$R))), tolerance = 1e-12)),
        isTRUE(all.equal(attr(d, 'true_S_tau'), true_survival(tau, s), tolerance = 1e-12)))
    x <- seq(0, tau, length.out = 201L)
    list(data = d, tau = tau, censoring_proportions = proportions, grid = x, truth = survival(x),
         quantiles = quantiles, probabilities = c(.75, .5, .25))
}

both_errors <- function(estimated, truth) {
    # Preserve EMPT's existing probability-increment RMSE definitions.
    c(SSE = sum((estimated - truth)^2),
      RMSE = sqrt(mean((abs(diff(c(estimated, 0))) - abs(diff(c(truth, 0))))^2)),
      RMSE_nolast = sqrt(mean((abs(diff(estimated)) - abs(diff(truth)))^2)))
}

check_curve <- function(x) {
    if (any(!is.finite(x))) stop('Non-finite point prediction (NA/NaN/Inf)')
    stopifnot(length(x) == 201L, all(is.finite(x)), all(x >= -1e-8 & x <= 1 + 1e-8),
              all(diff(x) <= 1e-8))
    invisible(TRUE)
}

normalize_probability <- function(x, label = 'probability', tolerance = 1e-8) {
    stopifnot(is.numeric(x), length(tolerance) == 1L, is.finite(tolerance), tolerance >= 0)
    finite <- is.finite(x)
    if (any(x[finite] < -tolerance | x[finite] > 1 + tolerance))
        stop(label, ' outside [0,1] by more than tolerance ', format(tolerance), call. = FALSE)
    # Assign only out-of-range entries: valid numbers and RNG state are unchanged.
    x[finite & x < 0] <- 0
    x[finite & x > 1] <- 1
    x
}

posterior_table <- function(input, mean, lower, upper, level = .95, allow_missing = FALSE, estimate = mean) {
    stopifnot(length(mean) == 3L, length(lower) == 3L, length(upper) == 3L, length(estimate) == 3L)
    missing <- !is.finite(mean) | !is.finite(lower) | !is.finite(upper)
    if ((any(missing) || any(!is.finite(estimate))) && !allow_missing) stop('Non-finite posterior prediction (NA/NaN/Inf)')
    valid <- !missing
    estimate[!is.finite(estimate)] <- NA_real_
    stopifnot(all(lower[valid] <= upper[valid]))
    estimate <- normalize_probability(estimate, 'estimate')
    mean[valid] <- normalize_probability(mean[valid], 'posterior mean')
    lower[valid] <- normalize_probability(lower[valid], 'posterior lower')
    upper[valid] <- normalize_probability(upper[valid], 'posterior upper')
    # Keep valid point estimates and targets when another posterior target is missing.
    out <- data.frame(probability = input$probabilities, time = input$quantiles,
        mean = NA_real_, lower = NA_real_, upper = NA_real_, absolute_bias = NA_real_,
        length = NA_real_, covered = NA_integer_, winkler_score = NA_real_)
    if (any(valid)) out[valid, ] <- data.frame(probability = input$probabilities[valid],
        time = input$quantiles[valid], mean = as.numeric(mean[valid]),
        lower = as.numeric(lower[valid]), upper = as.numeric(upper[valid]),
        absolute_bias = abs(as.numeric(mean[valid]) - input$probabilities[valid]),
        length = as.numeric(upper[valid] - lower[valid]),
        covered = as.integer(lower[valid] <= input$probabilities[valid] & input$probabilities[valid] <= upper[valid]),
        winkler_score = winkler_interval_score(input$probabilities[valid], lower[valid], upper[valid], level))
    out$estimate <- as.numeric(estimate)
    out$estimate_absolute_bias <- abs(out$estimate - input$probabilities)
    out
}

method_seed <- function(t, setting_id, method) {
    # Preserve the original 20 PT and five comparator IDs when extending the grid.
    spec <- pt_methods(TRUE)
    catalog <- c(spec$method[spec$prior != 3L], 'OPT', 'Fiducial', 'NPMLE', 'MPT', 'TBP',
                 spec$method[spec$prior == 3L])
    gamma_spec <- pt_methods(gamma_sensitivity = TRUE)
    catalog <- c(catalog, gamma_spec$method[gamma_spec$gamma_sensitivity])
    id <- match(canonical_method(method), catalog)
    stopifnot(!is.na(id), t <= 9999L, setting_id <= 999L)
    as.integer(100000000 + 1000000 * id + 10000 * setting_id + t)
}

expected_methods <- function(s, analysis) {
    stopifnot(analysis %in% c('all', 'main', 'sensitivity', 'gamma'))
    main <- c(pt_methods()$method, if (s$type == 'right') 'OPT', 'Fiducial', 'NPMLE', 'MPT', 'TBP')
    spec <- pt_methods(TRUE, TRUE)
    list(main = if (analysis %in% c('all', 'main')) main else character(),
         sensitivity = if (analysis %in% c('all', 'sensitivity')) spec$method[spec$sensitivity] else character(),
         gamma_sensitivity = if (analysis %in% c('all', 'gamma')) spec$method[spec$gamma_sensitivity] else character())
}

fit_method_strict <- function(method, input, s, t, config) {
    seed <- method_seed(t, s$setting_id, method)
    set.seed(seed)
    warnings <- character()
    start <- proc.time()[['elapsed']]
    result <- withCallingHandlers({
        posterior <- sampling <- prior <- convergence <- NULL
        issue <- list()
        spec <- pt_methods(TRUE, TRUE)
        spec <- spec[spec$method == method, , drop = FALSE]
        if (nrow(spec)) {
            # SurvPT 0.4.2 auto-estimates prior=2 only; prior=3 requires explicit MLEs.
            prior_params <- if (spec$prior == 3L) SurvPT::estimate(input$data$L, input$data$R) else NULL
            fit <- SurvPT::survPT(input$data[, c('L', 'R')], prior = spec$prior,
                prior_params = prior_params, flat = spec$flat, kappa = spec$kappa, gamma = spec$gamma,
                max_M = config$max_M, max_EM_iter = config$max_EM_iter,
                EM_tol = config$EM_tol, verbose = FALSE)
            prior <- fit$Prior
            convergence <- fit$Convergence
            curve <- as.numeric(predict(fit, input$grid, sample = FALSE)$estimate)
            # Default 2000 retained / 3000 burn-in; sample only at true quantiles.
            p <- predict(fit, input$quantiles, level = config$level)
            posterior <- posterior_table(input, p$mean, p$lower, p$upper,
                level = config$level, allow_missing = TRUE, estimate = p$estimate)
            if (anyNA(posterior[c('estimate', 'mean')])) issue <- posterior_failure_info(method, input, posterior)
            sampling <- attributes(p)[c('n_samples', 'burn_in', 'thin')]
        } else if (method == 'OPT') {
            stopifnot(s$type == 'right')
            fit <- SurvPT::survOPT(input$data$L, input$data$status, test_seq = input$grid,
                max_M = config$max_M, max_iter = config$max_EM_iter,
                iter_tol = config$EM_tol, verbose = FALSE)
            curve <- as.numeric(fit$Survival_Probs)
        } else if (method == 'Fiducial') {
            grid <- if (s$suite == 'fiducial') seq(0, if (s$type == '4') 3 else 5,
                                                  length.out = 100L) else input$grid
            times <- sort(unique(c(input$grid, input$quantiles)))
            fit <- fiducial(input$data$L, input$data$R, s$n, grid, times,
                mfid = config$fiducial$mfid, mburn = config$fiducial$mburn, level = config$level)
            curve <- as.numeric(1 - fit$point_li[match(input$grid, times)])
            ix <- match(input$quantiles, times)
            # The sampler returns CDF summaries: survival reverses the interval endpoints.
            posterior <- posterior_table(input, 1 - fit$mean[ix],
                1 - fit$upper[ix], 1 - fit$lower[ix], level = config$level, estimate = 1 - fit$point_li[ix])
            sampling <- config$fiducial
        } else if (method == 'NPMLE') {
            fit <- interval::icfit(input$data$L, input$data$R)
            curve <- as.numeric(interval::getsurv(input$grid, fit, 'interpolation')[[1L]]$S)
        } else if (method %in% c('MPT', 'TBP')) {
            d <- input$data[, c('L', 'R')]
            d$L[d$L == 0 & d$R > 0] <- NA_real_
            d$R[is.infinite(d$R)] <- NA_real_
            fit <- if (method == 'MPT') fit_mpt(d) else fit_tbp(d)
            times <- sort(unique(c(input$grid[input$grid > 0], input$quantiles)))
            p <- if (method == 'MPT') plot(fit, xpred = matrix(0, 1L, 1L),
                tgrid = times, CI = config$level, PLOT = FALSE) else
                plot(fit, tgrid = times, CI = config$level, PLOT = FALSE)
            curve <- c(1, as.numeric(p$Shat[match(input$grid[-1L], times), 1L]))
            ix <- match(input$quantiles, times)
            posterior <- posterior_table(input, p$Shat[ix, 1L], p$Shatlow[ix, 1L], p$Shatup[ix, 1L],
                                         level = config$level, estimate = p$Shat[ix, 1L])
            sampling <- fit$mcmc
        } else stop('Unknown method: ', method)
        check_curve(curve)
        c(list(curve = curve, errors = both_errors(curve, input$truth), posterior = posterior,
             sampling = sampling, prior = prior, convergence = convergence), issue)
    }, warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart('muffleWarning')
    })
    c(result, list(seed = seed, seconds = proc.time()[['elapsed']] - start, warnings = unique(warnings)))
}


# Protect method fitting, prediction and result validation independently of peers.
fit_method <- function(method, input, s, t, config) {
    if (!method %in% method_error_policy$methods)
        return(fit_method_strict(method, input, s, t, config))
    start <- proc.time()[['elapsed']]
    captured_warnings <- character()
    tryCatch(withCallingHandlers({
        v <- fit_method_strict(method, input, s, t, config)
        validate_fit(v, method, input, s, t, config)
        v
    }, warning = function(w) captured_warnings <<- c(captured_warnings, conditionMessage(w))),
    error = function(err) {
        v <- skipped_method_fit(method, input, s, t, err, proc.time()[['elapsed']] - start)
        v$warnings <- unique(captured_warnings)
        v
    })
}

simulate_replicate <- function(t, s, analysis, config) {
    input <- generate_input(t, s)
    methods <- expected_methods(s, analysis)
    values <- list()
    for (method in unique(unlist(methods, use.names = FALSE))) {
        # Reuse kappa=0.1 and gamma=0.5 main fits and posterior streams exactly.
        actual <- canonical_method(method)
        if (!actual %in% names(values)) values[[actual]] <- fit_method(actual, input, s, t, config)
        if (actual != method) values[[method]] <- values[[actual]]
    }
    list(replicate = t, setting_id = s$setting_id, input = input,
         data_hash = digest::digest(input, algo = 'sha256', serializeVersion = 3),
         main = values[methods$main], sensitivity = values[methods$sensitivity],
         gamma_sensitivity = values[methods$gamma_sensitivity])
}
