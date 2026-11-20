#!/usr/bin/env Rscript
self <- sub('^--file=', '', grep('^--file=', commandArgs(), value = TRUE)[1L])
e <- new.env(parent = globalenv()); source(file.path(dirname(normalizePath(self)), 'Pipeline.R'), local = e)
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[1L] else tempfile('verify-survpt-errors-')
s <- e$all_settings()[29L, ]; config <- e$scientific_config()
# Test mixed success, errors, nonfinite targets, and validation failures in all groups.
run <- file.path(out, 'run')
e$prepare_run(run, replicates = 3L, chunk_size = 2L, settings = s)
f <- new.env(parent = globalenv()); source(file.path(run, 'source', 'Pipeline.R'), local = f)
fixture <- function(method, input, s, t, config) {
    if (method == 'Unif' || method == 'Unif_kappa0.2' || method == 'Exp_gamma2') stop('injected SurvPT fit error')
    if (method == 'OPT' && t == 1L) stop('injected OPT failure')
    spec <- pt_methods(TRUE, TRUE); spec <- spec[spec$method == method, , drop = FALSE]
    estimate <- input$probabilities + c(.01, -.02, .03)
    mean <- input$probabilities + c(-.04, .05, -.06)
    if (method == 'Exp' && t == 1L) mean[3] <- NA_real_
    if (method == 'Exp' && t == 2L) estimate[2] <- Inf
    posterior <- if (method == 'OPT' || method == 'NPMLE') NULL else
        posterior_table(input, mean, mean - .1, mean + .1, allow_missing = TRUE, estimate = estimate)
    v <- list(curve = input$truth, errors = both_errors(input$truth, input$truth), posterior = posterior,
        sampling = if (nrow(spec)) config$survpt_sampling else if (method == 'Fiducial') config$fiducial else NULL,
        prior = if (nrow(spec)) list(family = c('uniform', 'exponential', 'weibull')[spec$prior],
            flat = spec$flat, kappa = spec$kappa, gamma = spec$gamma, parameters = c(.5, 2)) else NULL,
        seed = method_seed(t, s$setting_id, method), seconds = 0, warnings = character())
    if (!is.null(posterior) && anyNA(posterior[c('mean', 'estimate')]))
        v <- c(v, posterior_failure_info(method, input, posterior))
    if (method == 'Weibull_PT' && t == 1L) v$curve[10L] <- NA_real_
    if (method == 'Weibull_PT_Flat' && t == 2L) v$sampling$thin <- -1L
    v
}
environment(fixture) <- f; f$fit_method_strict <- fixture
m <- f$read_run(run)
f$report_run_errors(run)
stopifnot(sum(read.csv(file.path(run, 'fit_status_summary.csv'))$attempted) == 0L)
for (i in seq_len(nrow(m$tasks))) f$run_task(run, i, cores = 2L)
f$merge_run(run)
p <- read.csv(file.path(run, 'posterior', s$label, 'Exp.csv'))
summary <- read.csv(file.path(run, 'posterior', 'summary.csv'))
summary <- summary[summary$method == 'Exp', ]
point <- read.csv(file.path(run, 'point_estimate', 'summary.csv'))
stopifnot(all(summary$replicates == c(3, 3, 2)), all(summary$estimate_replicates == c(3, 2, 3)),
    all(abs(summary$average_estimate_absolute_bias - c(.01, .02, .03)) < 1e-12),
    all(abs(summary$average_absolute_bias - c(.04, .05, .06)) < 1e-12),
    is.na(p$S25_mean[1L]), is.finite(p$S25_estimate[1L]),
    is.na(p$S50_estimate[2L]), is.finite(p$S50_mean[2L]),
    point$replicates[point$method == 'Exp'] == 3L,
    point$replicates[point$method == 'Unif'] == 0L,
    point$replicates[point$method == 'NPMLE'] == 3L,
    point$replicates[point$method == 'Weibull_PT'] == 2L,
    point$replicates[point$method == 'Weibull_PT_Flat'] == 2L,
    point$replicates[point$method == 'OPT'] == 2L)
failures <- read.csv(file.path(run, 'fit_failures.csv'))
stopifnot(setequal(unique(failures$group), c('main', 'sensitivity', 'gamma_sensitivity')),
    any(failures$error_type == 'non_finite_point'), any(failures$status == 'partial_error'),
    all(c('S75_estimate', 'S75_estimate_absolute_bias', 'S75_mean', 'S75_absolute_bias') %in% names(p)))
paths <- list.files(run, '[.]csv$', recursive = TRUE, full.names = TRUE)
hashes <- tools::md5sum(paths); f$merge_run(run)
stopifnot(identical(hashes, tools::md5sum(paths)))
cat('PASS: protected PT/OPT errors, NA/Inf, validation errors, aliases, independent estimate/mean biases and denominators\n')
# Actual 0.4.2 tail NA: retain point curve and two finite posterior targets.
input <- e$generate_input(1L, e$all_settings()[1L, ])
input$quantiles <- c(input$tau / 4, input$tau / 2, input$tau + 1)
v <- e$fit_method('Unif', input, e$all_settings()[1L, ], 1L, config)
stopifnot(identical(v$status, 'partial_error'), identical(v$error_type, 'prediction_outside_tau'),
    all(is.finite(v$curve)), all(is.finite(v$posterior$estimate[1:2])),
    all(is.finite(v$posterior$mean[1:2])), is.na(v$posterior$estimate[3L]),
    is.na(v$posterior$mean[3L]), is.na(v$posterior$covered[3L]))
fit <- SurvPT::survPT(input$data[, c('L', 'R')], prior = 1L, flat = FALSE,
    kappa = .1, gamma = .5, max_M = config$max_M, max_EM_iter = config$max_EM_iter,
    EM_tol = config$EM_tol, verbose = FALSE)
raw <- predict(fit, input$quantiles, sample = FALSE)$estimate
stopifnot(identical(v$posterior$estimate, as.numeric(raw)))
cat('PASS: real 0.4.2 tail NA and exact MAP-EM estimates at posterior target times\n')
saveRDS(list(completed = Sys.time()), file.path(out, 'VERIFIED.rds'))
