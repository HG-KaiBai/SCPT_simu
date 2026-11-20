# Isolate SurvPT (including sensitivity variants and OPT) and MPT/TBP errors.
method_error_policy <- list(version = 2L,
    methods = c(pt_methods(TRUE, TRUE)$method, 'OPT', 'MPT', 'TBP'),
    action = 'retain NA failures; preserve valid PT point/quantile results; use metric-specific denominators')
method_failure_files <- c('fit_failures.csv', 'fit_failure_summary.csv', 'fit_status_summary.csv')
is_skipped_fit <- function(x) identical(x$status, 'skipped_error')
has_method_issue <- function(x) is_skipped_fit(x) || identical(x$status, 'partial_error')
failure_reason <- function(message) {
    if (grepl('outside [0,1] by more than tolerance', message, fixed = TRUE)) return('probability_out_of_bounds')
    if (grepl('survreg initializer', message, fixed = TRUE)) return('initializer_failure')
    if (grepl('singular', message, ignore.case = TRUE)) return('singular_matrix')
    if (grepl('Non-finite point', message, fixed = TRUE)) return('non_finite_point')
    if (grepl('Non-finite posterior', message, fixed = TRUE)) return('non_finite_posterior')
    if (grepl('all(c(mean', message, fixed = TRUE)) return('posterior_out_of_bounds')
    'other_error'
}
missing_posterior <- function(input) {
    data.frame(probability = input$probabilities, time = input$quantiles,
        mean = NA_real_, lower = NA_real_, upper = NA_real_, absolute_bias = NA_real_,
        length = NA_real_, covered = NA_integer_, winkler_score = NA_real_,
        estimate = NA_real_, estimate_absolute_bias = NA_real_)
}
posterior_failure_info <- function(method, input, posterior) {
    missing <- is.na(posterior$mean) | is.na(posterior$estimate)
    stopifnot(any(missing))
    outside <- all(input$quantiles[missing] > input$tau)
    list(status = 'partial_error', method = canonical_method(method),
        error = if (outside) 'Posterior target beyond tau: SurvPT returns NA' else
            'Non-finite posterior prediction (NA/NaN/Inf)',
        error_type = if (outside) 'prediction_outside_tau' else 'non_finite_posterior',
        error_class = c('survpt_prediction_missing', 'error', 'condition'),
        error_call = 'predict.survPT', error_stage = 'posterior',
        affected_probabilities = paste(input$probabilities[missing], collapse = ';'))
}
skipped_method_fit <- function(method, input, s, t, err, seconds) {
    stopifnot(method %in% method_error_policy$methods)
    list(status = 'skipped_error', method = canonical_method(method), error = conditionMessage(err),
        error_type = failure_reason(conditionMessage(err)),
        error_class = class(err), error_call = paste(deparse(conditionCall(err)), collapse = ' '),
        error_stage = 'fit_prediction_or_validation',
        affected_probabilities = if (method == 'OPT') '' else paste(input$probabilities, collapse = ';'),
        curve = rep(NA_real_, length(input$grid)),
        errors = setNames(rep(NA_real_, 3L), c('SSE', 'RMSE', 'RMSE_nolast')),
        posterior = if (method == 'OPT') NULL else missing_posterior(input),
        sampling = NULL, prior = NULL, convergence = NULL,
        seed = method_seed(t, s$setting_id, method), seconds = seconds, warnings = character())
}
validate_skipped_fit <- function(v, name, group = NULL, input, s, t) {
    actual <- canonical_method(name)
    stopifnot(is_skipped_fit(v), actual %in% method_error_policy$methods,
        identical(v$method, actual), is.character(v$error), length(v$error) == 1L, nzchar(v$error),
        identical(v$error_type, failure_reason(v$error)), is.character(v$error_class),
        'error' %in% v$error_class, is.character(v$error_call), length(v$error_call) == 1L,
        identical(v$error_stage, 'fit_prediction_or_validation'),
        identical(v$affected_probabilities, if (actual == 'OPT') '' else paste(input$probabilities, collapse = ';')),
        identical(v$curve, rep(NA_real_, length(input$grid))),
        identical(v$errors, setNames(rep(NA_real_, 3L), c('SSE', 'RMSE', 'RMSE_nolast'))),
        identical(v$posterior, if (actual == 'OPT') NULL else missing_posterior(input)),
        identical(v$seed, method_seed(t, s$setting_id, actual)),
        length(v$seconds) == 1L, is.finite(v$seconds), v$seconds >= 0)
    invisible(TRUE)
}
mean_available <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
sd_available <- function(x) if (sum(!is.na(x)) < 2L) NA_real_ else sd(x, na.rm = TRUE)
empty_failures <- function() data.frame(setting_id = integer(), setting = character(),
    replicate = integer(), group = character(), method = character(), canonical_method = character(),
    seed = integer(), status = character(), error_stage = character(), affected_probabilities = character(),
    tau = numeric(), error_type = character(), error_message = character(), error_class = character(),
    error_call = character(), warnings = character(), stringsAsFactors = FALSE)
collect_method_failures <- function(results, s) {
    rows <- list()
    for (z in results) for (group in c('main', 'sensitivity', 'gamma_sensitivity'))
        for (name in names(z[[group]])) {
            v <- z[[group]][[name]]
            if (has_method_issue(v)) rows[[length(rows) + 1L]] <- data.frame(
                setting_id = s$setting_id, setting = s$label, replicate = z$replicate,
                group = group, method = name, canonical_method = canonical_method(name), seed = v$seed,
                status = v$status, error_stage = v$error_stage, affected_probabilities = v$affected_probabilities,
                tau = z$input$tau, error_type = v$error_type, error_message = v$error,
                error_class = paste(v$error_class, collapse = ';'), error_call = v$error_call,
                warnings = paste(v$warnings, collapse = ';'), stringsAsFactors = FALSE)
        }
    if (length(rows)) do.call(rbind, rows) else empty_failures()
}
collect_method_status <- function(results, s, m) {
    rows <- list(); expected <- expected_methods(s, m$analysis)
    for (group in names(expected)) for (name in expected[[group]]) {
        values <- lapply(results, function(z) z[[group]][[name]])
        stopifnot(!any(vapply(values, is.null, logical(1))))
        failed <- sum(vapply(values, has_method_issue, logical(1)))
        skipped <- sum(vapply(values, is_skipped_fit, logical(1)))
        row <- data.frame(setting_id = s$setting_id, setting = s$label,
            group = group, method = name, expected = m$replicates, attempted = length(values),
            successful = length(values) - failed, failed = failed, skipped = skipped,
            partial = failed - skipped, point_successful = length(values) - skipped,
            unattempted = m$replicates - length(values))
        for (j in 1:3) row[[paste0('posterior_S', c('75', '50', '25')[j], '_successful')]] <-
            if (name %in% c('OPT', 'NPMLE')) NA_integer_ else
                sum(vapply(values, function(v) !is.na(v$posterior$mean[j]), logical(1)))
        for (j in 1:3) row[[paste0('estimate_S', c('75', '50', '25')[j], '_successful')]] <-
            if (name %in% c('OPT', 'NPMLE')) NA_integer_ else
                sum(vapply(values, function(v) !is.na(v$posterior$estimate[j]), logical(1)))
        rows[[length(rows) + 1L]] <- row
    }
    do.call(rbind, rows)
}
write_method_failure_reports <- function(failures, status, root) {
    stopifnot(!anyDuplicated(failures[c('setting_id', 'replicate', 'group', 'method')]),
        !anyDuplicated(status[c('setting_id', 'group', 'method')]), sum(status$failed) == nrow(failures))
    summary <- data.frame(method = character(), error_type = character(), failed_fits = integer())
    if (nrow(failures)) {
        summary <- aggregate(list(failed_fits = rep(1L, nrow(failures))),
            failures[c('method', 'error_type')], sum)
        summary <- summary[order(summary$method, summary$error_type), ]
        failures <- failures[order(failures$setting_id, failures$replicate, failures$method), ]
    }
    atomic_write(failures, file.path(root, method_failure_files[1L]), csv = TRUE)
    atomic_write(summary, file.path(root, method_failure_files[2L]), csv = TRUE)
    atomic_write(status, file.path(root, method_failure_files[3L]), csv = TRUE)
    invisible(summary)
}
