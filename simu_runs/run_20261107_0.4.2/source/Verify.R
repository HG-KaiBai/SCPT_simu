#!/usr/bin/env Rscript
self <- sub('^--file=', '', grep('^--file=', commandArgs(), value = TRUE)[1L])
source(file.path(dirname(normalizePath(self)), 'Pipeline.R'))
args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args)) args[1L] else tempfile('empt-verification-')
if (file.exists(output)) stop('Verification output already exists: ', output)
dir.create(output, recursive = TRUE)
output <- normalizePath(output)
must_fail <- function(expr) stopifnot(inherits(tryCatch({force(expr); NULL}, error = identity), 'error'))
scientific <- function(z) {
    for (group in c('main', 'sensitivity', 'gamma_sensitivity')) for (method in names(z[[group]]))
        z[[group]][[method]]$seconds <- NULL
    z
}
settings <- all_settings()
methods <- pt_methods(TRUE, TRUE)
main_methods <- c('Unif', 'Unif_Flat', 'Exp', 'Exp_Flat', 'Weibull_PT', 'Weibull_PT_Flat')
stopifnot(identical(pt_methods()$method, main_methods), nrow(methods) == 39L,
    !anyDuplicated(methods$method), all(methods$prior %in% c(1L, 2L, 3L)),
    sum(methods$sensitivity) == 24L,
    scientific_config()$survpt_defaults$kappa == 0.1)
for (kappa in c(0.05, 0.1, 0.2, 0.3)) {
    spec <- methods[methods$sensitivity & methods$kappa == kappa, ]
    stopifnot(identical(spec$method, paste0(main_methods, '_kappa', kappa)),
        identical(spec$prior, c(1L, 1L, 2L, 2L, 3L, 3L)),
        identical(spec$flat, rep(c(FALSE, TRUE), 3L)))
}
seeds <- vapply(unique(canonical_method(methods$method)), function(method)
    method_seed(1L, 1L, method), integer(1))
stopifnot(!anyDuplicated(seeds))
for (method in main_methods) stopifnot(
    method_seed(1L, 1L, method) == method_seed(1L, 1L, paste0(method, '_kappa0.1')))
cat('PASS: six main prior/flat combinations, 24 kappa variants, deterministic distinct seeds\n')
stopifnot(nrow(settings) == 64L, sum(settings$suite == 'fiducial') == 28L,
          sum(settings$suite == 'exp') == 12L, sum(settings$suite == 'weibull') == 12L, sum(settings$suite == 'mixture') == 12L)
for (i in seq_len(nrow(settings))) {
    s <- settings[i, ]; input <- generate_input(1L, s)
    stopifnot(length(input$grid) == 201L, input$grid[1L] == 0,
        tail(input$grid, 1L) == input$tau, all(diff(input$grid) > 0),
        identical(input, generate_input(1L, s)), nrow(input$data) == s$n)
    truth <- if (s$suite == 'fiducial') pgamma(input$quantiles,
        shape = switch(s$type, '1' = 1, '2' = 3, '3' = 2, '4' = 1), lower.tail = FALSE) else if (s$suite == 'mixture') mixture_survival(input$quantiles) else
        pweibull(input$quantiles, s$K, 1 / s$Lambda, lower.tail = FALSE)
    stopifnot(isTRUE(all.equal(truth, c(.75, .5, .25))))
    stopifnot(('OPT' %in% expected_methods(s, 'all')$main) == (s$type == 'right'))
    stopifnot(length(expected_methods(s, 'all')$main) == 10L + as.integer(s$type == 'right'),
        length(expected_methods(s, 'all')$sensitivity) == 24L,
        length(expected_methods(s, 'all')$gamma_sensitivity) == 9L,
        length(expected_methods(s, 'main')$sensitivity) == 0L,
        length(expected_methods(s, 'sensitivity')$main) == 0L)
}
cat('PASS: all 64 designs, shared data, grids, true quantiles, OPT eligibility\n')

# Observe actual predict calls: the point grid must consume no random numbers.
input <- generate_input(1L, settings[1L, ])
fit <- SurvPT::survPT(input$data, max_M = 8L, max_EM_iter = 2000L, EM_tol = 1e-6, verbose = FALSE)
set.seed(42L); before <- .Random.seed
p <- predict(fit, input$grid, sample = FALSE)
stopifnot(identical(before, .Random.seed), identical(names(p), c('x', 'estimate')))
observed <- new.env(parent = environment(fit_method_strict))
observed$calls <- list()
observed$predict <- function(object, x, ...) {
    args <- list(...)
    observed$calls[[length(observed$calls) + 1L]] <- list(
        sample = if (is.null(args$sample)) TRUE else args$sample, count = length(x))
    stats::predict(object, x, ...)
}
instrumented <- fit_method_strict
environment(instrumented) <- observed
checked <- instrumented('Unif', input, settings[1L, ], 1L, scientific_config())
stopifnot(identical(observed$calls, list(list(sample = FALSE, count = 201L),
    list(sample = TRUE, count = 3L))), identical(checked$curve, as.numeric(p$estimate)),
    checked$sampling$n_samples == 2000L, checked$sampling$burn_in == 3000L,
    checked$sampling$thin == 1L)
cat('PASS: point predictions sample=FALSE; posterior only at 3 quantiles with defaults\n')

# Exercise the new prior/kappa interface for Weibull current and interval censoring.
for (id in c(34L, 38L)) {
    s <- settings[settings$setting_id == id, ]
    z <- simulate_replicate(1L, s, 'sensitivity', scientific_config())
    validate_result(z, 1L, s, list(analysis = 'sensitivity',
        config = scientific_config(), methods = methods))
}
cat('PASS: real default posterior sampling for all 24 variants under current/interval censoring\n')

# Real, unchanged-default fits: fiducial and right-censored exponential, two seeds each.
run <- file.path(output, 'run')
prepare_run(run, replicates = 2L, chunk_size = 2L,
            settings = settings[settings$setting_id %in% c(1L, 29L), ], analysis = 'all')
frozen <- new.env(parent = globalenv())
source(file.path(run, 'source', 'Pipeline.R'), local = frozen)
must_fail(frozen$merge_run(run))
stopifnot(!dir.exists(file.path(run, 'point_estimate')))
for (i in seq_len(nrow(frozen$read_run(run)$tasks))) frozen$run_task(run, i, cores = 2L)
m <- frozen$read_run(run)
for (i in seq_len(nrow(m$tasks))) {
    task <- m$tasks[i, , drop = FALSE]
    chunk <- frozen$read_chunk(run, m, task)
    s <- m$settings[match(task$setting_id, m$settings$setting_id), ]
    for (z in chunk$results) {
        for (method in main_methods) stopifnot(
            identical(z$main[[method]], z$sensitivity[[paste0(method, '_kappa0.1')]]))
        for (method in c('Unif', 'Exp', 'Weibull_PT')) stopifnot(
            identical(z$main[[method]], z$gamma_sensitivity[[paste0(method, '_gamma0.5')]]))
    }
    # Serial main/sensitivity-only runs must match the combined parallel run.
    if (s$setting_id == 1L) {
        for (analysis in c('main', 'sensitivity', 'gamma')) {
            serial <- frozen$simulate_replicate(1L, s, analysis, m$config)
            frozen$validate_result(serial, 1L, s, modifyList(m, list(analysis = analysis)))
            group <- if (analysis == 'gamma') 'gamma_sensitivity' else analysis
            stopifnot(identical(scientific(serial)[[group]],
                scientific(chunk$results[['1']])[[group]]),
                identical(serial$input, chunk$results[['1']]$input))
        }
    }
    path <- frozen$chunk_path(run, i)
    hash <- unname(tools::md5sum(path))
    frozen$run_task(run, i, cores = 1L)
    stopifnot(identical(hash, unname(tools::md5sum(path))))
}
cat('PASS: all methods/default draws, serial-parallel and analysis-mode reproducibility, kappa=0.1 reuse, resume\n')

path <- frozen$chunk_path(run, 1L)
original <- readRDS(path)
bad <- original; bad$manifest_id <- 'wrong-run'
saveRDS(bad, path)
must_fail(frozen$merge_run(run))
saveRDS(original, path)
frozen$merge_run(run)
stopifnot(file.exists(file.path(run, 'COMPLETE.rds')))
for (group in c('main', 'sensitivity', 'gamma_sensitivity')) {
    root <- group_root(run, group)
    summary <- read.csv(file.path(root, 'point_estimate', 'summary.csv'))
    stopifnot(nrow(summary) == switch(group, main = 21L, sensitivity = 48L, gamma_sensitivity = 18L))
    for (i in seq_len(nrow(summary))) {
        row <- summary[i, ]
        data <- read.csv(file.path(root, 'point_estimate', row$setting, 'replicates.csv'))
        data <- data[data$method == row$method, ]
        for (metric in c('SSE', 'RMSE', 'RMSE_nolast')) {
            stopifnot(isTRUE(all.equal(row[[paste0(metric, '_mean')]], mean_available(data[[metric]]))),
                      isTRUE(all.equal(row[[paste0(metric, '_sd')]], sd_available(data[[metric]]))))
        }
    }
    summary <- read.csv(file.path(root, 'posterior', 'summary.csv'))
    stopifnot(nrow(summary) == switch(group, main = 54L, sensitivity = 144L, gamma_sensitivity = 54L))
    for (i in seq_len(nrow(summary))) {
        row <- summary[i, ]
        tab <- read.csv(file.path(root, 'posterior', row$setting, paste0(row$method, '.csv')))
        stopifnot(identical(tail(tab$replicate, 1L), 'Average'), nrow(tab) == 3L)
        prefix <- paste0('S', c('75', '50', '25')[match(row$probability, c(.75, .5, .25))], '_')
        p <- tab[tab$replicate != 'Average', ]
        stopifnot(isTRUE(all.equal(row$average_absolute_bias, mean_available(p[[paste0(prefix, 'absolute_bias')]]))),
            isTRUE(all.equal(row$average_length, mean_available(p[[paste0(prefix, 'length')]]))),
            isTRUE(all.equal(row$coverage_rate, mean_available(p[[paste0(prefix, 'covered')]]))))
        lower <- p[[paste0(prefix, 'lower')]]; upper <- p[[paste0(prefix, 'upper')]]
        expected_score <- upper - lower + 40 *
            (ifelse(row$probability < lower, lower - row$probability, 0) +
             ifelse(row$probability > upper, row$probability - upper, 0))
        expected_absolute_bias <- abs(p[[paste0(prefix, 'mean')]] - row$probability)
        stopifnot(isTRUE(all.equal(p[[paste0(prefix, 'absolute_bias')]], expected_absolute_bias)),
            isTRUE(all.equal(row$average_absolute_bias, mean_available(expected_absolute_bias))))
        score_column <- paste0(prefix, 'winkler_score')
        stopifnot(isTRUE(all.equal(p[[score_column]], expected_score)),
            isTRUE(all.equal(tab[tab$replicate == 'Average', score_column], mean_available(expected_score))),
            isTRUE(all.equal(row$average_winkler_score, mean_available(expected_score))))
    }
}
cat('PASS: incomplete/corrupt runs rejected; merged means, replication-only SD, posterior coverage and Winkler scores\n')
atomic_write(list(completed = Sys.time(), session = sessionInfo(),
    source_md5 = tools::md5sum(file.path(.simu_dir, .sources))), file.path(output, 'VERIFIED.rds'))
cat('ALL CHECKS PASSED:', output, '\n')
