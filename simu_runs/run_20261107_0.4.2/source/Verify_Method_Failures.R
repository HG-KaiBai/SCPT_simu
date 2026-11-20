#!/usr/bin/env Rscript
self <- sub('^--file=', '', grep('^--file=', commandArgs(), value = TRUE)[1L])
e <- new.env(parent = globalenv())
source(file.path(dirname(normalizePath(self)), 'Pipeline.R'), local = e)
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[1L] else tempfile('verify-method-errors-')
if (file.exists(out)) stop('Use a new verification directory: ', out)
dir.create(out, recursive = TRUE)
must_fail <- function(expr) stopifnot(inherits(tryCatch({force(expr); NULL}, error = identity), 'error'))
config <- e$scientific_config()
s <- e$all_settings()[56L, ]; input <- e$generate_input(701L, s)
a <- e$fit_method_strict('NPMLE', input, s, 701L, config)
b <- e$fit_method('NPMLE', input, s, 701L, config)
a$seconds <- b$seconds <- NULL
stopifnot(identical(a, b))
cat('PASS: successful NPMLE fits unchanged\n')

# Verify the default complete design can be frozen without historical tooling.
full <- file.path(out, 'prepared_full')
e$prepare_run(full)
f <- new.env(parent = globalenv()); source(file.path(full, 'source', 'Pipeline.R'), local = f)
m <- f$read_run(full)
stopifnot(nrow(m$settings) == 64L, m$replicates == 1000L, nrow(m$tasks) == 640L,
    nrow(m$methods) == 39L, identical(m$method_error_policy, e$method_error_policy),
    setequal(list.files(file.path(full, 'source')), e$.sources),
    !any(grepl('Supplement|Resume_Method|Backfill|Mixture_Gamma', e$.sources)))
cat('PASS: complete 64 x 1000 design is standalone and includes the error policy\n')

# Inject known failures to exercise the complete worker/merge pipeline, including
# a method for which every replicate fails. Other fitted metadata stay valid.
install_test_fits <- function(env) {
    f <- function(method, input, s, t, config) {
        if (s$setting_id == 56L && (method == 'TBP' || (method == 'MPT' && t == 1L)))
            stop(if (method == 'MPT') 'injected MPT failure' else 'injected TBP failure')
        estimate <- pmax(0, pmin(1, input$probabilities + c(0, .02, -.03)[t]))
        p <- if (method %in% c('NPMLE', 'OPT')) NULL else
            posterior_table(input, estimate, pmax(0, estimate - .025), pmin(1, estimate + .025), config$level)
        spec <- pt_methods(TRUE, TRUE); spec <- spec[spec$method == method, , drop = FALSE]
        prior <- sampling <- NULL
        if (nrow(spec)) {
            prior <- list(family = c('uniform', 'exponential', 'weibull')[spec$prior],
                flat = spec$flat, kappa = spec$kappa, gamma = spec$gamma,
                parameters = if (spec$prior == 3L) c(.5, 2) else 1)
            sampling <- config$survpt_sampling
        }
        if (method == 'Fiducial') sampling <- config$fiducial
        list(curve = input$truth, errors = both_errors(input$truth, input$truth), posterior = p,
            sampling = sampling, prior = prior, convergence = NULL,
            seed = method_seed(t, s$setting_id, method), seconds = 0, warnings = character())
    }
    environment(f) <- env; env$fit_method_strict <- f
}
scientific <- function(z) {
    for (group in c('main', 'sensitivity', 'gamma_sensitivity'))
        for (method in names(z[[group]])) z[[group]][[method]]$seconds <- NULL
    z
}
outputs <- list()
for (partition in c(1L, 3L)) {
    run <- file.path(out, paste0('partition_', partition))
    e$prepare_run(run, replicates = 3L, chunk_size = partition, settings = s, analysis = 'all')
    frozen <- new.env(parent = globalenv())
    source(file.path(run, 'source', 'Pipeline.R'), local = frozen)
    install_test_fits(frozen); m <- frozen$read_run(run)
    must_fail(frozen$merge_run(run))
    for (i in seq_len(nrow(m$tasks))) frozen$run_task(run, i, cores = if (partition == 1L) 1L else 2L)
    paths <- vapply(m$tasks$task_id, function(i) frozen$chunk_path(run, i), character(1))
    hashes <- tools::md5sum(paths)
    for (i in seq_len(nrow(m$tasks))) frozen$run_task(run, i, 1L)
    stopifnot(identical(hashes, tools::md5sum(paths)))
    frozen$merge_run(run)
    label <- s$label
    point <- read.csv(file.path(run, 'point_estimate', label, 'replicates.csv'))
    ps <- read.csv(file.path(run, 'point_estimate', label, 'summary.csv'))
    post <- read.csv(file.path(run, 'posterior', label, 'MPT.csv'))
    tbp <- read.csv(file.path(run, 'posterior', label, 'TBP.csv'))
    failure <- read.csv(file.path(run, 'fit_failures.csv'))
    counts <- read.csv(file.path(run, 'fit_failure_summary.csv'))
    status <- read.csv(file.path(run, 'fit_status_summary.csv'))
    stopifnot(identical(point$replicate[point$method == 'MPT'], 1:3),
        all(is.na(point[point$method == 'MPT' & point$replicate == 1L, c('SSE', 'RMSE', 'RMSE_nolast')])),
        ps$replicates[ps$method == 'MPT'] == 2L, ps$replicates[ps$method == 'TBP'] == 0L,
        ps$replicates[ps$method == 'NPMLE'] == 3L,
        all(is.na(ps[ps$method == 'TBP', grep('_(mean|sd)$', names(ps))])),
        identical(post$replicate, c('1', '2', '3', 'Average')),
        is.na(post$S75_mean[1L]), abs(post$S75_mean[4L] - .745) < 1e-14,
        post$S75_covered[4L] == .5, all(is.na(tbp$S75_mean)), all(is.na(tbp$S75_covered)),
        nrow(failure) == 4L, sum(counts$failed_fits) == 4L,
        counts$failed_fits[counts$method == 'MPT'] == 1L,
        counts$failed_fits[counts$method == 'TBP'] == 3L, all(status$unattempted == 0L))
    csvs <- list.files(run, pattern = '[.]csv$', recursive = TRUE, full.names = TRUE)
    hash <- tools::md5sum(csvs); frozen$merge_run(run)
    stopifnot(identical(hash, tools::md5sum(csvs)))
    results <- unlist(lapply(seq_len(nrow(m$tasks)), function(i)
        frozen$read_chunk(run, m, m$tasks[i, , drop = FALSE])$results), recursive = FALSE)
    outputs[[as.character(partition)]] <- lapply(results[as.character(1:3)], scientific)
    bad <- results[['1']]; bad$main$MPT$curve[1L] <- 0
    must_fail(frozen$validate_result(bad, 1L, s, m))
    bad <- results[['1']]; bad$main$NPMLE <- bad$main$MPT
    must_fail(frozen$validate_result(bad, 1L, s, m))
    frozen$fit_method_strict <- function(...) stop('unrelated method failure')
    must_fail(frozen$fit_method('NPMLE', input, s, 701L, config))
    bad <- readRDS(file.path(run, 'manifest.rds')); bad$method_error_policy$version <- -1L
    saveRDS(bad, file.path(run, 'manifest.rds')); must_fail(frozen$read_run(run))
    saveRDS(m[names(m) != 'id'], file.path(run, 'manifest.rds'))
}
stopifnot(identical(outputs[['1']], outputs[['3']]))
cat('PASS: method isolation, NA exports, valid denominators, all-failed NA, resume, repeatable merge\n')
cat('PASS: serial/parallel and different chunk sizes give identical scientific results\n')
cat('PASS: corrupted records, unsupported skips and mismatched error policy rejected\n')
saveRDS(list(completed = Sys.time(), output = normalizePath(out), checks = 'passed'), file.path(out, 'VERIFIED.rds'))
cat('ALL METHOD ERROR HANDLING CHECKS PASSED\n')
