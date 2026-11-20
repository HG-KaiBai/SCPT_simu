#!/usr/bin/env Rscript
self <- sub('^--file=', '', grep('^--file=', commandArgs(), value = TRUE)[1L])
source(file.path(dirname(normalizePath(self)), 'Time_Evaluation.R'))
must_fail <- function(expr) stopifnot(inherits(tryCatch({force(expr); NULL}, error = identity), 'error'))
cfg <- time_config(); settings <- time_settings()
stopifnot(nrow(settings) == 18L, !anyDuplicated(settings$setting_id),
    setequal(settings$n, c(200L, 500L, 1000L)), setequal(settings$q, c(.3, .5)),
    setequal(settings$censoring, c('RC', 'CS', 'IC')), all(settings$K == 2),
    cfg$pt_sampling$burn_in + cfg$pt_sampling$n_samples * cfg$pt_sampling$thin == 5000L,
    cfg$mcmc$nburn + cfg$mcmc$nsave * (cfg$mcmc$nskip + 1) == 5000L,
    cfg$fiducial$mburn == 100, cfg$fiducial$mfid == 1000, identical(cfg$InitParamMCMC, TRUE),
    identical(cfg$mcmc, eval(formals(spBayesSurv::survregbayes2)$mcmc)),
    identical(cfg$mcmc, eval(formals(spBayesSurv::survregbayes)$mcmc)),
    isTRUE(eval(formals(spBayesSurv::survregbayes2)$InitParamMCMC)),
    isTRUE(eval(formals(spBayesSurv::survregbayes)$InitParamMCMC)))
for (i in seq_len(nrow(settings))) {
    input <- generate_input(1L, settings[i, ])
    stopifnot(length(input$data$L) == settings$n[i], length(input$grid) == 201L,
        input$grid[1L] == 0, input$grid[201L] == input$tau,
        identical(input, generate_input(1L, settings[i, ])))
}
bad <- time_stage(function() stop('injected error'))
stopifnot(bad$status == 'error', bad$error == 'injected error', bad$elapsed_seconds >= 0)
bad_validation <- time_stage(function() 1, function(x) stop('invalid prediction'))
stopifnot(bad_validation$status == 'error', bad_validation$error == 'invalid prediction')
warn <- time_stage(function() { warning('record me'); 5 })
stopifnot(warn$status == 'ok', warn$warnings == 'record me', warn$value == 5)
# Quantiles exclude failed fits and include explicit empty groups.
raw <- data.frame(setting_id = settings$setting_id[1L], method = 'MPT', stage = 'total',
    status = c(rep('ok', 4), 'error'), elapsed_seconds = c(1, 2, 3, 4, 1000))
m <- list(settings = settings, replicates = 50L, config = cfg)
summary <- time_summary(raw, m)
z <- summary[summary$setting_id == settings$setting_id[1L] & summary$method == 'MPT', ]
stopifnot(nrow(summary) == 216L, z$median_seconds == 2.5, z$q1_seconds == 1.75,
    z$q3_seconds == 3.25, z$iqr_seconds == 1.5, z$successful == 4L, z$failed == 1L,
    all(is.na(summary$median_seconds[summary$attempted == 0])))
must_fail(assert_time_runtime(list(job_id = '', cpu_model = 'test', affinity = '0-1', threads = '1')))
# Synthetic report/checkpoint checks are isolated from every real run directory.
check_reports <- function() {
    root <- tempfile('empt-time-report-'); dir.create(root)
    on.exit(unlink(root, recursive = TRUE), add = TRUE)
    s <- settings[1L, ]; input <- generate_input(1L, s)
    task <- data.frame(replicate = 1L, setting_id = s$setting_id, task_id = 1L)
    m <- list(id = 'test-only', settings = s, tasks = task, replicates = 1L,
        config = cfg, R_version = as.character(getRversion()), packages = package_versions())
    raw <- do.call(rbind, lapply(cfg$methods, function(method) {
        stages <- if (method %in% c('Unif', 'Exp', 'Weibull_PT')) c('MAP-EM', 'DA-Gibbs', 'total') else 'total'
        data.frame(setting_id = s$setting_id, method = method, stage = stages,
            status = 'ok', elapsed_seconds = seq_along(stages), grid_points = 201L,
            host = 'synthetic-host', cpu_model = 'synthetic-test-data')
    }))
    z <- list(manifest_id = m$id, task = task, input = input,
        data_hash = digest::digest(input, algo = 'sha256', serializeVersion = 3), rows = raw)
    path <- time_checkpoint_path(root, 1L)
    atomic_write(z, path)
    stopifnot(identical(read_time_checkpoint(root, task, m), z))
    report_time(root, m, complete = TRUE)
    first <- read.csv(file.path(root, 'summary.csv'))
    report_time(root, m, complete = TRUE)
    stopifnot(identical(first, read.csv(file.path(root, 'summary.csv'))), nrow(first) == 12L,
        file.exists(file.path(root, 'COMPLETE.rds')), all(first$successful == 1L))
    z$input$tau <- z$input$tau + 1
    atomic_write(z, path)
    must_fail(read_time_checkpoint(root, task, m))
}
check_reports()
cat('PASS: design, deterministic datasets, sweep budgets, failures, quartiles, checkpoint integrity and repeatable reports\n')
if ('--real' %in% commandArgs(trailingOnly = TRUE)) {
    runtime <- time_runtime()
    # Validation may use a shared compute node; measured production work may not.
    stopifnot(nzchar(runtime$job_id), grepl('^[0-9]+$', runtime$affinity), runtime$threads == '1')
    # Real package calls across all censoring types, with the full required budgets.
    for (type in c('right', 'current', 'interval')) {
        s <- settings[settings$type == type & settings$n == 200L & settings$q == .5, ]
        s$n <- 20L
        input <- generate_input(123L, s)
        values <- lapply(cfg$methods, function(method) time_one_method(method, input, s, 123L))
        rows <- do.call(rbind, lapply(values, `[[`, 'rows'))
        print(rows[, c('method', 'stage', 'status', 'elapsed_seconds', 'error')], row.names = FALSE)
        stopifnot(nrow(rows) == 12L, all(rows$status == 'ok'), all(rows$grid_points == 201L))
        for (method in c('Unif', 'Exp', 'Weibull_PT')) {
            z <- rows[rows$method == method, ]
            stopifnot(z$elapsed_seconds[z$stage == 'total'] ==
                sum(z$elapsed_seconds[z$stage != 'total']))
        }
    }
    cat('PASS: real serial fits, common grids, stage totals and actual package sampling budgets\n')
}
if ('--large' %in% commandArgs(trailingOnly = TRUE)) {
    stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
    s <- settings[settings$type == 'interval' & settings$n == 1000L & settings$q == .5, ]
    input <- generate_input(1L, s)
    rows <- do.call(rbind, lapply(cfg$methods, function(method)
        time_one_method(method, input, s, 1L)$rows))
    print(rows[, c('method', 'stage', 'status', 'elapsed_seconds', 'error')], row.names = FALSE)
    stopifnot(all(rows$status == 'ok'))
    cat('PASS: n=1000 IC resource-sizing check (not a benchmark result)\n')
}
