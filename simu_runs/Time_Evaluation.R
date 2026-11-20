#!/usr/bin/env Rscript
# A separate, frozen, serial benchmark; never reuses multi-core simulation times.
.time_sources <- Filter(Negate(is.null), lapply(sys.frames(), function(f) f$ofile))
.time_self <- if (length(.time_sources)) tail(.time_sources, 1L)[[1L]] else
    sub('^--file=', '', grep('^--file=', commandArgs(), value = TRUE)[1L])
.time_dir <- dirname(normalizePath(.time_self, mustWork = TRUE))
source(file.path(.time_dir, 'Pipeline.R'), local = TRUE)
.time_files <- c('Time_Evaluation.R', 'Time_Run.sbatch', 'Time_Submit.sh', 'Verify_Time_Evaluation.R')

time_config <- function() list(methods = c('Unif', 'Exp', 'Weibull_PT', 'MPT', 'TBP', 'Fiducial'),
    n = c(200L, 500L, 1000L), q = c(.3, .5), shape = 2, inverse_scale = .12,
    grid_points = 201L, max_M = 8L, max_EM_iter = 2000L, EM_tol = 1e-6,
    flat = FALSE, kappa = .1, gamma = .5, level = .95,
    pt_sampling = list(n_samples = 2000L, burn_in = 3000L, thin = 1L),
    mcmc = list(nburn = 3000, nsave = 2000, nskip = 0, ndisplay = 500),
    InitParamMCMC = TRUE,
    fiducial = list(mfid = 1000, mburn = 100), order_seed = 20261108L,
    scope = 'elapsed seconds; fitting plus prediction at seq(0,tau,length.out=201); excludes data generation, post-fit validation, explicit pre-fit GC and benchmark file writes',
    pt_stages = 'MAP-EM includes centering-parameter estimation, tree construction, fitting and sample=FALSE prediction; DA-Gibbs is predict(sample=TRUE), including sampling and grid summaries; total is their per-dataset sum',
    other_stages = 'MPT/TBP call the original simulation fit_mpt()/fit_tbp() including all parameter fits and default InitParamMCMC=TRUE: the extra 5000+5000 parametric initialization chain is included in total time, in addition to the 3000+2000 main chain and grid summaries; Fiducial includes its original sampler, interpolation and summaries',
    quantile_type = 7L)

time_settings <- function() {
    base <- all_settings()
    rows <- list()
    next_id <- max(base$setting_id)
    for (n in time_config()$n) for (type in c('right', 'current', 'interval'))
        for (q in time_config()$q) {
            s <- base[base$suite == 'weibull' & base$n == 200L & base$type == type & base$q == q, ]
            old <- base[base$suite == 'weibull' & base$n == n & base$type == type & base$q == q, ]
            if (nrow(old)) s <- old else {
                next_id <- next_id + 1L
                s$n <- n; s$setting_id <- next_id
                s$label <- sub('_n200_', paste0('_n', n, '_'), s$label, fixed = TRUE)
            }
            s$censoring <- c(right = 'RC', current = 'CS', interval = 'IC')[[type]]
            rows[[length(rows) + 1L]] <- s
        }
    out <- do.call(rbind, rows); rownames(out) <- NULL
    out
}

prepare_time <- function(parent, replicates = 50L) {
    parent <- normalizePath(parent, mustWork = TRUE)
    replicates <- positive_integer(replicates, 'replicates')
    stopifnot(replicates <= 9999L)
    original <- readRDS(file.path(parent, 'manifest.rds'))
    stopifnot(identical(original$R_version, as.character(getRversion())),
        identical(original$packages, package_versions()),
        identical(original$archive_md5, installed_survpt_archive()),
        identical(unname(original$source_md5), unname(tools::md5sum(
            file.path(parent, 'source', names(original$source_md5))))))
    root <- file.path(parent, 'time_evaluation')
    if (file.exists(root)) stop('Timing directory already exists: ', root)
    dir.create(root)
    for (d in c('source', 'checkpoints', 'logs', 'sessions')) dir.create(file.path(root, d))
    stopifnot(all(file.copy(file.path(parent, 'source', names(original$source_md5)), file.path(root, 'source'))),
        all(file.copy(file.path(.time_dir, .time_files), file.path(root, 'source'))))
    settings <- time_settings()
    tasks <- expand.grid(replicate = seq_len(replicates), setting_id = settings$setting_id)
    set.seed(time_config()$order_seed)
    tasks <- tasks[sample.int(nrow(tasks)), ]; rownames(tasks) <- NULL
    tasks$task_id <- seq_len(nrow(tasks))
    files <- sort(unique(c(names(original$source_md5), .time_files)))
    hashes <- tools::md5sum(file.path(root, 'source', files)); names(hashes) <- files
    m <- list(version = 1L, parent = parent,
        parent_manifest_md5 = unname(tools::md5sum(file.path(parent, 'manifest.rds'))),
        settings = settings, tasks = tasks, replicates = replicates, config = time_config(),
        R_version = as.character(getRversion()), packages = package_versions(), RNGkind = RNGkind(),
        archive_md5 = installed_survpt_archive(), source_md5 = hashes, created = Sys.time())
    atomic_write(m, file.path(root, 'manifest.rds'))
    atomic_write(settings, file.path(root, 'settings.csv'), csv = TRUE)
    atomic_write(tasks, file.path(root, 'tasks.csv'), csv = TRUE)
    writeLines(c('PREPARED', paste(nrow(settings), 'settings x', replicates, 'datasets'),
        'No measured results yet. Submit the frozen Time_Run.sbatch on an exclusive compute node.'),
        file.path(root, 'STATUS.txt'))
    writeLines(c('# Timing protocol', '',
        sprintf('18 combinations x %d datasets: Weibull(shape=2, inverse scale=0.12), RC/CS/IC, n=200/500/1000, q=0.3/0.5.', replicates),
        'Methods: Unif, Exp, Weibull_PT (all non-flat), MPT, TBP, Fiducial.',
        'Every evaluation uses 201 equally spaced nodes from 0 to dataset tau.',
        'Bayesian budgets: 3000 burn-in + 2000 retained, thin=1 (spBayesSurv nskip=0).',
        'Fiducial: original 100 burn-in + 1000 retained samples.', '',
        time_config()$scope, time_config()$pt_stages, time_config()$other_stages, '',
        'Execution: exclusive Slurm node, one CPU affinity, all numerical library threads=1.',
        'The exclusive allocation reserves the whole node although only one CPU does benchmark work.',
        'Warm-up is excluded. Data generation and explicit pre-fit GC are excluded; GC triggered during a fit remains included.',
        'Median, Q1, Q3 and IQR (Q3-Q1), in seconds, use R quantile(type=7). Failures are counted separately.',
        'Source and parent-run hashes are frozen in manifest.rds. Full measured environment is saved when the worker starts.',
        'See STATUS.txt and job_id.txt for progress; report.md and summary.csv appear after the first completed dataset.',
        'Only COMPLETE.rds certifies all datasets have been attempted. Partial summaries are explicitly labelled.',
        'Recovery must stay on the original host to avoid pooling timings across nodes.'), file.path(root, 'protocol.md'))
    message('Prepared timing benchmark: ', root)
    invisible(m)
}

read_time <- function(root, environment = TRUE) {
    m <- readRDS(file.path(root, 'manifest.rds'))
    stopifnot(identical(m$version, 1L), identical(m$config, time_config()),
        identical(unname(m$source_md5), unname(tools::md5sum(file.path(root, 'source', names(m$source_md5))))),
        identical(m$parent_manifest_md5, unname(tools::md5sum(file.path(m$parent, 'manifest.rds')))))
    if (environment) stopifnot(identical(m$R_version, as.character(getRversion())),
        identical(m$packages, package_versions()), identical(m$archive_md5, installed_survpt_archive()),
        identical(unname(m$source_md5), unname(tools::md5sum(file.path(.time_dir, names(m$source_md5))))))
    m$id <- unname(tools::md5sum(file.path(root, 'manifest.rds')))
    m
}

proc_field <- function(key) {
    lines <- readLines('/proc/self/status')
    trimws(sub('^[^:]+:', '', lines[startsWith(lines, paste0(key, ':'))]))
}
time_runtime <- function() {
    cpu <- readLines('/proc/cpuinfo')
    list(host = Sys.info()[['nodename']], cpu_model = unique(trimws(sub('^[^:]+:', '',
        cpu[grepl('^model name[[:space:]]*:', cpu)]))),
        affinity = proc_field('Cpus_allowed_list'), threads = proc_field('Threads'),
        job_id = Sys.getenv('SLURM_JOB_ID'),
        thread_environment = Sys.getenv(c('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
            'MKL_NUM_THREADS', 'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'RCPP_PARALLEL_NUM_THREADS')),
        R = R.version.string, packages = package_versions(), session = sessionInfo(),
        load = readLines('/proc/loadavg'), started = Sys.time())
}
assert_time_runtime <- function(runtime) {
    stopifnot(nzchar(runtime$job_id), Sys.getenv('EMPT_TIME_EXCLUSIVE') == '1',
        length(runtime$cpu_model) == 1L, grepl('^[0-9]+$', runtime$affinity),
        runtime$threads == '1', all(runtime$thread_environment == '1'))
    invisible(TRUE)
}

time_stage <- function(f, validate = function(x) invisible(TRUE)) {
    invisible(gc()) # GC explicitly requested before the timed interval, not within it.
    warnings <- character(); error <- ''; value <- NULL
    start <- proc.time()
    value <- tryCatch(withCallingHandlers(f(), warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w)); invokeRestart('muffleWarning')
    }), error = function(e) { error <<- conditionMessage(e); NULL })
    delta <- proc.time() - start
    # Validation cost is excluded, but invalid fits are not reported as successes.
    if (!nzchar(error)) tryCatch(validate(value), error = function(e) error <<- conditionMessage(e))
    list(value = value, status = if (nzchar(error)) 'error' else 'ok', error = error,
        warnings = paste(unique(warnings), collapse = ' | '),
        elapsed_seconds = unname(delta[['elapsed']]),
        user_seconds = unname(delta[['user.self']]), system_seconds = unname(delta[['sys.self']]))
}

validate_time_prediction <- function(p) {
    stopifnot(is.data.frame(p), nrow(p) == 201L,
        identical(names(p), c('estimate', 'mean', 'lower', 'upper')))
    for (name in names(p)) {
        p[[name]] <- normalize_probability(p[[name]], name)
        stopifnot(all(is.finite(p[[name]])))
    }
    check_curve(p$estimate); check_curve(p$mean)
    stopifnot(all(p$lower <= p$upper + 1e-8))
    invisible(TRUE)
}

time_one_method <- function(method, input, s, replicate, cfg = time_config()) {
    seed <- method_seed(replicate, s$setting_id, method)
    set.seed(seed)
    grid <- input$grid
    stopifnot(identical(grid, seq(0, input$tau, length.out = 201L)))
    rows <- list()
    row <- function(stage, z, burn = NA_integer_, keep = NA_integer_) data.frame(
        method = method, stage = stage, method_seed = seed, status = z$status,
        elapsed_seconds = z$elapsed_seconds, user_seconds = z$user_seconds,
        system_seconds = z$system_seconds, burn_in = burn, retained = keep,
        grid_points = length(grid), warnings = z$warnings, error = z$error,
        stringsAsFactors = FALSE)
    if (method %in% c('Unif', 'Exp', 'Weibull_PT')) {
        prior <- match(method, c('Unif', 'Exp', 'Weibull_PT'))
        map <- time_stage(function() {
            parameters <- if (prior == 3L) SurvPT::estimate(input$data$L, input$data$R) else NULL
            fit <- SurvPT::survPT(input$data[, c('L', 'R')], prior = prior,
                prior_params = parameters, flat = cfg$flat, kappa = cfg$kappa, gamma = cfg$gamma,
                max_M = cfg$max_M, max_EM_iter = cfg$max_EM_iter, EM_tol = cfg$EM_tol, verbose = FALSE)
            list(fit = fit, prediction = predict(fit, grid, sample = FALSE))
        }, function(v) check_curve(v$prediction$estimate))
        rows[[1L]] <- row('MAP-EM', map)
        if (map$status == 'ok') {
            gibbs <- time_stage(function() predict(map$value$fit, grid, sample = TRUE,
                n_samples = cfg$pt_sampling$n_samples, burn_in = cfg$pt_sampling$burn_in,
                thin = cfg$pt_sampling$thin, level = cfg$level), function(p) {
                    stopifnot(identical(attributes(p)[names(cfg$pt_sampling)], cfg$pt_sampling))
                    validate_time_prediction(p[, c('estimate', 'mean', 'lower', 'upper')])
                })
        } else gibbs <- list(status = 'dependency_error', error = paste('MAP-EM failed:', map$error),
            warnings = '', elapsed_seconds = NA_real_, user_seconds = NA_real_, system_seconds = NA_real_)
        rows[[2L]] <- row('DA-Gibbs', gibbs, 3000L, 2000L)
        total <- gibbs
        total$status <- if (map$status == 'ok' && gibbs$status == 'ok') 'ok' else 'error'
        total$warnings <- paste(unique(c(map$warnings, gibbs$warnings)), collapse = ' | ')
        for (field in c('elapsed_seconds', 'user_seconds', 'system_seconds'))
            total[[field]] <- if (total$status == 'ok') map[[field]] + gibbs[[field]] else NA_real_
        rows[[3L]] <- row('total', total, 3000L, 2000L)
        diagnostics <- if (map$status == 'ok') map$value$fit$Convergence else NULL
    } else {
        z <- time_stage(function() {
            if (method == 'Fiducial') {
                fit <- fiducial(input$data$L, input$data$R, s$n, grid = grid, testgrid = grid,
                    mfid = cfg$fiducial$mfid, mburn = cfg$fiducial$mburn, level = cfg$level)
                return(data.frame(estimate = 1 - fit$point_li, mean = 1 - fit$mean,
                    lower = 1 - fit$upper, upper = 1 - fit$lower))
            }
            d <- input$data[, c('L', 'R')]
            d$L[d$L == 0 & d$R > 0] <- NA_real_; d$R[is.infinite(d$R)] <- NA_real_
            # Exactly the main simulation adapters: include every parameter fit,
            # initialization chain, model validation and default package operation.
            fit <- if (method == 'MPT') fit_mpt(d) else if (method == 'TBP') fit_tbp(d) else
                stop('Unknown timing method: ', method)
            stopifnot(identical(fit$mcmc, cfg$mcmc), ncol(fit$theta.scaled) == cfg$mcmc$nsave)
            # The exact S(0)=1 boundary avoids log(0) in the package plotting routines.
            p <- if (method == 'MPT') plot(fit, xpred = matrix(0, 1L, 1L),
                tgrid = grid[-1L], CI = cfg$level, PLOT = FALSE) else
                plot(fit, tgrid = grid[-1L], CI = cfg$level, PLOT = FALSE)
            data.frame(estimate = c(1, p$Shat[, 1L]), mean = c(1, p$Shat[, 1L]),
                lower = c(1, p$Shatlow[, 1L]), upper = c(1, p$Shatup[, 1L]))
        }, validate_time_prediction)
        rows[[1L]] <- row('total', z, if (method == 'Fiducial') 100L else 3000L,
            if (method == 'Fiducial') 1000L else 2000L)
        diagnostics <- NULL
    }
    list(rows = do.call(rbind, rows), diagnostics = diagnostics)
}

time_summary <- function(raw, m) {
    rows <- list()
    for (i in seq_len(nrow(m$settings))) for (method in m$config$methods)
        for (stage in if (method %in% c('Unif', 'Exp', 'Weibull_PT')) c('MAP-EM', 'DA-Gibbs', 'total') else 'total') {
            s <- m$settings[i, ]
            z <- raw[raw$setting_id == s$setting_id & raw$method == method & raw$stage == stage, ]
            ok <- z$status == 'ok'
            x <- z$elapsed_seconds[ok]
            qs <- if (length(x)) quantile(x, c(.25, .5, .75), names = FALSE, type = 7) else rep(NA_real_, 3L)
            rows[[length(rows) + 1L]] <- data.frame(setting_id = s$setting_id,
                censoring = s$censoring, n = s$n, q = s$q, method = method, stage = stage,
                expected = m$replicates, attempted = nrow(z), successful = sum(ok),
                failed = sum(!ok), median_seconds = qs[2L], q1_seconds = qs[1L],
                q3_seconds = qs[3L], iqr_seconds = qs[3L] - qs[1L])
        }
    do.call(rbind, rows)
}

time_checkpoint_path <- function(root, task_id) file.path(root, 'checkpoints', sprintf('dataset_%04d.rds', task_id))
read_time_checkpoint <- function(root, task, m) {
    path <- time_checkpoint_path(root, task$task_id)
    if (!file.exists(path)) return(NULL)
    z <- readRDS(path)
    stopifnot(identical(z$manifest_id, m$id), identical(z$task, task),
        identical(z$data_hash, digest::digest(z$input, algo = 'sha256', serializeVersion = 3)),
        nrow(z$rows) == 12L, !anyDuplicated(z$rows[c('method', 'stage')]),
        setequal(z$rows$method, m$config$methods), all(z$rows$grid_points == 201L))
    z
}

report_time <- function(root, m = read_time(root, environment = FALSE), complete = FALSE) {
    saved <- lapply(seq_len(nrow(m$tasks)), function(i) read_time_checkpoint(root, m$tasks[i, , drop = FALSE], m))
    saved <- Filter(Negate(is.null), saved)
    if (!length(saved)) { message('No timing checkpoints yet.'); return(invisible(NULL)) }
    if (complete && length(saved) != nrow(m$tasks)) stop('Cannot finalize an incomplete benchmark')
    raw <- do.call(rbind, lapply(saved, `[[`, 'rows'))
    # Do not silently combine different machines or changed execution environments.
    stopifnot(length(unique(raw$host)) == 1L, length(unique(raw$cpu_model)) == 1L)
    summary <- time_summary(raw, m)
    atomic_write(raw, file.path(root, 'replicates.csv'), csv = TRUE)
    atomic_write(summary, file.path(root, 'summary.csv'), csv = TRUE)
    atomic_write(raw[raw$status != 'ok', ], file.path(root, 'failures.csv'), csv = TRUE)
    packages <- data.frame(package = names(m$packages), version = unname(m$packages))
    atomic_write(packages, file.path(root, 'package_versions.csv'), csv = TRUE)
    status <- if (length(saved) == nrow(m$tasks)) 'COMPLETE' else 'PARTIAL'
    lines <- c('# Computational time evaluation', '', paste('Status:', status), '',
        sprintf('Weibull(shape=2, scale=1/0.12); RC, CS, IC; n=200, 500, 1000; q=0.3, 0.5; %d datasets per combination.', m$replicates),
        sprintf('Completed datasets: %d/%d. All methods share each dataset. Seeds and hashes are saved.', length(saved), nrow(m$tasks)), '',
        paste('CPU:', unique(raw$cpu_model), '; host:', unique(raw$host)),
        paste('R:', m$R_version),
        'One pinned CPU, one R process, numerical-library threads=1, exclusive Slurm node allocation.', '',
        'All grids are seq(0,tau,length.out=201). MPT/TBP use exact S(0)=1 and evaluate the remaining 200 points with the package.',
        'Bayesian main chains: 3000 burn-in + 2000 retained sweeps, no thinning. MPT/TBP also retain the simulation default preliminary 5000+5000 parametric chain, included in time. Fiducial: 100 burn-in + 1000 retained samples.',
        'PT priors: Unif, Exp, Weibull_PT; flat=FALSE, kappa=0.1, gamma=0.5, max_M=8, max_EM_iter=2000, EM_tol=1e-6.', '',
        m$config$scope, m$config$pt_stages, m$config$other_stages, '',
        'Elapsed wall time is the primary metric. Raw user/system CPU time is also saved.',
        'Each process warms up all methods before measurement; warm-ups are excluded. Dataset order is seeded; method order rotates across datasets.',
        'Quartiles use R quantile(type=7). IQR=Q3-Q1. Totals are summed per dataset before computing quantiles.',
        'Errors retain their elapsed time in replicates.csv but are excluded from successful-fit summaries; no replacement datasets or retries.', '',
        '| Censoring | n | q | Method | Stage | Success/attempted | Median (s) | Q1 (s) | Q3 (s) | IQR (s) |',
        '|---|---:|---:|---|---|---:|---:|---:|---:|---:|')
    format_seconds <- function(x) ifelse(is.na(x), 'NA', sprintf('%.3f', x))
    for (i in seq_len(nrow(summary))) {
        s <- summary[i, ]
        lines <- c(lines, sprintf('| %s | %d | %.1f | %s | %s | %d/%d | %s | %s | %s | %s |',
            s$censoring, s$n, s$q, s$method, s$stage, s$successful, s$attempted,
            format_seconds(s$median_seconds), format_seconds(s$q1_seconds),
            format_seconds(s$q3_seconds), format_seconds(s$iqr_seconds)))
    }
    lines <- c(lines, '', 'Package versions:', '', '| Package | Version |', '|---|---|',
        sprintf('| %s | %s |', packages$package, packages$version), '',
        'Full sessionInfo(), CPU topology, Slurm allocation, affinity, load and thread settings are in sessions/.')
    writeLines(lines, file.path(root, 'report.md'))
    writeLines(c(status, sprintf('%d/%d datasets completed', length(saved), nrow(m$tasks)),
        paste('Updated:', format(Sys.time(), tz = 'America/New_York'))), file.path(root, 'STATUS.txt'))
    if (complete) atomic_write(list(manifest_id = m$id, completed = Sys.time(),
        datasets = length(saved), failed_stages = sum(raw$status != 'ok')),
        file.path(root, 'COMPLETE.rds'))
    invisible(summary)
}

worker_time <- function(root, limit = Inf) {
    root <- normalizePath(root, mustWork = TRUE)
    m <- read_time(root)
    do.call(RNGkind, as.list(m$RNGkind))
    runtime <- time_runtime(); assert_time_runtime(runtime)
    lock <- file.path(root, 'WORKER.lock')
    if (!dir.create(lock, showWarnings = FALSE)) stop('Timing worker already locked: ', lock)
    on.exit(unlink(lock, recursive = TRUE), add = TRUE)
    writeLines(c(paste('job', runtime$job_id), paste('pid', Sys.getpid())), file.path(lock, 'owner'))
    session_id <- paste0('job_', runtime$job_id, '_', Sys.getpid())
    session <- file.path(root, 'sessions', session_id)
    dir.create(session)
    atomic_write(runtime, file.path(session, 'runtime.rds'))
    writeLines(capture.output(runtime$session), file.path(session, 'sessionInfo.txt'))
    writeLines(capture.output(str(runtime)), file.path(session, 'runtime.txt'))
    writeLines(readLines('/proc/cpuinfo'), file.path(session, 'cpuinfo.txt'))
    # Check all previous records before doing any new work on a resumed job.
    for (i in seq_len(nrow(m$tasks))) {
        old <- read_time_checkpoint(root, m$tasks[i, , drop = FALSE], m)
        if (!is.null(old)) stopifnot(all(old$rows$host == runtime$host),
            all(old$rows$cpu_model == runtime$cpu_model))
    }
    # Warm every method with the real sweep budget, outside the measured results.
    warm_s <- m$settings[1L, ]; warm_s$n <- 20L
    warm_input <- generate_input(9999L, warm_s)
    warm <- lapply(m$config$methods, function(method) time_one_method(method, warm_input, warm_s, 9999L, m$config))
    names(warm) <- m$config$methods
    atomic_write(warm, file.path(session, 'warmup.rds'))
    assert_time_runtime(time_runtime())
    done <- 0L
    for (i in seq_len(nrow(m$tasks))) {
        task <- m$tasks[i, , drop = FALSE]
        old <- read_time_checkpoint(root, task, m)
        if (!is.null(old)) {
            stopifnot(all(old$rows$host == runtime$host), all(old$rows$cpu_model == runtime$cpu_model))
            next
        }
        if (done >= limit) break
        s <- m$settings[match(task$setting_id, m$settings$setting_id), ]
        input <- generate_input(task$replicate, s)
        data_hash <- digest::digest(input, algo = 'sha256', serializeVersion = 3)
        shift <- (i - 1L) %% length(m$config$methods)
        order <- m$config$methods[(seq_along(m$config$methods) + shift - 1L) %% length(m$config$methods) + 1L]
        values <- list(); loads <- list()
        for (method in order) {
            message('Dataset ', i, '/', nrow(m$tasks), ': ', s$label, ', replicate ', task$replicate, ', ', method)
            loads[[method]] <- list(before = readLines('/proc/loadavg'))
            values[[method]] <- time_one_method(method, input, s, task$replicate, m$config)
            loads[[method]]$after <- readLines('/proc/loadavg')
            stopifnot(proc_field('Threads') == '1', identical(proc_field('Cpus_allowed_list'), runtime$affinity))
            message('  ', paste(values[[method]]$rows$stage, values[[method]]$rows$status,
                sprintf('%.3fs', values[[method]]$rows$elapsed_seconds), collapse = '; '))
        }
        rows <- do.call(rbind, lapply(values, `[[`, 'rows')); rownames(rows) <- NULL
        rows$setting_id <- s$setting_id; rows$censoring <- s$censoring; rows$n <- s$n; rows$q <- s$q
        rows$replicate <- task$replicate; rows$data_seed <- task$replicate; rows$data_hash <- data_hash
        rows$task_id <- i; rows$tau <- input$tau; rows$host <- runtime$host
        rows$cpu_model <- runtime$cpu_model; rows$session_id <- session_id
        z <- list(manifest_id = m$id, task = task, input = input, data_hash = data_hash,
            rows = rows, diagnostics = lapply(values, `[[`, 'diagnostics'), load = loads)
        atomic_write(z, time_checkpoint_path(root, i))
        done <- done + 1L
        if (done == 1L || done %% 10L == 0L) report_time(root, m)
    }
    all_done <- all(file.exists(vapply(m$tasks$task_id, function(id) time_checkpoint_path(root, id), character(1))))
    report_time(root, m, complete = all_done)
}

time_main <- function(args = commandArgs(trailingOnly = TRUE)) {
    if (!length(args) || '--help' %in% args) {
        cat('Time_Evaluation.R prepare --run-dir=PARENT [--replicates=50]\n',
            'Time_Evaluation.R worker|report|status --time-dir=PARENT/time_evaluation\n')
        return(invisible(NULL))
    }
    mode <- args[1L]; opts <- list()
    allowed <- switch(mode, prepare = c('run-dir', 'replicates'), worker = 'time-dir',
        report = 'time-dir', status = 'time-dir', stop('Unknown mode: ', mode))
    for (arg in args[-1L]) {
        if (!grepl('^--[^=]+=.+$', arg)) stop('Use --name=value')
        key <- sub('^--([^=]+)=.*$', '\\1', arg)
        if (!key %in% allowed || key %in% names(opts)) stop('Unknown or duplicate option: ', key)
        opts[[key]] <- sub('^--[^=]+=', '', arg)
    }
    if (mode == 'prepare') {
        if (is.null(opts[['run-dir']])) stop('--run-dir is required')
        prepare_time(opts[['run-dir']], if (is.null(opts$replicates)) 50L else opts$replicates)
    } else {
        if (is.null(opts[['time-dir']])) stop('--time-dir is required')
        root <- opts[['time-dir']]
        if (mode == 'worker') worker_time(root)
        if (mode == 'report') report_time(root)
        if (mode == 'status') {
            m <- read_time(root, environment = FALSE)
            n <- sum(file.exists(vapply(m$tasks$task_id, function(id) time_checkpoint_path(root, id), character(1))))
            cat(n, '/', nrow(m$tasks), 'datasets checkpointed\n')
        }
    }
}
if (sys.nframe() == 0L) time_main()
