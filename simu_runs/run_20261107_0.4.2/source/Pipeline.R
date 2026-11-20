#!/usr/bin/env Rscript
.source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(f) f$ofile))
.self <- if (length(.source_files)) tail(.source_files, 1L)[[1L]] else
    sub('^--file=', '', grep('^--file=', commandArgs(), value = TRUE)[1L])
.simu_dir <- dirname(normalizePath(.self, mustWork = TRUE))
source(file.path(.simu_dir, 'Simulation.R'), local = TRUE)
source(file.path(.simu_dir, 'Method_Failures.R'), local = TRUE)
source(file.path(.simu_dir, 'Dataset_Summary.R'), local = TRUE)
# Only the standalone unified workflow and its validation travel with a new run.
.sources <- c('Pipeline.R', 'Simulation.R', 'Fiducial.R', 'LinInterpolation.R',
    'Simulation_MPT.R', 'Simulation_TBP.R', 'Method_Failures.R', 'Dataset_Summary.R', 'Winkler.R',
    'Environment.sh', 'Install_SurvPT.sh', 'Submit.sh', 'Run.sbatch',
    'README.md', 'simu_guideline.md', 'Verify.R', 'Verify_Method_Failures.R',
    'Verify_Posterior.R', 'Verify_Fiducial_Posterior.R', 'Verify_SurvPT_Errors.R', 'Verify_Datasets.R', 'Verify_Roundoff.R')

positive_integer <- function(x, label) {
    x <- suppressWarnings(as.numeric(x))
    if (length(x) != 1L || !is.finite(x) || x < 1 || x != floor(x) || x > .Machine$integer.max)
        stop(label, ' must be a positive integer')
    as.integer(x)
}
package_versions <- function() {
    packages <- c('SurvPT', 'spBayesSurv', 'Rcpp', 'RcppArmadillo', 'interval', 'Icens',
                  'pracma', 'quadprog', 'survival', 'digest', 'coda', 'fields', 'spam', 'dotCall64', 'maps')
    setNames(vapply(packages, function(p) as.character(packageVersion(p)), character(1)), packages)
}
installed_survpt_archive <- function() {
    path <- file.path(find.package('SurvPT'), 'EMPT_ARCHIVE_MD5')
    if (!file.exists(path)) stop('Install the revised archive with Install_SurvPT.sh before running.')
    hash <- readLines(path, warn = FALSE)
    stopifnot(length(hash) == 1L, grepl('^[a-f0-9]{32}$', hash))
    hash
}
group_root <- function(root, group) if (group == 'main') root else file.path(root, group)


atomic_write <- function(x, path, csv = FALSE) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    tmp <- paste0(path, '.', Sys.getpid(), '.tmp')
    on.exit(unlink(tmp), add = TRUE)
    if (csv) write.csv(x, tmp, row.names = FALSE) else saveRDS(x, tmp, version = 3)
    if (!file.rename(tmp, path)) stop('Cannot save: ', path)
}

prepare_run <- function(run_dir, replicates = 1000L, chunk_size = 100L,
                        settings = all_settings(), analysis = 'all') {
    replicates <- positive_integer(replicates, 'replicates')
    chunk_size <- positive_integer(chunk_size, 'chunk-size')
    stopifnot(replicates <= 9999L, analysis %in% c('all', 'main', 'sensitivity', 'gamma'),
              nrow(settings) > 0L, !anyDuplicated(settings$setting_id))
    versions <- package_versions()
    config <- scientific_config()
    stopifnot(config$survpt_sampling$n_samples == 2000L, config$survpt_sampling$burn_in == 3000L,
              config$survpt_sampling$thin == 1L, config$survpt_defaults$kappa == 0.1)
    archive <- file.path(dirname(.simu_dir), 'SurvPT_0.4.2.tar.gz')
    if (!file.exists(archive)) archive <- file.path(dirname(.simu_dir), 'packages', 'SurvPT_0.4.2.tar.gz')
    if (!file.exists(archive)) stop('Missing package archive: ', archive)
    stopifnot(identical(installed_survpt_archive(), unname(tools::md5sum(archive))))
    if (file.exists(run_dir)) stop('Run directory already exists: ', run_dir)
    if (!dir.create(run_dir, recursive = TRUE)) stop('Cannot create run: ', run_dir)
    run_dir <- normalizePath(run_dir)
    for (sub in c('source', 'packages', 'chunks', 'logs')) dir.create(file.path(run_dir, sub))
    stopifnot(all(file.copy(file.path(.simu_dir, .sources), file.path(run_dir, 'source'))),
              file.copy(archive, file.path(run_dir, 'packages')))
    tasks <- do.call(rbind, lapply(settings$setting_id, function(id) {
        first <- seq.int(1L, replicates, by = chunk_size)
        data.frame(setting_id = id, first = first, last = pmin(first + chunk_size - 1L, replicates))
    }))
    tasks <- tasks[order(-settings$n[match(tasks$setting_id, settings$setting_id)],
                         tasks$first, tasks$setting_id), ]
    rownames(tasks) <- NULL
    tasks$task_id <- seq_len(nrow(tasks))
    hashes <- tools::md5sum(file.path(run_dir, 'source', .sources)); names(hashes) <- .sources
    manifest <- list(version = 6L, method_error_policy = method_error_policy, replicates = replicates, chunk_size = chunk_size,
        settings = settings, tasks = tasks, analysis = analysis, config = config,
        methods = pt_methods(TRUE, TRUE), RNGkind = RNGkind(), R_version = as.character(getRversion()),
        packages = versions, library_paths = .libPaths(), source_md5 = hashes,
        archive_md5 = unname(tools::md5sum(archive)), created = Sys.time(), session = sessionInfo())
    atomic_write(manifest, file.path(run_dir, 'manifest.rds'))
    for (name in c('settings', 'tasks', 'methods'))
        atomic_write(manifest[[name]], file.path(run_dir, paste0(name, '.csv')), csv = TRUE)
    writeLines(as.character(nrow(tasks)), file.path(run_dir, 'task_count.txt'))
    message('Prepared ', nrow(settings), ' settings x ', replicates, ' replicates; ',
            nrow(tasks), ' tasks; analysis=', analysis, ': ', run_dir)
    invisible(manifest)
}

read_run <- function(run_dir) {
    path <- file.path(run_dir, 'manifest.rds')
    m <- readRDS(path)
    stopifnot(identical(m$version, 6L), identical(m$R_version, as.character(getRversion())),
        identical(m$packages, package_versions()), identical(m$config, scientific_config()),
        identical(m$archive_md5, installed_survpt_archive()),
        identical(m$methods, pt_methods(TRUE, TRUE)),
        identical(unname(m$source_md5), unname(tools::md5sum(file.path(run_dir, 'source', names(m$source_md5))))),
        identical(m$method_error_policy, method_error_policy),
        identical(unname(m$source_md5), unname(tools::md5sum(file.path(.simu_dir, names(m$source_md5))))),
        identical(m$archive_md5, unname(tools::md5sum(file.path(run_dir, 'packages', 'SurvPT_0.4.2.tar.gz')))))
    m$id <- unname(tools::md5sum(path))
    do.call(RNGkind, as.list(m$RNGkind))
    m
}

validate_fit <- function(v, name, input, s, t, config) {
    actual <- canonical_method(name)
    if (is_skipped_fit(v)) {
        validate_skipped_fit(v, name, input = input, s = s, t = t)
        return(invisible(TRUE))
    }
    check_curve(v$curve)
    stopifnot(identical(v$errors, both_errors(v$curve, input$truth)),
        identical(v$seed, method_seed(t, s$setting_id, actual)),
        length(v$seconds) == 1L, is.finite(v$seconds), v$seconds >= 0)
    spec <- pt_methods(TRUE, TRUE); spec <- spec[spec$method == name, , drop = FALSE]
    if (name %in% c('OPT', 'NPMLE')) {
        stopifnot(is.null(v$posterior), is.null(v$status))
    } else {
        p <- v$posterior
        stopifnot(identical(p, posterior_table(input, p$mean, p$lower, p$upper,
            level = config$level, allow_missing = nrow(spec) > 0L, estimate = p$estimate)))
        if (anyNA(p[c('estimate', 'mean')])) {
            stopifnot(nrow(spec) == 1L)
            issue <- posterior_failure_info(actual, input, p)
            stopifnot(identical(v[names(issue)], issue))
        } else stopifnot(is.null(v$status))
        if (name == 'Fiducial') stopifnot(identical(v$sampling, config$fiducial))
        if (nrow(spec)) stopifnot(identical(v$sampling, config$survpt_sampling),
            identical(v$prior$family, c('uniform', 'exponential', 'weibull')[spec$prior]),
            identical(v$prior$flat, spec$flat), identical(v$prior$kappa, spec$kappa),
            identical(v$prior$gamma, spec$gamma))
        if (nrow(spec) && spec$prior == 3L) stopifnot(length(v$prior$parameters) == 2L,
            all(is.finite(v$prior$parameters)), all(v$prior$parameters > 0))
    }
    invisible(TRUE)
}
validate_result <- function(z, t, s, m, expected = expected_methods(s, m$analysis)) {
    stopifnot(is.list(z), identical(z$replicate, t), identical(z$setting_id, s$setting_id),
        identical(z$data_hash, digest::digest(z$input, algo = 'sha256', serializeVersion = 3)),
        identical(z$input, generate_input(t, s)))
    for (group in names(expected)) {
        stopifnot(identical(as.character(names(z[[group]])), expected[[group]]))
        for (name in expected[[group]]) validate_fit(z[[group]][[name]], name, z$input, s, t, m$config)
    }
    invisible(TRUE)
}

chunk_path <- function(run_dir, task_id) file.path(run_dir, 'chunks', sprintf('task_%04d.rds', task_id))
read_chunk <- function(run_dir, m, task) {
    path <- chunk_path(run_dir, task$task_id)
    if (!file.exists(path)) return(list(manifest_id = m$id, task = task, results = list()))
    chunk <- readRDS(path)
    stopifnot(identical(chunk$manifest_id, m$id), identical(chunk$task, task),
        !anyDuplicated(names(chunk$results)),
        all(names(chunk$results) %in% as.character(seq.int(task$first, task$last))))
    s <- m$settings[match(task$setting_id, m$settings$setting_id), ]
    for (id in names(chunk$results)) validate_result(chunk$results[[id]], as.integer(id), s, m)
    chunk
}

run_task <- function(run_dir, task_id, cores = 1L) {
    task_id <- positive_integer(task_id, 'task'); cores <- positive_integer(cores, 'cores')
    m <- read_run(run_dir)
    stopifnot(task_id <= nrow(m$tasks))
    task <- m$tasks[task_id, , drop = FALSE]
    s <- m$settings[match(task$setting_id, m$settings$setting_id), ]
    lock <- paste0(chunk_path(run_dir, task_id), '.lock')
    if (!dir.create(lock, showWarnings = FALSE)) stop('Task is locked (check for a live worker): ', lock)
    on.exit(unlink(lock, recursive = TRUE), add = TRUE)
    writeLines(c(paste('pid', Sys.getpid()), paste('job', Sys.getenv('SLURM_JOB_ID'))), file.path(lock, 'owner'))
    chunk <- read_chunk(run_dir, m, task)
    remaining <- setdiff(seq.int(task$first, task$last), as.integer(names(chunk$results)))
    for (start in seq_along(remaining)[(seq_along(remaining) - 1L) %% cores == 0L]) {
        wave <- remaining[start:min(length(remaining), start + cores - 1L)]
        values <- parallel::mclapply(wave, function(t) tryCatch({
            z <- simulate_replicate(t, s, m$analysis, m$config)
            validate_result(z, t, s, m)
            z
        }, error = function(e) list(failure = conditionMessage(e))),
            mc.cores = min(cores, length(wave)), mc.set.seed = FALSE)
        ok <- vapply(values, function(v) is.list(v) && !is.null(v$replicate), logical(1))
        for (j in which(ok)) chunk$results[[as.character(wave[j])]] <- values[[j]]
        atomic_write(chunk, chunk_path(run_dir, task_id))
        atomic_write(collect_method_failures(chunk$results, s),
            file.path(run_dir, 'method_errors', sprintf('task_%04d.csv', task_id)), csv = TRUE)
        message('Task ', task_id, ': saved ', length(chunk$results), '/', task$last - task$first + 1L)
        if (!all(ok)) stop('Failed replicates ', paste(wave[!ok], collapse = ','), ': ',
                          paste(vapply(values[!ok], function(v) paste(capture.output(str(v)), collapse = ' '), character(1)), collapse = '; '))
    }
    invisible(chunk)
}

run_status <- function(run_dir) {
    m <- read_run(run_dir)
    status <- do.call(rbind, lapply(seq_len(nrow(m$tasks)), function(i) {
        task <- m$tasks[i, , drop = FALSE]
        chunk <- read_chunk(run_dir, m, task)
        data.frame(task_id = i, completed = length(chunk$results), expected = task$last - task$first + 1L)
    }))
    cat('Completed replicates:', sum(status$completed), '/', sum(status$expected), '\n')
    cat('Incomplete task IDs:', paste(status$task_id[status$completed < status$expected], collapse = ','), '\n')
    invisible(status)
}

write_posterior <- function(posteriors, s, method, root, replicate_ids = seq_along(posteriors)) {
    rows <- lapply(posteriors, function(p) {
        z <- c()
        for (j in 1:3) {
            fields <- unlist(p[j, c('time', 'estimate', 'mean', 'lower', 'upper', 'estimate_absolute_bias', 'absolute_bias', 'length', 'covered')])
            names(fields) <- paste0('S', c('75', '50', '25')[j], '_', names(fields))
            z <- c(z, fields)
        }
        c(z, setNames(p$winkler_score, paste0('S', c('75', '50', '25'), '_winkler_score')))
    })
    tab <- do.call(rbind, rows)
    stopifnot(length(replicate_ids) == length(posteriors), !anyDuplicated(replicate_ids))
    out <- data.frame(replicate = c(as.character(replicate_ids), 'Average'),
                      rbind(tab, apply(tab, 2L, mean_available)), row.names = NULL)
    atomic_write(out, file.path(root, 'posterior', s$label, paste0(method, '.csv')), csv = TRUE)
    summary <- lapply(1:3, function(j) {
        p <- do.call(rbind, lapply(posteriors, function(p) p[j, ]))
        data.frame(setting_id = s$setting_id, setting = s$label, method = method,
            probability = p$probability[1L], time = p$time[1L], replicates = sum(!is.na(p$mean)),
            estimate_replicates = sum(!is.na(p$estimate)), estimate = mean_available(p$estimate),
            average_estimate_absolute_bias = mean_available(p$estimate_absolute_bias),
            posterior_mean = mean_available(p$mean), average_absolute_bias = mean_available(p$absolute_bias),
            average_length = mean_available(p$length), coverage_rate = mean_available(p$covered),
            average_winkler_score = mean_available(p$winkler_score))
    })
    do.call(rbind, summary)
}

write_group <- function(results, s, group, root) {
    methods <- names(results[[1L]][[group]])
    if (!length(methods)) return(NULL)
    point_rows <- list(); point_summary <- list(); posterior_summary <- list()
    for (method in methods) {
        values <- lapply(results, function(z) z[[group]][[method]])
        errors <- do.call(rbind, lapply(values, `[[`, 'errors'))
        point_rows[[method]] <- data.frame(replicate = vapply(results, `[[`, integer(1), 'replicate'), method = method, errors,
            seconds = vapply(values, `[[`, numeric(1), 'seconds'))
        summary <- data.frame(setting_id = s$setting_id, setting = s$label, method = method,
                              replicates = sum(!vapply(values, is_skipped_fit, logical(1))))
        for (metric in colnames(errors)) {
            summary[[paste0(metric, '_mean')]] <- mean_available(errors[, metric])
            summary[[paste0(metric, '_sd')]] <- sd_available(errors[, metric])
        }
        point_summary[[method]] <- summary
        if (!is.null(values[[1L]]$posterior)) {
            posterior_summary[[method]] <- write_posterior(
                lapply(values, `[[`, 'posterior'), s, method, root,
                replicate_ids = vapply(results, `[[`, integer(1), 'replicate'))
        }
    }
    atomic_write(do.call(rbind, point_rows), file.path(root, 'point_estimate', s$label, 'replicates.csv'), csv = TRUE)
    atomic_write(do.call(rbind, point_summary), file.path(root, 'point_estimate', s$label, 'summary.csv'), csv = TRUE)
    list(point = do.call(rbind, point_summary), posterior = do.call(rbind, posterior_summary))
}

merge_run <- function(run_dir) {
    m <- read_run(run_dir)
    # Refuse incomplete or invalid runs before publishing any summaries.
    for (i in seq_len(nrow(m$tasks))) {
        task <- m$tasks[i, , drop = FALSE]
        chunk <- read_chunk(run_dir, m, task)
        if (!setequal(names(chunk$results), as.character(seq.int(task$first, task$last))))
            stop('Incomplete task ', i)
    }
    summaries <- list(main = list(), sensitivity = list(), gamma_sensitivity = list())
    failures <- status <- dataset_rows <- list()
    for (i in seq_len(nrow(m$settings))) {
        s <- m$settings[i, ]
        tasks <- m$tasks[m$tasks$setting_id == s$setting_id, , drop = FALSE]
        results <- unlist(lapply(seq_len(nrow(tasks)), function(j)
            read_chunk(run_dir, m, tasks[j, , drop = FALSE])$results), recursive = FALSE)
        stopifnot(!anyDuplicated(names(results)), setequal(names(results), as.character(seq_len(m$replicates))))
        results <- results[as.character(seq_len(m$replicates))]
        dataset_rows[[i]] <- dataset_replicate_rows(results, s)
        failures[[i]] <- collect_method_failures(results, s)
        status[[i]] <- collect_method_status(results, s, m)
        for (group in names(summaries)) {
            root <- group_root(run_dir, group)
            summaries[[group]][[as.character(s$setting_id)]] <- write_group(results, s, group, root)
        }
    }
    for (group in names(summaries)) if (length(summaries[[group]])) {
        root <- group_root(run_dir, group)
        for (kind in c('point', 'posterior')) {
            tab <- do.call(rbind, lapply(summaries[[group]], `[[`, kind))
            if (!is.null(tab)) atomic_write(tab,
                file.path(root, if (kind == 'point') 'point_estimate' else 'posterior', 'summary.csv'), csv = TRUE)
        }
    }
    write_dataset_tables(do.call(rbind, dataset_rows), m, run_dir)
    failures <- do.call(rbind, failures); status <- do.call(rbind, status)
    write_method_failure_reports(failures, status, run_dir)
    atomic_write(list(manifest_id = m$id, completed = Sys.time(),
        method_error_policy = method_error_policy, probability_bounds = m$config$probability_bounds,
        failed_method_fits = nrow(failures),
        replicates = nrow(m$settings) * m$replicates), file.path(run_dir, 'COMPLETE.rds'))
    message('Merged: ', run_dir)
    invisible(TRUE)
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
    if (!length(args) || '--help' %in% args) {
        cat('Rscript Pipeline.R prepare|worker|merge|status|datasets|errors --run-dir=PATH\n',
            'prepare: --replicates=1000 --chunk-size=100 --suite=all|exp|weibull|fiducial|mixture\n',
            '         --analysis=all|main|sensitivity|gamma [--settings=1,29]\n',
            'worker: --task=N --cores=N (defaults: Slurm task ID; 1 core)\n',
            'Use RUN/source/Pipeline.R for worker, merge and status.\n')
        return(invisible(NULL))
    }
    mode <- args[1L]
    allowed <- switch(mode, prepare = c('run-dir', 'replicates', 'chunk-size', 'suite', 'analysis', 'settings'),
        worker = c('run-dir', 'task', 'cores'), merge = 'run-dir', status = 'run-dir',
        datasets = 'run-dir', errors = 'run-dir', stop('Unknown mode'))
    opts <- list()
    for (arg in args[-1L]) {
        if (!grepl('^--[^=]+=.+$', arg)) stop('Use --name=value')
        key <- sub('^--([^=]+)=.*$', '\\1', arg)
        if (!key %in% allowed || key %in% names(opts)) stop('Unknown or duplicate option: ', key)
        opts[[key]] <- sub('^--[^=]+=', '', arg)
    }
    value <- function(key, default) if (is.null(opts[[key]])) default else opts[[key]]
    if (is.null(opts[['run-dir']])) stop('--run-dir is required')
    run_dir <- opts[['run-dir']]
    switch(mode, prepare = {
        settings <- all_settings()
        suite <- value('suite', 'all')
        stopifnot(suite %in% c('all', 'exp', 'weibull', 'fiducial', 'mixture'))
        if (suite != 'all') settings <- settings[settings$suite == suite, ]
        if (!is.null(opts$settings)) {
            ids <- vapply(strsplit(opts$settings, ',', fixed = TRUE)[[1L]], positive_integer, integer(1), label = 'settings')
            stopifnot(all(ids %in% settings$setting_id))
            settings <- settings[settings$setting_id %in% ids, ]
        }
        prepare_run(run_dir, value('replicates', 1000L), value('chunk-size', 100L),
                    settings, value('analysis', 'all'))
    }, worker = run_task(run_dir, value('task', Sys.getenv('SLURM_ARRAY_TASK_ID')), value('cores', 1L)),
       merge = merge_run(run_dir), status = run_status(run_dir),
       datasets = generate_dataset_report(run_dir), errors = report_run_errors(run_dir))
}
if (sys.nframe() == 0L) main()
