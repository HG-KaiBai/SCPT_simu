# Package attributes are averaged over datasets, independently of fit success.
dataset_replicate_rows <- function(results, s) {
    do.call(rbind, lapply(results, function(z) {
        input <- z$input
        data.frame(setting_id = s$setting_id, setting = s$label, replicate = z$replicate,
            tau = input$tau, as.list(input$censoring_proportions),
            total_censored = 1 - unname(input$censoring_proportions['Exact']))
    }))
}
write_dataset_tables <- function(rows, m, run_dir) {
    stopifnot(nrow(rows) == nrow(m$settings) * m$replicates,
        !anyDuplicated(rows[c('setting_id', 'replicate')]))
    rows <- rows[order(rows$setting_id, rows$replicate), ]; rownames(rows) <- NULL
    summaries <- lapply(seq_len(nrow(m$settings)), function(i) {
        s <- m$settings[i, ]; r <- rows[rows$setting_id == s$setting_id, ]
        stopifnot(setequal(r$replicate, seq_len(m$replicates)), all(r$setting == s$label))
        average_tau <- mean(r$tau)
        data.frame(s[c('setting_id', 'label', 'suite', 'n', 'type', 'K', 'Lambda', 'q')],
            replicates = nrow(r), average_tau = average_tau,
            true_S_at_average_tau = true_survival(average_tau, s),
            as.list(setNames(vapply(r[c('Exact', 'Left', 'Interval', 'Right', 'total_censored')],
                mean, numeric(1)), paste0('average_', c('Exact', 'Left', 'Interval', 'Right', 'total_censored')))))
    })
    atomic_write(rows, file.path(run_dir, 'datasets', 'replicates.csv'), csv = TRUE)
    atomic_write(do.call(rbind, summaries), file.path(run_dir, 'datasets', 'summary.csv'), csv = TRUE)
    invisible(do.call(rbind, summaries))
}
generate_dataset_report <- function(run_dir) {
    m <- read_run(run_dir)
    rows <- lapply(seq_len(nrow(m$settings)), function(i) {
        s <- m$settings[i, ]
        # Reproduce exactly the seeds used by workers; no fitting or selection.
        results <- lapply(seq_len(m$replicates), function(t)
            list(replicate = t, input = generate_input(t, s)))
        out <- dataset_replicate_rows(results, s)
        message('Datasets: ', i, '/', nrow(m$settings), ' ', s$label)
        out
    })
    write_dataset_tables(do.call(rbind, rows), m, run_dir)
    message('Dataset report: ', file.path(run_dir, 'datasets', 'summary.csv'))
    invisible(TRUE)
}
report_run_errors <- function(run_dir) {
    m <- read_run(run_dir)
    failures <- status <- list()
    for (i in seq_len(nrow(m$settings))) {
        s <- m$settings[i, ]; tasks <- m$tasks[m$tasks$setting_id == s$setting_id, , drop = FALSE]
        results <- unlist(lapply(seq_len(nrow(tasks)), function(j)
            read_chunk(run_dir, m, tasks[j, , drop = FALSE])$results), recursive = FALSE)
        failures[[i]] <- collect_method_failures(results, s)
        status[[i]] <- collect_method_status(results, s, m)
    }
    failures <- do.call(rbind, failures)
    write_method_failure_reports(failures, do.call(rbind, status), run_dir)
    message('Reported ', nrow(failures), ' method issues from completed checkpoints')
    invisible(TRUE)
}
