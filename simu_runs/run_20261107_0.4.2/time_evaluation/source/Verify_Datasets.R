#!/usr/bin/env Rscript
self <- sub('^--file=', '', grep('^--file=', commandArgs(), value = TRUE)[1L])
source(file.path(dirname(normalizePath(self)), 'Pipeline.R'))
settings <- all_settings(); rows <- list()
for (i in seq_len(nrow(settings))) {
    s <- settings[i, ]
    z <- lapply(1:3, function(t) list(replicate = t, input = generate_input(t, s)))
    rows[[i]] <- dataset_replicate_rows(z, s)
    for (v in z) stopifnot(identical(v$input$tau, attr(v$input$data, 'tau')),
        identical(v$input$censoring_proportions, attr(v$input$data, 'censoring_proportions')))
}
out <- tempfile('datasets-'); dir.create(out)
summary <- write_dataset_tables(do.call(rbind, rows), list(settings = settings, replicates = 3L), out)
for (i in seq_len(nrow(settings))) {
    s <- settings[i, ]; tau <- mean(rows[[i]]$tau)
    truth <- if (s$suite == 'mixture') .5 * exp(-(tau / 3)^3) + .5 * exp(-(tau / 10)^4) else
        if (s$suite == 'fiducial') switch(s$type, '1' = exp(-tau),
            '2' = exp(-tau) * (1 + tau + tau^2 / 2), '3' = exp(-tau) * (1 + tau), '4' = exp(-tau)) else
        exp(-(s$Lambda * tau)^s$K)
    stopifnot(abs(summary$true_S_at_average_tau[i] - truth) < 1e-12,
        summary$average_tau[i] == tau,
        abs(sum(unlist(summary[i, paste0('average_', c('Exact', 'Left', 'Interval', 'Right'))])) - 1) < 1e-12)
}
s <- settings[1L, ]; r <- rows[[1L]]
stopifnot(abs(summary$true_S_at_average_tau[1L] - mean(exp(-r$tau))) > 1e-8)
unlink(out, recursive = TRUE)
cat('PASS: all 64 dataset types, package tau/proportion attributes, analytic S(mean(tau)), not mean(S(tau))\n')
