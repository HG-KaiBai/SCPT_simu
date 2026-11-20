#!/usr/bin/env Rscript
self <- sub('^--file=', '', grep('^--file=', commandArgs(), value = TRUE)[1L])
source(file.path(dirname(normalizePath(self)), 'Pipeline.R'))
must_fail <- function(expr) stopifnot(inherits(tryCatch({force(expr); NULL}, error = identity), 'error'))

# Hand-calculated cases: covered, missed below/above, both endpoints and a point interval.
stopifnot(isTRUE(all.equal(winkler_interval_score(
    c(.5, .1, .9, .2, .8, .6), rep(.2, 6), c(rep(.8, 5), .2)),
    c(.6, 4.6, 4.6, .6, .6, 16))),
    isTRUE(all.equal(winkler_interval_score(.1, .2, .8, level = .9), 2.6)))
must_fail(winkler_interval_score(.5, .8, .2))
must_fail(winkler_interval_score(.5, .2, .8, level = 1))
must_fail(winkler_interval_score(.5, NA_real_, .8))
must_fail(winkler_interval_score(c(.2, .3), c(.1, .2, .3), c(.5, .6, .7)))
cat('PASS: Winkler scores, boundaries, non-default level and invalid inputs\n')

output <- tempfile('empt-posterior-verification-')
dir.create(output)
saveRDS(list(config = list(level = .95)), file.path(output, 'manifest.rds'))
input <- list(probabilities = c(.75, .5, .25), quantiles = c(1, 2, 3))
make_value <- function(lower, upper) list(
    posterior = posterior_table(input, (lower + upper) / 2, lower, upper),
    errors = c(SSE = 1, RMSE = .1, RMSE_nolast = .2), seconds = 1)
values <- list(make_value(c(.5, .6, .1), c(.9, .8, .2)),
               make_value(c(.8, .3, .1), c(.9, .7, .4)))
stopifnot(isTRUE(all.equal(values[[1L]]$posterior$winkler_score, c(.4, 4.2, 2.1))),
          isTRUE(all.equal(values[[1L]]$posterior$absolute_bias, c(.05, .2, .1))))
results <- lapply(seq_along(values), function(i) list(replicate = as.integer(i),
    main = list(Unif = values[[i]], MPT = values[[i]], TBP = values[[i]]),
    sensitivity = list(Unif_kappa0.1 = values[[i]])))
s <- data.frame(setting_id = 1L, label = 'test_setting')
for (group in c('main', 'sensitivity')) {
    root <- if (group == 'main') output else file.path(output, 'sensitivity')
    summary <- write_group(results, s, group, root)
    atomic_write(summary$posterior, file.path(root, 'posterior', 'summary.csv'), csv = TRUE)
}
files <- list.files(output, pattern = '\\.csv$', full.names = TRUE, recursive = TRUE)
files <- files[grepl('/posterior/', files, fixed = TRUE)]
expected <- setNames(lapply(files, read.csv), files)
for (path in files) {
    tab <- expected[[path]]
    if (basename(path) == 'summary.csv') {
        stopifnot(isTRUE(all.equal(tab$average_winkler_score,
            rep(c(1.25, 2.3, 1.2), nrow(tab) / 3L))),
            isTRUE(all.equal(tab$average_absolute_bias,
                rep(c(.075, .1, .05), nrow(tab) / 3L))))
    } else {
        stopifnot(identical(tab$replicate, c('1', '2', 'Average')),
            isTRUE(all.equal(unname(as.numeric(tail(tab, 1)[, paste0('S', c('75', '50', '25'), '_winkler_score')])), c(1.25, 2.3, 1.2))))
    }
}
unlink(output, recursive = TRUE)
cat('PASS: main/sensitivity exports, absolute bias without sign cancellation and score averages\n')
cat('ALL POSTERIOR METRIC CHECKS PASSED\n')
