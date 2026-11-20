#!/usr/bin/env Rscript
self <- sub('^--file=', '', grep('^--file=', commandArgs(), value = TRUE)[1L])
source(file.path(dirname(normalizePath(self)), 'Pipeline.R'))
must_fail <- function(expr) stopifnot(inherits(tryCatch({force(expr); NULL}, error = identity), 'error'))
input <- list(probabilities = c(.75, .5, .25), quantiles = c(1, 2, 3))
set.seed(123); before <- .Random.seed
v <- c(0, .2, 1, NA_real_, NaN, Inf, -Inf)
stopifnot(identical(normalize_probability(v), v), identical(.Random.seed, before),
    identical(normalize_probability(c(-1e-8, 1 + 1e-8, 1 + .Machine$double.eps)), c(0,1,1)))
must_fail(normalize_probability(c(-2e-8, .5)))
must_fail(normalize_probability(c(.5, 1 + 2e-8)))
p <- posterior_table(input, c(1 + 1e-9, .5, -1e-9),
    c(.9, .4, -1e-9), c(1 + 1e-9, .6, .1), estimate = c(1 + .Machine$double.eps, .6, -1e-9))
stopifnot(identical(p$estimate, c(1,.6,0)), identical(p$mean, c(1,.5,0)),
    identical(p$lower, c(.9,.4,0)), identical(p$upper, c(1,.6,.1)),
    identical(p$estimate_absolute_bias, abs(c(1,.6,0)-input$probabilities)),
    identical(p$absolute_bias, c(.25,0,.25)), identical(.Random.seed,before))
must_fail(posterior_table(input, c(1+2e-8,.5,.25), c(.9,.4,.1), c(1,.6,.4)))
must_fail(posterior_table(input, c(.75,.5,.25), c(.9,.4,.1), c(.8,.6,.4)))
p <- posterior_table(input, c(.75,NA,.25), c(.6,NA,.1), c(.9,NA,.4),
    estimate=c(1+.Machine$double.eps,.5,NA),allow_missing=TRUE)
stopifnot(p$estimate[1]==1,is.na(p$mean[2]),p$estimate[2]==.5,is.na(p$estimate[3]),p$mean[3]==.25)
cat('PASS: 1e-8 boundary tolerance, clipping, strict real-error rejection, missing targets and unchanged RNG\n')
# The previous successful outputs remain byte-for-byte identical after table rebuilding.
args <- commandArgs(trailingOnly=TRUE)
if (length(args)) {
    parent <- args[1L]
    files <- list.files(file.path(parent,'chunks'),pattern='[.]rds$',full.names=TRUE)
    for (path in files[unique(round(seq(1,length(files),length.out=min(6,length(files)))) )]) {
        z <- readRDS(path)$results[[1L]]
        for (group in c('main','sensitivity','gamma_sensitivity')) for (v in z[[group]]) {
            if (is_skipped_fit(v) || is.null(v$posterior)) next
            p <- v$posterior
            stopifnot(identical(p,posterior_table(z$input,p$mean,p$lower,p$upper,
                allow_missing=TRUE,estimate=p$estimate)))
        }
    }
    cat('PASS: previously successful posterior tables remain exactly unchanged\n')
}
