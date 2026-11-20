#!/usr/bin/env Rscript
self <- sub('^--file=', '', grep('^--file=', commandArgs(), value = TRUE)[1L])
source(file.path(dirname(normalizePath(self)), 'Pipeline.R'))
must_fail <- function(expr) stopifnot(inherits(tryCatch({force(expr); NULL}, error = identity), 'error'))
s <- all_settings()[1L, ]
input <- generate_input(1L, s)
grid <- seq(0, 5, length.out = 100L)
times <- sort(unique(c(input$grid, input$quantiles)))
set.seed(method_seed(1L, s$setting_id, 'Fiducial'))
raw <- fiducial(input$data$L, input$data$R, s$n, grid, times, return_samples = TRUE)
v <- fit_method('Fiducial', input, s, 1L, scientific_config())
ix <- match(input$quantiles, times)
stopifnot(identical(dim(raw$cdf_samples), c(1000L, length(times))),
    isTRUE(all.equal(raw$mean, colMeans(raw$cdf_samples))),
    identical(v$curve, as.numeric(1 - raw$point_li[match(input$grid, times)])),
    identical(v$posterior$estimate, as.numeric(1 - raw$point_li[ix])),
    isTRUE(all.equal(v$posterior$mean, colMeans(1 - raw$cdf_samples[, ix]))),
    isTRUE(all.equal(v$posterior$lower, 1 - apply(raw$cdf_samples[, ix], 2L, sort)[975L, ])),
    isTRUE(all.equal(v$posterior$upper, 1 - apply(raw$cdf_samples[, ix], 2L, sort)[25L, ])),
    identical(v$sampling, list(mfid = 1000, mburn = 100)))
cat('PASS: fiducial mean from retained draws, interval endpoint reversal, original median point curve\n')

cat('ALL FIDUCIAL POSTERIOR CHECKS PASSED\n')
