# Baseline transformed Bernstein polynomial: intercept-only model, with the
# checked initializer recovery shared with Simulation_MPT.R (source it first).
fit_tbp <- function(data) {
    f <- adapt_bayes_initializers(spBayesSurv::survregbayes, bayes_survreg)
    fit <- f(survival::Surv(L, R, type = 'interval2') ~ 1, data = data)
    stopifnot(identical(fit$mcmc, eval(formals(spBayesSurv::survregbayes)$mcmc)),
              fit$maxL == 15, fit$prior$a0 == 1, fit$prior$b0 == 1)
    fit
}

curve_tbp <- function(fit, x) {
    ans <- plot(fit, tgrid = x[x > 0], PLOT = FALSE)
    c(1, as.numeric(ans$Shat[, 1L]))
}
