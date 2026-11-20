# Unweighted Winkler interval score for a central (1 - alpha) interval.
# Gneiting & Raftery (2007), equation (43). Smaller scores are better.
winkler_interval_score <- function(truth, lower, upper, level = .95) {
    stopifnot(is.numeric(level), length(level) == 1L, is.finite(level),
              level > 0, level < 1, is.numeric(truth), is.numeric(lower),
              is.numeric(upper), length(lower) == length(upper),
              length(truth) %in% c(1L, length(lower)),
              all(is.finite(c(truth, lower, upper))), all(lower <= upper))
    as.numeric((upper - lower) + 2 / (1 - level) *
        (pmax(lower - truth, 0) + pmax(truth - upper, 0)))
}
