# Run from Data_analysis: Rscript HIV.R
RNGkind('Mersenne-Twister', 'Inversion', 'Rejection')
set.seed(123)
source('analysis_helpers.R')

# Original interval endpoints, preserved verbatim from HIV.R.
# Heavily treated group
L = c(15, 15, 10, 12, 9, 15, 5, 1, 15, 1, 9, 11, 12, 7, 11, 10, 12, 10, 1, 5, 10, 8,
           16, 16, 16, 1, 14, 14, 15, 15, 14, 9, 14, 14, 10, 13, 1, 11, 12, 9, 14, 3, 9, 13, 9,
           17, 17, 17, 12, 13, 1, 1, 14, 10, 10, 3, 15, 1, 10, 9, 11, 7, 13, 7, 10, 1, 1, 10,
           13, 14, 12, 10, 12, 10, 1, 12, 1, 9, 13, 7, 1, 9, 12, 13, 9, 10, 10,
           14, 13, 11, 5, 12, 8, 1, 14, 1, 1, 1, 14, 8, 3, 13, 10, 5, 10)

R = c(Inf, Inf, 11, 14, 12, 16, 7, 14, 16, 13, 10, 13, 13, 10, 13, 11, 13, 11, 7, 7, 11, 10,
      Inf, Inf, Inf, 16, 15, 15, 16, 15, 15, 10, 15, 15, 11, 14, 15, 13, 13, 13, 15, 7, 11, 15, 12,
      Inf, Inf, Inf, 13, 16, 11, 13, 15, 11, 12, 15, 16, 7, 15, 12, 13, 9, 14, 9, 11, 10, 7, 12,
      15, 15, 14, 11, 13, 11, 14, 13, 15, 12, 15, 15, 6, 12, 13, 15, 13, 12, 14,
      16, 15, 12, 7, 13, 10, 15, 15, 14, 11, 12, 16, 15, 14, 14, 12, 8, 12)
high <- data.frame(L = L, R = R)

# Lightly treated group
L = c(1, 10, 10, 8, 1, 1, 1, 12, 1, 6, 12, 1, 11, 11, 1, 1, 1, 15, 11, 1, 7,
      rep(15, 19), 12, 1, 1, 9, 1, 13, 9, 14, 13, 1, 13, 13, 13, 12, 12, 1, 13, 1, 12,
      rep(16, 31), 1, 15, 13, 8, 11, 9, 13, 1, 15, 12, 8, 3, 12, 1, 9, 12, 13, 11, 6,
      rep(17, 10), 1, 3, 5, 11, 8, 1, 4, 14, 7, 13, 10, 6, 7, 10, 12, 14, 10, 1, 11,
      18, 1, 8, 14, 1, 9, 1, 1, 1, 12, 15, 14, 15, 1, 11, 11, 1, 6, 10)

R = c(Inf, 15, 12, 13, 11, 18, 15, 15, 15, 13, 14, 15, 15, 14, 15, 14, 11, 16, 12, 12, 13,
      rep(Inf, 19), 14, 16, 12, 14, 15, 14, 12, 15, 14, 13, 15, 15, 14, 13, 13, 16, 13, 3, 13,
      rep(Inf, 31), 15, 16, 14, 16, 16, 14, 15, 12, 16, 13, 16, 16, 14, 15, 14, 13, 14, 14, 13,
      rep(Inf, 10), 15, 10, 11, 12, 12, 5, 11, 15, 12, 15, 12, 8, 10, 16, 14, 15, 12, 5, 14,
      Inf, 15, 15, 16, 17, 13, 16, 16, 14, 14, 16, 15, 16, 12, 14 ,12, 15, 12, 11)
low <- data.frame(L = L, R = R)

hiv_data <- list(high = high, low = low)
rm(L, R, high, low)

output_dir <- 'results/HIV'
cfg <- analysis_config()
record_environment(output_dir, cfg)
results <- list()
for (group in c('high', 'low')) {
    results[[group]] <- run_analysis_group(hiv_data[[group]], 'HIV', group, cfg)
    export_analysis(results, output_dir)
    p <- plot_survival(results[[group]]$predictions,
        paste0('HIV: ', group, ' treatment (n = ', nrow(hiv_data[[group]]), ')'), 'Time')
    save_survival_plot(p, output_dir, paste0('HIV_', group))
}
