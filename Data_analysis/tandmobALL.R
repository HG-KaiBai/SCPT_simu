# Run from Data_analysis: Rscript tandmobALL.R
# Read tandmobAll directly from the installed icensBKL package.
RNGkind('Mersenne-Twister', 'Inversion', 'Rejection')
set.seed(123)
source('analysis_helpers.R')
utils::data('tandmobAll', package = 'icensBKL')
output_dir <- 'results/tandmobALL'
cfg <- analysis_config()
record_environment(output_dir, cfg)
results <- list()
for (tooth in c(14L, 24L, 34L, 44L)) {
    tooth_results <- list()
    for (group in c('boy', 'girl')) {
        selected <- tandmobAll$GENDERNum == if (group == 'boy') 0 else 1
        d <- tandmobAll[which(selected), paste0(c('EBEG.', 'EEND.'), tooth)]
        names(d) <- c('L', 'R')
        z <- run_analysis_group(d, 'tandmobALL', group, cfg, tooth)
        results[[paste(tooth, group, sep = '_')]] <- tooth_results[[group]] <- z
        export_analysis(results, output_dir)
    }
    p <- plot_survival(do.call(rbind, lapply(tooth_results, `[[`, 'predictions')),
        paste('Tooth', tooth, ': emergence in boys and girls'), 'Age (years)', facet = TRUE)
    save_survival_plot(p, output_dir, paste0('tooth_', tooth), paired = TRUE)
}
