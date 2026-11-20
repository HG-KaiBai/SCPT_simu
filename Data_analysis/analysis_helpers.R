# Run from Data_analysis. Use the same implementations and parameters as simu_runs.
.simu_dir <- normalizePath('../simu_runs')
source(file.path(.simu_dir, 'Simulation.R'))
analysis_methods <- c('SurvPT', 'Fiducial', 'MPT', 'TBP', 'Turnbull')

analysis_config <- function() {
    cfg <- scientific_config()
    cfg$exp_spec <- subset(pt_methods(), method == 'Exp')
    cfg$grid_points <- 201L
    cfg
}

record_environment <- function(output_dir, cfg) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    cpu <- readLines('/proc/cpuinfo')
    cpu <- unique(trimws(sub('^[^:]+:', '', cpu[grepl('^model name[[:space:]]*:', cpu)])))
    packages <- c('SurvPT', 'spBayesSurv', 'interval', 'ggplot2', 'pracma', 'survival')
    versions <- sapply(packages, function(p) as.character(utils::packageVersion(p)))
    writeLines(c(paste('Started:', Sys.time()), paste('Slurm job:', Sys.getenv('SLURM_JOB_ID')),
        paste('Node:', Sys.info()[['nodename']]), paste('CPU:', paste(cpu, collapse = '; ')),
        paste('CPUs:', Sys.getenv('SLURM_CPUS_PER_TASK')), 'Numerical library threads: 1; seed: 123',
        'Runtime includes fitting, sampling and interval prediction; excludes plotting and file writes.',
        'Means average one fit per tooth/sex group (8 groups) or treatment group (2 groups).',
        capture.output(print(versions)), capture.output(utils::sessionInfo()),
        'Parameters from simu_runs:', capture.output(dput(cfg[c('exp_spec', 'max_M', 'max_EM_iter',
            'EM_tol', 'level', 'grid_points', 'survpt_defaults', 'survpt_sampling', 'fiducial', 'MPT', 'TBP')]))),
        file.path(output_dir, 'environment.txt'))
}

fit_analysis_method <- function(method, input, cfg) {
    d <- input$data
    x <- seq(0, input$tau, length.out = cfg$grid_points)
    if (method == 'SurvPT') {
        spec <- cfg$exp_spec
        fit <- SurvPT::survPT(d, prior = spec$prior, prior_params = NULL,
            flat = spec$flat, kappa = spec$kappa, gamma = spec$gamma,
            max_M = cfg$max_M, max_EM_iter = cfg$max_EM_iter, EM_tol = cfg$EM_tol, verbose = FALSE)
        p <- do.call(predict, c(list(object = fit, x = x, level = cfg$level, sample = TRUE),
                               cfg$survpt_sampling))
    } else if (method == 'Fiducial') {
        fit <- fiducial(d$L, d$R, nrow(d), grid = x, testgrid = x,
                       mfid = cfg$fiducial$mfid, mburn = cfg$fiducial$mburn, level = cfg$level)
        p <- data.frame(x = x, estimate = 1 - fit$point_li, mean = 1 - fit$mean,
                        lower = 1 - fit$upper, upper = 1 - fit$lower)
    } else if (method %in% c('MPT', 'TBP')) {
        d$L[d$L == 0 & d$R > 0] <- NA_real_
        d$R[is.infinite(d$R)] <- NA_real_
        fit <- if (method == 'MPT') fit_mpt(d) else fit_tbp(d)
        curves <- if (method == 'MPT') plot(fit, xpred = matrix(0, 1L, 1L),
            tgrid = x[-1L], CI = cfg$level, PLOT = FALSE) else
            plot(fit, tgrid = x[-1L], CI = cfg$level, PLOT = FALSE)
        p <- data.frame(x = x, estimate = c(1, curves$Shat[, 1L]), mean = c(1, curves$Shat[, 1L]),
                        lower = c(1, curves$Shatlow[, 1L]), upper = c(1, curves$Shatup[, 1L]))
    } else if (method == 'Turnbull') {
        fit <- interval::icfit(d$L, d$R, conf.int = FALSE)
        curve <- as.numeric(interval::getsurv(x, fit, 'interpolation')[[1L]]$S)
        p <- data.frame(x = x, estimate = curve, mean = curve, lower = NA_real_, upper = NA_real_)
    }
    p[, c('x', 'estimate', 'mean', 'lower', 'upper')]
}

run_analysis_group <- function(data, dataset, group, cfg, tooth = NA_integer_) {
    data <- data[!(is.na(data$L) & is.na(data$R)), c('L', 'R')]
    data$L[is.na(data$L)] <- 0
    data$R[is.na(data$R)] <- Inf
    endpoints <- unlist(data, use.names = FALSE)
    input <- list(data = data, tau = max(endpoints[is.finite(endpoints)]))
    rows <- timings <- list()
    for (method in analysis_methods) {
        message(dataset, ' / ', tooth, ' / ', group, ' / ', method)
        invisible(gc())
        started <- proc.time()[['elapsed']]
        p <- fit_analysis_method(method, input, cfg)
        elapsed <- proc.time()[['elapsed']] - started
        p$dataset <- dataset; p$tooth <- tooth; p$group <- group; p$method <- method; p$n <- nrow(data)
        p$level <- if (method == 'Turnbull') NA_real_ else cfg$level
        rows[[method]] <- p
        timings[[method]] <- data.frame(dataset = dataset, tooth = tooth, group = group,
            method = method, n = nrow(data), elapsed_seconds = elapsed)
        message(method, ': ', round(elapsed, 3), ' s')
    }
    list(predictions = do.call(rbind, rows), runtimes = do.call(rbind, timings))
}

export_analysis <- function(results, output_dir) {
    predictions <- do.call(rbind, lapply(results, `[[`, 'predictions'))
    times <- do.call(rbind, lapply(results, `[[`, 'runtimes'))
    summary <- do.call(rbind, lapply(analysis_methods, function(method) {
        x <- times$elapsed_seconds[times$method == method]
        data.frame(method = method, runs = length(x), total_seconds = sum(x), mean_seconds = mean(x))
    }))
    write.csv(predictions, file.path(output_dir, 'survival_intervals.csv'), row.names = FALSE, na = '')
    write.csv(times, file.path(output_dir, 'runtimes.csv'), row.names = FALSE, na = '')
    write.csv(summary, file.path(output_dir, 'runtime_summary.csv'), row.names = FALSE)
}

plot_survival <- function(predictions, title, x_label, facet = FALSE) {
    palette <- c(SurvPT = '#0072B2', Fiducial = '#D55E00', MPT = '#009E73', TBP = '#CC79A7', Turnbull = '#222222')
    predictions$method <- factor(predictions$method, levels = analysis_methods)
    summaries <- unique(predictions[, c('group', 'n')])
    groups <- summaries$group
    predictions$group <- factor(predictions$group, levels = groups)
    bands <- predictions[predictions$method != 'Turnbull', ]
    p <- ggplot2::ggplot(predictions, ggplot2::aes(x = x, y = estimate, colour = method)) +
        ggplot2::geom_line(data = bands, ggplot2::aes(y = lower), linetype = 'dashed', linewidth = 0.35) +
        ggplot2::geom_line(data = bands, ggplot2::aes(y = upper), linetype = 'dashed', linewidth = 0.35) +
        ggplot2::geom_line(linewidth = 0.75) +
        ggplot2::scale_colour_manual(values = palette, breaks = analysis_methods,
            labels = c('SurvPT (Exp)', 'Fiducial', 'MPT', 'TBP', 'Turnbull (icfit)'), name = 'Method') +
        ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(linetype = 'solid'), nrow = 1)) +
        ggplot2::coord_cartesian(ylim = c(0, 1)) +
        ggplot2::labs(title = title, subtitle = 'Solid: estimates; dashed: pointwise 95% intervals',
            x = x_label, y = 'Survival probability',
            caption = paste('SurvPT: MAP-EM; Fiducial: interpolated median; MPT/TBP: posterior mean.',
                'Bayesian credible intervals; Fiducial uses fiducial intervals. Turnbull: no interval.', sep = '\n')) +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::theme(legend.position = 'bottom', plot.caption = ggplot2::element_text(hjust = 0))
    if (facet) p <- p + ggplot2::facet_wrap(~group, nrow = 1L,
        labeller = ggplot2::as_labeller(setNames(paste0(groups, ' (n = ', summaries$n, ')'), groups)))
    p
}

save_survival_plot <- function(plot, output_dir, name, paired = FALSE) {
    width <- if (paired) 11 else 8
    ggplot2::ggsave(file.path(output_dir, paste0(name, '.png')),
        plot, width = width, height = 5.5, units = 'in', dpi = 220, bg = 'white')
}
