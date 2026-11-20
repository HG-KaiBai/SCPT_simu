# SCPT_simu

SCPT_simu compares methods for estimating survival curves from censored event times. It includes an R simulation study of **SurvPT 0.4.2**, a single-core timing benchmark, and two applications to interval-censored data. SurvPT fits a finite Pólya tree by maximum a posteriori estimation with EM and uses a data-augmentation Gibbs sampler for pointwise credible intervals.

The companion [SurvPT](https://github.com/HG-KaiBai/SurvPT) repository contains the R package source code, installation instructions, and usage documentation.

## Explore the results

| Study | What was measured | Saved results |
| --- | --- | --- |
| Simulation | Survival-curve error and 95% interval coverage across 64 settings and 1,000 replications per setting | [Point estimates](simu_runs/run_20261107_0.4.2/point_estimate/summary.csv) · [Interval summaries](simu_runs/run_20261107_0.4.2/posterior/summary.csv) |
| Timing | Fitting and prediction time for six methods on 900 datasets using one CPU | [Timing report](simu_runs/run_20261107_0.4.2/time_evaluation/report.md) |
| HIV | Survival by treatment group, using the interval endpoints included in the analysis script | [Heavily treated group](Data_analysis/results/HIV/HIV_high.png) · [Lightly treated group](Data_analysis/results/HIV/HIV_low.png) · [Estimates](Data_analysis/results/HIV/survival_intervals.csv) |
| Tooth emergence | Age at emergence for four teeth, by sex, using `tandmobAll` from the R package `icensBKL` | [Plots and tables](Data_analysis/results/tandmobALL/) |

The CSVs and plots can be read directly. The [completed simulation](simu_runs/run_20261107_0.4.2/) occupies about 17 GB including worker checkpoints. Its manifest records the package versions, random-number settings and checksums of the source used for the run.

## Reproduce the simulation

Run commands from the repository root. The study used R 4.4.0, a C++17 compiler to build SurvPT, and Slurm for batch jobs. The [simulation guide](simu_runs/README.md) covers environment setup, a one-setting pilot, the full study, and how to interpret the outputs.

On the original cluster or with the scripts configured for your site, this command prepares the full design in a new directory and prints its path. It stops before submitting any jobs.

```bash
bash Submit_Simulation.sh --prepare-only
```

## Reproduce the data analyses

Both analyses compare SurvPT with exponential centering, fiducial inference, two Bayesian survival models (MPT and TBP), and the Turnbull estimator. Each results directory contains the survival estimates, pointwise 95% intervals for the first four methods, plots, per-fit runtimes and the recorded R environment. The plots used `ggplot2` 4.0.3, and the tooth analysis requires `icensBKL`. [HIV script](Data_analysis/HIV.R) · [Tooth-emergence script](Data_analysis/tandmobALL.R) · [Shared fitting code](Data_analysis/analysis_helpers.R)

The scripts write to `results/` beneath their working directory. To keep the archived outputs available for comparison, create a fresh working directory on a filesystem shared with the compute nodes. Run the following from the repository root after configuring the environment in the simulation guide.

```bash
work=$(mktemp -d "$PWD/data-rerun.XXXXXX")
mkdir "$work/Data_analysis"
cp Data_analysis/*.R Data_analysis/Submit.sh Data_analysis/Run.sbatch "$work/Data_analysis/"
ln -s "$PWD/simu_runs" "$work/simu_runs"
bash "$work/Data_analysis/Submit.sh"
```

This submits one Slurm job per dataset and writes the outputs under `$work/Data_analysis/results/`. To run one dataset, pass `HIV` or `tandmobALL` to the submission script. With the same R libraries available on a Linux machine, you can also run `Rscript HIV.R` or `Rscript tandmobALL.R` from the fresh `Data_analysis/` working directory.
