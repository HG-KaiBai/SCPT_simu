# Simulation and timing guide

The reference run, `run_20261107_0.4.2`, contains 64,000 simulated datasets, the fitted results and the code that produced them. Shell examples below start at the repository root. To explore the study without fitting a model, open the [setting definitions](run_20261107_0.4.2/settings.csv), [point-estimate summary](run_20261107_0.4.2/point_estimate/summary.csv) or [posterior summary](run_20261107_0.4.2/posterior/summary.csv).

## Environment setup

The reference environment used **R 4.4.0**, **SurvPT 0.4.2** and **spBayesSurv 1.1.9**. Use the bundled SurvPT source archive and a C++17 compiler. Install the other R packages and their dependencies at the versions recorded for the study. The pipeline checks package versions, the installed SurvPT archive checksum and the source checksums when it reads a saved run. [Package versions](run_20261107_0.4.2/time_evaluation/package_versions.csv) · [SurvPT source archive](../SurvPT_0.4.2.tar.gz)

**On the original cluster**, the installer loads the R module, installs SurvPT in the study library and records the archive checksum. Load the environment in each shell where you use `Rscript` directly.

```bash
bash simu_runs/Install_SurvPT.sh
source simu_runs/Environment.sh
```

**On another Linux system**, use a writable R library for the pinned packages. For example, once the dependencies are installed and visible to R 4.4.0, the following installs SurvPT and records the checksum used by the pipeline. Run the local pilot below in the same shell.

```bash
export R_LIBS_USER="$HOME/R/empt-4.4"
mkdir -p "$R_LIBS_USER"
R CMD INSTALL --library="$R_LIBS_USER" SurvPT_0.4.2.tar.gz
md5sum SurvPT_0.4.2.tar.gz | cut -d ' ' -f 1 > "$R_LIBS_USER/SurvPT/EMPT_ARCHIVE_MD5"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
```

For submission on another Slurm cluster, adapt the module and R library paths in the working copies of [Environment.sh](Environment.sh) and [Install_SurvPT.sh](Install_SurvPT.sh), and the account in the applicable [simulation](Run.sbatch), [timing](Time_Run.sbatch) or [data-analysis](../Data_analysis/Run.sbatch) batch script before preparing a new run. The archived `source/` directories participate in checksum verification.

## Run a pilot

This example uses the archived pipeline, one dataset from the first fiducial scenario, and one CPU. The worker fits the main methods and the merge produces the same table layout as a full run. Choose a fresh directory name each time.

```bash
reference=simu_runs/run_20261107_0.4.2
pilot=simu_runs/pilot_01
Rscript "$reference/source/Pipeline.R" prepare --run-dir="$pilot" \
  --settings=1 --analysis=main --replicates=1 --chunk-size=1
Rscript "$pilot/source/Pipeline.R" worker --run-dir="$pilot" --task=1 --cores=1
Rscript "$pilot/source/Pipeline.R" merge --run-dir="$pilot"
```

Inspect `$pilot/point_estimate/summary.csv` and `$pilot/posterior/summary.csv`. With one replication, coverage at each target is the outcome of a single interval and the between-replication standard deviation is undefined.

## Run the full study on Slurm

On the configured cluster, use a new shared directory. The first command prepares the run and the second submits that same run. The submission launches a dataset summary job, a worker array and a merge job that starts after the array succeeds.

```bash
run=simu_runs/full_01
EMPT_RUN_DIR="$run" bash Submit_Simulation.sh --prepare-only
bash Submit_Simulation.sh --submit-run="$run"
```

To replay the archived source on the original cluster, use `bash simu_runs/run_20261107_0.4.2/source/Submit.sh` in place of `bash Submit_Simulation.sh` in both commands, keeping the fresh output path. For a new study with the working source and an automatically named directory, `bash Submit_Simulation.sh` prepares and submits in one step.

| Control | Default | Purpose |
| --- | --- | --- |
| `EMPT_REPLICATES` / `EMPT_CHUNK_SIZE` | 1,000 / 100 | Replications per setting and per array task. The full design creates 640 tasks. |
| `EMPT_SUITE` | `all` | Select `exp`, `weibull`, `fiducial` or `mixture`. The preparation option `--settings=1,29,53` selects specific IDs. |
| `EMPT_ANALYSIS` | `all` | Select `main`, `sensitivity` (kappa), or `gamma`. |
| `EMPT_CPUS_PER_TASK` / `EMPT_MEMORY` / `EMPT_WALLTIME` | 16 / 32G / 12 hours | Resources for each worker. The merge uses one CPU, 8G and 12 hours. |
| `EMPT_CONCURRENT_TASKS` | 0 | Set a positive integer to limit simultaneous array tasks. |

Set design options at preparation and Slurm resource options at submission. Workers checkpoint completed replications. To check a run, use `Rscript "$run/source/Pipeline.R" status --run-dir="$run"` in the matching R environment. This validates saved replications and can take time on the full archive. A quick read of `COMPLETE.rds` gives the recorded completion count. To resume an interrupted task, invoke the run's frozen worker with the original task ID after the previous worker has stopped. The merge validates all expected tasks before writing summaries.

## Study design

| Dataset family | Settings | Parameter grid |
| --- | ---: | --- |
| Fiducial scenarios | 28 | Four scenarios at sample sizes 50, 75, 100, 200, 500, 800 and 1,000 |
| Exponential and Weibull | 24 | Shape 1 or 2, inverse scale 0.12, sample size 200 or 500, and `q` of 0.3 or 0.5 under right, current-status and interval censoring |
| Weibull mixture | 12 | Equal weights on Weibull(shape 3, scale 3) and Weibull(shape 4, scale 10), with the same sample sizes, censoring types and `q` grid |

SurvPT uses uniform, exponential and Weibull centering, each with a flat variant. The comparators are fiducial inference, the nonparametric maximum likelihood estimator (NPMLE), the mixture of Pólya trees (MPT), and the transformed Bernstein polynomial (TBP) model. Right-censored settings also include the optional Pólya tree (OPT). The 201-point prediction grid runs from zero to the largest finite observed endpoint, `tau`. The three interval targets are the true times at which survival equals 0.75, 0.50 and 0.25.

The main SurvPT fits use `kappa = 0.1`, `gamma = 0.5`, maximum tree depth 8, at most 2,000 EM iterations and tolerance `1e-6`. The Weibull centering parameters are estimated once from the censored observations. Point curves use the MAP-EM fit. Posterior intervals use 3,000 burn-in draws and 2,000 retained draws with thinning 1. Fiducial fits use 100 burn-in and 1,000 retained draws. The kappa sensitivity study evaluates 0.05, 0.1, 0.2 and 0.3 across all six variants. The gamma study evaluates 0.5, 1 and 2 across the three non-flat variants.

Methods share the dataset generated with the replication number as its seed. Each fitting seed is `100000000 + 1000000 × method_id + 10000 × setting_id + replicate`. Matching sensitivity settings reuse the main fit. The manifest records the random-number generator and fitted-model configuration. [Method definitions](run_20261107_0.4.2/methods.csv) · [Scientific implementation](run_20261107_0.4.2/source/Simulation.R)

## Read the output

| File or directory inside a run | Interpretation |
| --- | --- |
| `point_estimate/summary.csv` | Mean and SD by setting and method. SSE sums squared survival errors across the grid. RMSE compares absolute changes between adjacent survival values, including the final drop to zero. `RMSE_nolast` uses only the observed grid. |
| `posterior/summary.csv` | Mean interval length, 95% coverage and Winkler score at each target. `average_estimate_absolute_bias` uses the method's point prediction. `average_absolute_bias` uses its posterior or fiducial mean. Both average the absolute error per replication. |
| `sensitivity/`, `gamma_sensitivity/` | Point and posterior tables for the kappa and gamma grids, respectively. |
| `datasets/summary.csv` | Censoring proportions and mean `tau` by setting. `true_S_at_average_tau` evaluates the population survival function at that mean. |
| `fit_status_summary.csv`, `fit_failures.csv` | Expected and successful fit counts, target-specific counts and details of skipped or partially missing fits. |
| `manifest.rds`, `source/`, `packages/`, `chunks/`, `COMPLETE.rds` | Run configuration, frozen source and package, worker checkpoints, and the completion record. |

Point predictions at the interval targets are the SurvPT MAP-EM estimate, the interpolated fiducial median, or the MPT/TBP posterior mean. Coverage counts only available intervals. In the posterior summary, `replicates` counts valid posterior summaries and `estimate_replicates` counts valid point predictions for each target. A missing SurvPT prediction beyond `tau` can therefore coexist with valid results at other targets. Finite posterior probabilities within `1e-8` of a boundary are clipped if they fall outside [0, 1]. Larger excursions trigger validation errors.

## Single-core timing benchmark

The benchmark measures fitting and prediction for six methods on 18 Weibull settings, with 50 datasets per setting. One R process uses one pinned CPU on an exclusive node. SurvPT reports MAP-EM and Gibbs stages separately. MPT and TBP timings include their preliminary parametric chains. The reported medians and interquartile ranges use successful measurements and exclude warm-up and data generation. [Protocol](run_20261107_0.4.2/time_evaluation/protocol.md) · [Timing summary](run_20261107_0.4.2/time_evaluation/summary.csv) · [Full report](run_20261107_0.4.2/time_evaluation/report.md)

To run a new benchmark, use the new parent directory from the full-study example. These commands create its `time_evaluation/` directory and submit the prepared benchmark.

```bash
bash simu_runs/Time_Submit.sh --prepare-only --run-dir="$run"
bash simu_runs/Time_Submit.sh --submit-time="$run/time_evaluation"
```

The timing worker requests 16G and 24 hours. Set `EMPT_TIME_REPLICATES` during preparation for a smaller benchmark or `EMPT_TIME_WALLTIME` at submission to adjust the allocation. The timing manifest records the absolute path of its parent run. A fresh preparation records the path on your system. The archived report and CSVs can be read from any checkout.

## Validation

The archived source includes checks for dataset generation, posterior summaries, method failures, numerical roundoff, checkpoint handling, and serial versus parallel execution. In the configured R environment, run the focused boundary check below. The broader integration check is `Verify.R` in the same directory. [Archived validation scripts](run_20261107_0.4.2/source/)

```bash
Rscript simu_runs/run_20261107_0.4.2/source/Verify_Roundoff.R
```
