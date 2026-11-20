# Unified EMPT simulation — SurvPT 0.4.2

## Complete run

From the project directory:

```bash
bash Submit_Simulation.sh --prepare-only
bash Submit_Simulation.sh
```

Each invocation creates a new run directory. To submit an already prepared run:

```bash
bash simu_runs/Submit.sh --submit-run=simu_runs/run_ID
```

The default complete run includes 64 settings x 1000 replications (640 tasks at
100 replications/task). It produces point estimates, posterior intervals,
kappa sensitivity, gamma sensitivity and failure reports in one submission.
The merge job starts after every worker succeeds. No supplement or backfill is needed.

## Scientific specification

- Fiducial: scenarios 1–4, n=50,75,100,200,500,800,1000 (28 settings).
- Exponential/Weibull: shape=1/2, inverse scale=.12, n=200,500, q=.3,.5,
  right/current/interval censoring (24 settings).
- Mixture: fixed equal mixture of Weibull(shape=3,scale=3) and
  Weibull(shape=4,scale=10), n=200,500, q=.3,.5 and the same censoring grid
  (12 settings). `SurvPT::data_mixture` defines q as the right-open interval
  probability. Mixture settings remain IDs 53–64.
- Main methods: Unif, Unif_Flat, Exp, Exp_Flat, Weibull_PT, Weibull_PT_Flat,
  Fiducial, NPMLE, MPT, TBP; OPT additionally for right censoring.
- PT: max_M=8, max_EM_iter=2000, EM_tol=1e-6, main kappa=.1, gamma=.5.
  Weibull priors use censored-data `SurvPT::estimate(L,R)` and fix its parameters
  during fitting and posterior sampling. Other parameters retain package defaults.
- Kappa: all six PT methods x .05,.1,.2,.3 (24 variants). Gamma: the three
  non-flat PT methods x .5,1,2 at kappa=.1 (nine variants).
- PT posterior: 2000 retained draws, 3000 burn-in, thin=1; point predictions
  explicitly use `sample=FALSE`. MPT/TBP retain spBayesSurv 1.1.9 defaults,
  including the existing baseline adapters. Fiducial retains mfid=1000, mburn=100.
- Point grid: 201 equally spaced times from 0 to the largest finite endpoint.
  Posterior targets: the three true times where S(t)=.75,.5,.25. They are not
  truncated at the observed maximum. Interval level=.95.
- Fiducial point curves use the original interpolated sample median; posterior
  tables use the fiducial sample mean and equal-tailed fiducial interval.

Every method uses the same dataset within a replication. Data seed=replication;
method seed=100000000+1000000*method_id+10000*setting_id+replication. Historical
IDs are retained. Kappa=.1 and gamma=.5 reuse the matching main fit and posterior.
Chunk size, CPU count and execution order do not change method seeds.

For reproducibility, keep the same R/package versions, installed SurvPT archive
hash and scientific settings. Each run freezes these plus source hashes and RNGkind.
Numerical results and failure locations are reproducible under that environment;
runtime, run IDs, timestamps and manifest hashes naturally differ across runs.

## Error handling and summaries

SurvPT (all six priors/flat variants and every kappa/gamma variant), OPT,
MPT and TBP errors in fitting, prediction or result validation are isolated to that
method and replication. The original condition, call, method seed, tau and affected
targets are recorded. Warnings alone do not skip a fit. NPMLE/Fiducial errors still
stop a task after successful replication peers are checkpointed.

For PT nonfinite posterior targets, valid point estimates and other posterior
targets are retained (`partial_error`). Beyond-tau NA is reported as
`prediction_outside_tau`; other NA/NaN/Inf as `non_finite_posterior`. A failed fit
or invalid point curve becomes `skipped_error`, with missing results for that method.
No replacement data, changed seeds, tail extrapolation or retries alter the analysis.

Means, SDs and coverage exclude missing values. Point summary `replicates` counts
valid point curves. Posterior `replicates` counts valid posterior summaries separately
for each quantile, while `estimate_replicates` counts valid estimates there. All-failed
metrics are NA with count 0. A single success has SD=NA. Coverage NA is never treated
as a missed interval. Known target times remain available.

Every posterior CSV contains `S75_estimate`, `S75_mean`, `S75_estimate_absolute_bias`,
`S75_absolute_bias` and corresponding S50/S25 fields, plus intervals and scores.
`estimate` is PT's exact MAP-EM prediction at the true quantile (not interpolation
from the 201-point grid), the original Fiducial median, or MPT/TBP's Shat point
prediction (which is their posterior mean). `absolute_bias` retains its existing
posterior-mean definition; `estimate_absolute_bias=abs(estimate-p)` is separate.
Both are averaged across replications, with independent effective sample counts.

- SSE=sum((S_hat-S_true)^2).
- RMSE compares probability increments after appending zero to both survival curves.
  RMSE_nolast omits that appended zero. These are not survival-value RMSEs.
- Two absolute biases are abs(estimate-p) and abs(posterior_mean-p), each averaged per replication.
- Winkler score=(U-L)+40*(max(L-p,0)+max(p-U,0)) for 95% intervals,
  averaged per replication. Coverage is the proportion of successful intervals
  covering p; a failed fit is missing rather than an uncovered interval.

## Outputs

- `point_estimate/<setting>/replicates.csv`, `summary.csv`, and the global summary.
- `posterior/<setting>/<method>.csv`: one row per replication plus Average.
  Fields include target time, estimate, posterior mean, interval endpoints, both absolute biases, length,
  coverage and Winkler score. `posterior/summary.csv` includes effective counts.
- `sensitivity/` and `gamma_sensitivity/` contain the corresponding point/posterior outputs.
- `datasets/summary.csv`: one row per setting, 1000 datasets, `average_tau`,
  `true_S_at_average_tau`, `average_Exact`, `average_Left`, `average_Interval`,
  `average_Right` and `average_total_censored`. Tau and the four proportions come
  directly from the package attributes. True S is independently evaluated at the
  average tau; this is not the average of `attr(d, "true_S_tau")`.
- `datasets/replicates.csv`: the dataset attributes for each setting/seed. A separate
  dataset job produces this report early; merge recomputes it from the actual fit
  checkpoints, including datasets where methods failed.
- `fit_failures.csv`: every failed or partially missing method/replication, original
  error, category, stage, affected targets, tau, canonical method and seed.
- `method_errors/task_NNNN.csv`: errors exported after each worker checkpoint wave.
- `fit_failure_summary.csv`: failed-fit counts by method and error category.
- `fit_status_summary.csv`: expected, attempted, successful, failed and unattempted
  counts for each setting/group/method, skipped/partial counts, and separate point,
  estimate and posterior target success counts. Repeated merges do not double-count.
- `source/`, `packages/`, `manifest.rds`, design CSVs, `chunks/`, `logs/`, and
  `COMPLETE.rds` retain reproducibility and completion records.

## Controls and recovery

Defaults: 16 CPUs/task, 32G, 12 hours, no array concurrency throttle; merge uses
1 CPU, 8G and 12 hours. Override `EMPT_CPUS_PER_TASK`, `EMPT_MEMORY`,
`EMPT_WALLTIME`, `EMPT_CONCURRENT_TASKS`, `EMPT_CHUNK_SIZE`, `EMPT_REPLICATES`,
`EMPT_RUN_DIR`, `EMPT_SUITE=all|exp|weibull|fiducial|mixture`,
`EMPT_ANALYSIS=all|main|sensitivity|gamma`, or `EMPT_AFTER_JOB` as needed.
A full run uses suite=all and analysis=all. A small pilot is:

```bash
EMPT_REPLICATES=2 bash simu_runs/Submit.sh --settings=1,29,53
```

Workers checkpoint each wave and reuse validated saved replications on restart.
Always use the run's frozen source for status, workers and merge:

```bash
source simu_runs/Environment.sh
Rscript simu_runs/run_ID/source/Pipeline.R status --run-dir=simu_runs/run_ID
Rscript simu_runs/run_ID/source/Pipeline.R worker --run-dir=simu_runs/run_ID --task=N --cores=16
Rscript simu_runs/run_ID/source/Pipeline.R merge --run-dir=simu_runs/run_ID
Rscript simu_runs/run_ID/source/Pipeline.R errors --run-dir=simu_runs/run_ID
Rscript simu_runs/run_ID/source/Pipeline.R datasets --run-dir=simu_runs/run_ID
```

A task lock prevents concurrent writers. After a killed worker, confirm it has
stopped before removing only its stale `chunks/task_NNNN.rds.lock` directory.
Historical sources and results remain in their run directories. The already
submitted 2026-11-04 mixture recovery uses
`run_20261104_0.4.1/mixture_gamma_supplement/method_error_skip/source/Resume_Method_Errors.R`;
its frozen recovery code continues to support worker, merge, check and report.

## Validation

```bash
source simu_runs/Environment.sh
Rscript simu_runs/Verify_Method_Failures.R /tmp/empt-error-check-NEW
Rscript simu_runs/Verify_SurvPT_Errors.R /tmp/empt-survpt-check-NEW
Rscript simu_runs/Verify_Datasets.R
Rscript simu_runs/Verify_Posterior.R
Rscript simu_runs/Verify_Fiducial_Posterior.R
Rscript simu_runs/Verify.R /tmp/empt-full-check-NEW
```

The checks cover default 64-setting preparation, injected method/validation errors, real PT tail NA,
method isolation, NA rows, summary denominators, all-failed cases, serial/parallel
and chunk partition parity, checkpoint corruption, repeated merge, posterior
metrics, fiducial sample transformations and real fits with default sampling.


## Probability roundoff policy (2026-11-08)

Posterior `estimate`, `mean`, `lower` and `upper` permit finite boundary deviations
up to 1e-8. Only values in [-1e-8, 0) or (1, 1+1e-8] are clipped to the nearest
probability boundary. Larger deviations remain errors; the existing missing-target
policy is unchanged. Valid values in [0,1] are never changed. Biases, interval
lengths, coverage and Winkler scores use the corrected posterior probabilities.
The point-grid values retain their original 1e-8 validation and are not clipped.
This policy is frozen in every new manifest (version 6), along with source hashes.
No fitting parameters, posterior draws, data seeds or method seeds change.

Completed runs use the standard `source/`, `chunks/`, `manifest.rds` and
`COMPLETE.rds` layout. Use the frozen `source/Pipeline.R` entry point above for
status, workers, reports and merge. No supplemental entry point is required.

For a new full study, use `bash Submit_Simulation.sh`. The same probability policy
and seeds are frozen automatically. Under the same R/package/build/RNG environment,
scientific results are reproducible; elapsed time and timestamps can differ.

`Verify_Roundoff.R` checks tolerance boundaries, rejection beyond the tolerance,
unchanged valid values/RNG state and missing targets.
