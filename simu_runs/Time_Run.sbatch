#!/usr/bin/env bash
#SBATCH --account=def-anthree
#SBATCH --job-name=empt-time
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --hint=nomultithread
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --no-requeue
set -euo pipefail
time_dir=${1:?Prepared time_evaluation directory required}
source "$time_dir/source/Environment.sh"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 RCPP_PARALLEL_NUM_THREADS=1 OMP_DYNAMIC=FALSE MKL_DYNAMIC=FALSE
export TZ=America/New_York
allocation="$time_dir/logs/allocation_${SLURM_JOB_ID:?}.txt"
scontrol show job "$SLURM_JOB_ID" > "$allocation"
# Fail closed if launched without node exclusivity.
if ! grep -Eq 'OverSubscribe=NO|Shared=0' "$allocation"; then
    echo 'Timing requires an exclusive Slurm node allocation.' >&2
    exit 1
fi
export EMPT_TIME_EXCLUSIVE=1
scontrol show node "$SLURMD_NODENAME" > "$time_dir/logs/node_${SLURM_JOB_ID}.txt"
lscpu > "$time_dir/logs/lscpu_${SLURM_JOB_ID}.txt"
trap 'rc=$?; if (( rc != 0 )); then printf "Job %s failed with exit %s; inspect logs and retained checkpoints.\n" "$SLURM_JOB_ID" "$rc" > "$time_dir/JOB_FAILED.txt"; fi' EXIT
# --exact prevents the exclusive allocation from expanding the benchmark step.
srun --exact --nodes=1 --ntasks=1 --cpus-per-task=1 --cpu-bind=threads \
    Rscript "$time_dir/source/Verify_Time_Evaluation.R" --real
srun --exact --nodes=1 --ntasks=1 --cpus-per-task=1 --cpu-bind=threads \
    Rscript "$time_dir/source/Time_Evaluation.R" worker --time-dir="$time_dir"
