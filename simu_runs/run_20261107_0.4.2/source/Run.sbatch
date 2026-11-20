#!/usr/bin/env bash
#SBATCH --account=def-anthree
#SBATCH --job-name=empt-simu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --requeue
set -euo pipefail
run_dir=${1:?Prepared run directory required}
mode=${2:-worker}
source "$run_dir/source/Environment.sh"
args=()
if [[ $mode == worker ]]; then
    args=(--task="${SLURM_ARRAY_TASK_ID:?}" --cores="${SLURM_CPUS_PER_TASK:?}")
fi
exec Rscript "$run_dir/source/Pipeline.R" "$mode" --run-dir="$run_dir" "${args[@]}"
