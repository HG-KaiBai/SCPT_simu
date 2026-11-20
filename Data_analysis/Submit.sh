#!/usr/bin/env bash
# Submit both datasets, or one: bash Submit.sh tandmobALL
set -e
analysis_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
if [ "$#" -eq 0 ]; then set -- tandmobALL HIV; fi
for dataset in "$@"; do
    mkdir -p "$analysis_dir/results/$dataset"
    sbatch --job-name="empt-$dataset" --chdir="$analysis_dir" \
        --output="$analysis_dir/results/$dataset/slurm.log" \
        --error="$analysis_dir/results/$dataset/slurm.log" \
        "$analysis_dir/Run.sbatch" "$dataset"
done
