#!/usr/bin/env bash
# All new runs are written beneath simu_runs. --prepare-only does not submit jobs.
set -euo pipefail
simu_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
source "$simu_dir/Environment.sh"
prepare_only=0
prepared_run=
if [[ ${1:-} == --prepare-only ]]; then prepare_only=1; shift; fi
if [[ ${1:-} == --submit-run=* ]]; then
    prepared_run=${1#--submit-run=}; shift
    [[ -n $prepared_run && $# == 0 && $prepare_only == 0 ]] || {
        echo 'Use --submit-run=PATH by itself to submit an existing prepared run.' >&2; exit 1;
    }
fi
cores=${EMPT_CPUS_PER_TASK:-16}
concurrent=${EMPT_CONCURRENT_TASKS:-0}
chunk=${EMPT_CHUNK_SIZE:-100}
for value in "$cores" "$chunk"; do
    [[ $value =~ ^[1-9][0-9]*$ ]] || { echo 'CPU and chunk values must be positive integers' >&2; exit 1; }
done
[[ $concurrent =~ ^(0|[1-9][0-9]*)$ ]] || { echo 'Concurrency must be a nonnegative integer (0 = unlimited)' >&2; exit 1; }
if [[ -n $prepared_run ]]; then
    run_dir=$prepared_run
    [[ ! -e $run_dir/array_job_id.txt ]] || {
        echo 'This run was already submitted; use its recorded job IDs and frozen worker for recovery.' >&2; exit 1;
    }
    source "$run_dir/source/Environment.sh"
    Rscript "$run_dir/source/Pipeline.R" status --run-dir="$run_dir"
else
    run_dir=${EMPT_RUN_DIR:-$simu_dir/run_$(date +%Y%m%d_%H%M%S)_$$}
    Rscript "$simu_dir/Pipeline.R" prepare --run-dir="$run_dir" --chunk-size="$chunk" \
        --replicates="${EMPT_REPLICATES:-1000}" --suite="${EMPT_SUITE:-all}" \
        --analysis="${EMPT_ANALYSIS:-all}" "$@"
fi
run_dir=$(cd -- "$run_dir" && pwd -P)
read -r task_count < "$run_dir/task_count.txt"
array_spec="1-${task_count}"
max_cpus=unlimited
if (( concurrent > 0 )); then
    array_spec+="%${concurrent}"
    max_cpus=$((cores * concurrent))
fi
printf 'cpus_per_task=%s\nconcurrent_tasks=%s\nmax_active_cpus=%s\narray_tasks=%s\nwalltime=%s\nmemory=%s\nmerge_walltime=12:00:00\n' \
    "$cores" "$concurrent" "$max_cpus" "$task_count" \
    "${EMPT_WALLTIME:-12:00:00}" "${EMPT_MEMORY:-32G}" > "$run_dir/resources.txt"
printf 'Run directory: %s\nTasks: %s\n' "$run_dir" "$task_count"
if (( prepare_only )); then exit 0; fi
queued=$(squeue -h -r -u "$(id -un)" -o '%i' | wc -l)
if (( task_count + queued + 2 > 800 )); then
    echo 'Insufficient reserved submit slots; increase EMPT_CHUNK_SIZE. Prepared run retained.' >&2
    exit 1
fi
dependency=()
if [[ -n ${EMPT_AFTER_JOB:-} ]]; then
    [[ $EMPT_AFTER_JOB =~ ^[0-9]+$ ]] || { echo 'EMPT_AFTER_JOB must be a job ID' >&2; exit 1; }
    dependency=(--dependency="afterok:${EMPT_AFTER_JOB}" --kill-on-invalid-dep=yes)
fi
dataset_job=$(sbatch --parsable --cpus-per-task=1 --mem=4G --time=01:00:00 --job-name=empt-datasets \
    --output="$run_dir/logs/datasets_%j.out" --error="$run_dir/logs/datasets_%j.err" \
    "$run_dir/source/Run.sbatch" "$run_dir" datasets)
dataset_job=${dataset_job%%;*}
printf '%s\n' "$dataset_job" > "$run_dir/dataset_job_id.txt"
array=$(sbatch "${dependency[@]}" --parsable --array="$array_spec" \
    --cpus-per-task="$cores" --mem="${EMPT_MEMORY:-32G}" --time="${EMPT_WALLTIME:-12:00:00}" \
    --output="$run_dir/logs/task_%A_%a.out" --error="$run_dir/logs/task_%A_%a.err" \
    "$run_dir/source/Run.sbatch" "$run_dir" worker)
array=${array%%;*}
printf '%s\n' "$array" > "$run_dir/array_job_id.txt"
merge=$(sbatch --parsable --dependency="afterok:$array" --kill-on-invalid-dep=yes \
    --cpus-per-task=1 --mem=8G --time=12:00:00 --job-name=empt-simu-merge \
    --output="$run_dir/logs/merge_%j.out" --error="$run_dir/logs/merge_%j.err" \
    "$run_dir/source/Run.sbatch" "$run_dir" merge)
merge=${merge%%;*}
printf '%s\n' "$merge" > "$run_dir/merge_job_id.txt"
printf 'Dataset report job: %s\nArray job: %s\nMerge job: %s\nMaximum active CPUs: %s\n' "$dataset_job" "$array" "$merge" "$max_cpus"
