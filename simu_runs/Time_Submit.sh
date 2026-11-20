#!/usr/bin/env bash
set -euo pipefail
simu_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
source "$simu_dir/Environment.sh"
prepare_only=0
if [[ ${1:-} == --prepare-only ]]; then prepare_only=1; shift; fi
if [[ $# != 1 ]]; then
    echo 'Usage: Time_Submit.sh [--prepare-only] --run-dir=PARENT_RUN' >&2
    echo '       Time_Submit.sh --submit-time=PARENT_RUN/time_evaluation' >&2
    exit 1
fi
case $1 in
    --run-dir=*)
        parent=${1#--run-dir=}
        Rscript "$simu_dir/Time_Evaluation.R" prepare --run-dir="$parent" \
            --replicates="${EMPT_TIME_REPLICATES:-50}"
        time_dir="$parent/time_evaluation"
        ;;
    --submit-time=*)
        [[ $prepare_only == 0 ]] || { echo 'Cannot combine --prepare-only and --submit-time' >&2; exit 1; }
        time_dir=${1#--submit-time=}
        ;;
    *) echo 'Expected --run-dir=PATH or --submit-time=PATH' >&2; exit 1 ;;
esac
time_dir=$(cd -- "$time_dir" && pwd -P)
Rscript "$time_dir/source/Time_Evaluation.R" status --time-dir="$time_dir"
if (( prepare_only )); then exit 0; fi
[[ ! -e $time_dir/job_id.txt ]] || { echo 'Already submitted; inspect job_id.txt before recovery.' >&2; exit 1; }
[[ ! -e $time_dir/COMPLETE.rds ]] || { echo 'Benchmark already complete.' >&2; exit 1; }
walltime=${EMPT_TIME_WALLTIME:-1-00:00:00}
extra=()
if [[ -n ${EMPT_TIME_PARTITION:-} ]]; then extra+=(--partition="$EMPT_TIME_PARTITION"); fi
printf 'nodes=1\nexclusive=yes\nbenchmark_processes=1\nbenchmark_cpus=1\nthreads=1\nwalltime=%s\nmemory=16G\n' \
    "$walltime" > "$time_dir/resources.txt"
job=$(sbatch --parsable --exclusive --nodes=1 --ntasks=1 --cpus-per-task=1 \
    --time="$walltime" "${extra[@]}" \
    --output="$time_dir/logs/time_%j.out" --error="$time_dir/logs/time_%j.err" \
    "$time_dir/source/Time_Run.sbatch" "$time_dir")
job=${job%%;*}
printf '%s\n' "$job" > "$time_dir/job_id.txt"
printf 'SUBMITTED\nSlurm job %s; measured results will be written after execution.\n' "$job" > "$time_dir/STATUS.txt"
printf 'Timing job: %s\nResults: %s\n' "$job" "$time_dir"
