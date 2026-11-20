#!/usr/bin/env bash
#SBATCH --account=def-anthree
#SBATCH --job-name=empt-data
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --no-requeue
set -e
source ../simu_runs/Environment.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 RCPP_PARALLEL_NUM_THREADS=1 OMP_DYNAMIC=FALSE MKL_DYNAMIC=FALSE
export TZ=America/New_York
srun --exact --nodes=1 --ntasks=1 --cpus-per-task=1 --cpu-bind=threads Rscript "$1.R"
