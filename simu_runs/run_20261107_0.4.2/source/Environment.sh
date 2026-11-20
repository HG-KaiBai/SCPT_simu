#!/usr/bin/env bash
if type module >/dev/null 2>&1; then
    module load StdEnv/2023 r/4.4.0
fi
export R_LIBS_USER=/home/y679zhao/R/empt-survpt-0.4.2:/home/y679zhao/R/empt-spbayessurv-20261025:/home/y679zhao/R/x86_64-pc-linux-gnu-library/4.4
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
