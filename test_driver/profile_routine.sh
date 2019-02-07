#!/bin/bash

# Binary options: cudenseTest cudenseTest_static
BINARY="./test_zhegvdx"
M_DIMENSION="800 1000 2000 3000 4000 5000 6000 7000"
ITERS="1 2 3 4 5"
PROFILER="/usr/local/cuda-10.0/bin/nvprof"

export MKL_NUM_THREADS="4"
export MKL_DYNAMIC="FALSE"
export OMP_NUM_THREADS=1

for M in ${M_DIMENSION}
do
    for i in ${ITERS}
    do
        #COMMAND_LINE="${PROFILER} --log-file nvtx_files/nvtx_${ROUTINE}_${M}.csv --csv ${BINARY} ${M}"
        COMMAND_LINE="./test_zhegvdx ${M}"
#        echo "Running ${COMMAND_LINE}"
        RESULT="$(${COMMAND_LINE})"
        echo $RESULT
    done
    echo ""
done
