#!/bin/bash

set -eo pipefail
if [ -e "run_py.srun" ]
then
    PYTHON="sbatch run_py.srun"
else
    PYTHON="python"
fi


for OUTBASE in one_pop constant_pop
do
    OUTDIR="${OUTBASE}_results"
    mkdir -p $OUTDIR

    for L in 1e6 5e6 1e7 5e7 1e8
    do
        for N in 100 500 1000 5000 10000
        do
            SEED=$RANDOM$RANDOM
            OUT=$OUTDIR/run_${L}_${N}_${SEED}
            echo "$OUT  ....."
            $PYTHON run_experiment.py --num_samples $N --length $L \
                --seed $SEED ${OUTBASE}.trees $OUT
            SEED=$((SEED + 1))
        done
    done
done
