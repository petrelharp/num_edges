#!/bin/bash

set -eo pipefail

OUTDIR="one_pop_results"
CSVFILE=${OUTDIR}.csv
if [ ! -e $OUTDIR ]
then
    mkdir -p $OUTDIR

    for L in 1e6 5e6 1e7 5e7 1e8
    do
        for N in 100 500 1000 5000 10000
        do
            SEED=$RANDOM$RANDOM
            OUT=$OUTDIR/run_${L}_${N}_${SEED}
            echo "$OUT  ....."
            python run_experiment.py --num_samples $N --length $L \
                --seed $SEED one_pop.trees $OUT
            SEED=$((SEED + 1))
        done
    done
else
    echo "${OUTDIR} already exists; doing nothing."
fi

python jsons-to-csv.py $OUTDIR $CSVFILE
python plot_results.py $CSVFILE
