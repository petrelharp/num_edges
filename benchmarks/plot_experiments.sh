#!/bin/bash

set -eo pipefail

for OUTBASE in one_pop constant_pop
do
    OUTDIR="${OUTBASE}_results"
    CSVFILE=${OUTDIR}.csv
    if [ ! -e $CSVFILE ]
    then
        python jsons-to-csv.py $OUTDIR $CSVFILE
    else
        echo "$CSVFILE already exists, not regenerating."
    fi
    python plot_results.py $CSVFILE
done
