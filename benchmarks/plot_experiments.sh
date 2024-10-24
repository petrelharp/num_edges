#!/bin/bash

set -eo pipefail

for OUTBASE in one_pop constant_pop
do
    OUTDIR="${OUTBASE}_results"
    CSVFILE=${OUTDIR}.csv
    python jsons-to-csv.py $OUTDIR $CSVFILE
    python plot_results.py $CSVFILE
done
