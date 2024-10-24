#!/bin/bash

set -eo pipefail

OUTDIR="one_pop_results"
CSVFILE=${OUTDIR}.csv
python jsons-to-csv.py $OUTDIR $CSVFILE
python plot_results.py $CSVFILE
