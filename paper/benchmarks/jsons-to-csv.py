#!/usr/bin/env python
import sys, os
import json
import glob
import pandas as pd

assert len(sys.argv) == 3, "Usage: json-to-csv.py <directory name> <output csv name>"
basedir = sys.argv[1]
outfile = sys.argv[2]

jj = []
for jf in glob.glob(os.path.join(basedir, "*.json")):
    with open(jf, "r") as f:
        jj.append(json.load(f))

df = pd.DataFrame(jj)

df.to_csv(outfile)
