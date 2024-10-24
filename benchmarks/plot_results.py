import sys
import json
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys, os

assert len(sys.argv) == 2, "Usage: python plot_results.py <directory name>"
basedir = sys.argv[1]

jj = []
for jf in glob.glob(os.path.join(basedir, "*.json")):
    with open(jf, "r") as f:
        jj.append(json.load(f))

df = pd.DataFrame(jj)

# combined number of edges & runtime for Tajimas D
fname = os.path.join(basedir, "absolute_values.pdf")
fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(9, 5))

ax0.set_xlabel("number of samples")
ax1.set_xlabel("number of samples")
ax0.set_ylabel("number of edges")
ax1.set_ylabel("runtime")
for L, sub_df in df.groupby("length"):
    ax0.plot(sub_df.num_samples, sub_df.num_edges_before, marker=".", label="before")
    ax0.plot(sub_df.num_samples, sub_df.num_edges_after, marker="o", label="after")
    ax1.plot(sub_df.num_samples, sub_df.runtime_before, marker=".", label="before")
    ax1.plot(sub_df.num_samples, sub_df.runtime_after, marker=".", label="after")

ax0.legend()

plt.tight_layout()
plt.savefig(fname)


# combined number of edges & runtime for Tajimas D
fname = os.path.join(basedir, "ratios.pdf")
fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(9, 5))

ax0.set_xlabel("number of samples")
ax1.set_xlabel("number of samples")
ax0.set_ylabel("number of edges, after/before")
ax1.set_ylabel("runtime, after/before")
for L, sub_df in df.groupby("length"):
    ax0.plot(sub_df.num_samples, sub_df.num_edges_after / sub_df.num_edges_before, marker=".", label="before")
    ax1.plot(sub_df.num_samples, sub_df.runtime_after / sub_df.runtime_before, marker=".", label="before")

ax0.legend()

plt.tight_layout()
plt.savefig(fname)


# timing plots
fname = os.path.join(basedir, "timing.pdf")
fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(9, 5), sharey=True)

ax0.set_xlabel("number of samples")
ax0.set_ylabel("runtime (s)")
ax1.set_xlabel("sequence length")
ax1.set_ylabel("runtime (s)")
for L, sub_df in df.groupby("length"):
    ax0.plot(sub_df.num_samples, sub_df.extend_time, marker="o")
    ax1.plot(sub_df.length, sub_df.extend_time, marker="o")

plt.tight_layout()
plt.savefig(fname)
