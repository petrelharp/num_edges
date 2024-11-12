import sys, os
import json
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

assert len(sys.argv) == 2, "Usage: python plot_results.py <csv file>"
infile = sys.argv[1]

df = pd.read_csv(infile)
basename = ".".join(infile.split(".")[:-1])

# average across replicates:
sdf = df.groupby(["length", "num_samples"]).agg({
        "length": "mean",
        "num_samples": "mean",
        "num_edges_before": "mean",
        "num_edges_after": "mean",
        "runtime_before": "mean",
        "runtime_after": "mean",
        "extend_time": "mean"
}).astype({'length': int, 'num_samples': int})
sdf.reset_index(drop=True, inplace=True)

# from Halley
colors = {'blue': (46,37,133),
          'red': (194,106,119),
          'lgreen': (93,168,153),
          'gold': (220,205,125),
          'green': (51, 117,56),
          'lblue': (148,203,236),
          'magenta': (159,74,150),
          'wine': (126,41,84),
         }
for x in colors:
    colors[x] = [u/256 for u in colors[x]] + [1.0]

length_vals = np.unique(df.length)
# length_cols = plt.colormaps["viridis"](np.linspace(0, 1, len(length_vals)))
length_cols = np.array([colors[n] for n in ['blue', 'green', 'lblue', 'lgreen', 'gold']])
ns_vals = np.unique(df.num_samples)
# ns_cols = plt.colormaps["plasma"](np.linspace(0, 1, len(ns_vals)))
ns_cols = np.array([colors[n] for n in ['wine', 'magenta', 'red', 'gold', 'blue']])

# combined number of edges & runtime for Tajimas D
fname = f"{basename}_absolute_values.pdf"
fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(6.5, 3))

ax0.set_xlabel("sequence length")
ax1.set_xlabel("sequence length")
ax0.set_ylabel("number of edges")
ax1.set_ylabel("runtime, Tajima's D")
for ns, this_df in df.groupby("num_samples"):
    sub_df = this_df.sort_values("length")
    j,  = np.where(ns == ns_vals)
    ax0.scatter(sub_df.length, sub_df.num_edges_before, marker=".",
             c=ns_cols[j])
    ax0.scatter(sub_df.length, sub_df.num_edges_after, marker="x",
             c=ns_cols[j])
    ax1.scatter(sub_df.length, sub_df.runtime_before, marker=".",
             c=ns_cols[j])
    ax1.scatter(sub_df.length, sub_df.runtime_after, marker="x",
             c=ns_cols[j])

for ns, this_df in sdf.groupby("num_samples"):
    sub_df = this_df.sort_values("length")
    j,  = np.where(ns == ns_vals)
    ax0.plot(sub_df.length, sub_df.num_edges_before,
             c=ns_cols[j])
    ax0.plot(sub_df.length, sub_df.num_edges_after,
             c=ns_cols[j])
    ax1.plot(sub_df.length, sub_df.runtime_before, label=f"samples={ns}",
             c=ns_cols[j])
    ax1.plot(sub_df.length, sub_df.runtime_after,
             c=ns_cols[j])

ax1.legend(fontsize="x-small")

plt.tight_layout()
plt.savefig(fname)


# combined number of edges & runtime for Tajimas D
fname = f"{basename}_ratios.pdf"
fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(6.5, 3))

ax0.set_xlabel("sequence length")
ax1.set_xlabel("sequence length")
ax0.set_ylabel("number of edges, after/before")
ax1.set_ylabel("runtime, after/before")
for ns, this_df in df.groupby("num_samples"):
    sub_df = this_df.sort_values("length")
    j,  = np.where(ns == ns_vals)
    ax0.scatter(sub_df.length, sub_df.num_edges_after / sub_df.num_edges_before,
             marker=".",
             c=ns_cols[j])
    ax1.scatter(sub_df.length, sub_df.runtime_after / sub_df.runtime_before,
             marker=".",
             c=ns_cols[j])

for ns, this_df in sdf.groupby("num_samples"):
    sub_df = this_df.sort_values("length")
    j,  = np.where(ns == ns_vals)
    ax0.plot(sub_df.length, sub_df.num_edges_after / sub_df.num_edges_before,
             label=f"samples={ns}",
             c=ns_cols[j])
    ax1.plot(sub_df.length, sub_df.runtime_after / sub_df.runtime_before,
             label=f"samples={ns}",
             c=ns_cols[j])

ax0.set_ylim((0, max(1, 1.05 * np.max(df.num_edges_after / df.num_edges_before))))
ax1.set_ylim((0, max(1, 1.05 * np.max(df.runtime_after / df.runtime_before))))
ax0.legend(fontsize="x-small")

plt.tight_layout()
plt.savefig(fname)


# timing plots
fname = f"{basename}_timing.pdf"
fig, (ax0, ax1, ax2) = plt.subplots(1, 3, figsize=(6.5, 3), sharey=True)

ax0.set_xlabel("number of samples")
ax0.set_ylabel("runtime (s)")
ax1.set_xlabel("sequence length")
ax1.set_ylabel("runtime (s)")
ax2.set_xlabel("number of edges")
ax2.set_ylabel("runtime (s)")
ax0.set_yscale("log")
ax1.set_yscale("log")
ax2.set_yscale("log")
for L, this_df in df.groupby("length"):
    sub_df = this_df.sort_values("num_samples")
    j,  = np.where(L == length_vals)
    ax0.scatter(sub_df.num_samples, sub_df.extend_time, marker="o",
             c=length_cols[j])

for ns, this_df in df.groupby("num_samples"):
    sub_df = this_df.sort_values("length")
    j,  = np.where(ns == ns_vals)
    ax1.scatter(sub_df.length, sub_df.extend_time, marker="o",
             c=ns_cols[j])

for L, this_df in df.groupby("length"):
    sub_df = this_df.sort_values("num_samples")
    j,  = np.where(L == length_vals)
    ax2.scatter(sub_df.num_edges_before, sub_df.extend_time, marker="o", label=f"length={L:g}",
                c=length_cols[j])

for L, this_df in sdf.groupby("length"):
    sub_df = this_df.sort_values("num_samples")
    j,  = np.where(L == length_vals)
    ax0.plot(sub_df.num_samples, sub_df.extend_time, label=f"length={L:g}",
             c=length_cols[j])

for ns, this_df in sdf.groupby("num_samples"):
    sub_df = this_df.sort_values("length")
    j,  = np.where(ns == ns_vals)
    ax1.plot(sub_df.length, sub_df.extend_time, label=f"samples={ns}",
             c=ns_cols[j])

ax0.legend(fontsize='xx-small')
ax1.legend(fontsize='xx-small')
ax2.legend(fontsize='xx-small')

plt.tight_layout()
plt.savefig(fname)
