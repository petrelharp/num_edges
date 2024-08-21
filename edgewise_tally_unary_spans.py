import msprime
import tskit
import numpy as np
from itertools import groupby

def edgewise_tally_unary_spans(ts):
    """
    Tally spans for each node into three classes, corresponding to the columns of the output:
        0. Region is not part of a contiguous span in which it is coalescent somewhere.
        1. Region is part of a contiguous span in which it is coalescent somewhere,
            and is unary.
        2. Region is coalescent.
    """
    out = np.zeros((ts.num_nodes, 3))
    p = ts.edge(0).parent
    pi = 0
    for e in ts.edges():
        if e.parent != p or e.id + 1 == ts.num_edges:
            # done with parent p
            pj = e.id + (e.id + 1 == ts.num_edges)
            lefts = ts.edges_left[pi:pj]
            rights = ts.edges_right[pi:pj]
            breaks = np.sort(np.unique(np.concatenate([[0], lefts, rights])))
            overlaps = np.zeros(len(breaks))
            for j, x in enumerate(breaks):
                overlaps[j] = np.sum(lefts <= x) - np.sum(rights <= x)
            not_isolated = np.full(len(breaks) - 1, False)
            coalescent_region = False
            for j in range(0, len(breaks) - 1, 1):
                if overlaps[j] == 0:
                    coalescent_region = False
                else:
                    if overlaps[j] == 1:
                        not_isolated[j] = coalescent_region
                    else:
                        coalescent_region = True
            coalescent_region = False
            for j in range(len(breaks) - 2, -1, -1):
                if overlaps[j] == 0:
                    coalescent_region = False
                else:
                    if overlaps[j] == 1:
                        not_isolated[j] = not_isolated[j] or coalescent_region
                    else:
                        coalescent_region = True
            for j in range(0, len(breaks) - 1, 1):
                if overlaps[j] > 1:
                    out[p, 2] += breaks[j + 1] - breaks[j]
                elif overlaps[j] == 1:
                    if not_isolated[j]:
                        out[p, 1] += breaks[j + 1] - breaks[j]
                    else:
                        out[p, 0] += breaks[j + 1] - breaks[j]
            # on to the next p
            p = e.parent
            pi = e.id
    return out

ts = msprime.sim_ancestry(10, sequence_length=1e7, recombination_rate=1e-8, population_size=1e4, coalescing_segments_only=False, random_seed=1)
ets = ts.simplify().extend_paths()

unary_spans = []
for x in [ts, ets]:
    out = edgewise_tally_unary_spans(x)
    sums = out.sum(axis=0)
    total = np.sum(sums[:2])
    props = 100 * np.round(sums / total, 2)
    print(f"Total unary span bordering coalescing: {sums[1]} ({props[1]}%)")
    print(f"Total trapped unary span: {sums[0]} ({props[0]}%)")
    print(f"Total coalescing span: {sums[2]}")
    unary_spans.append(out)

import matplotlib.pyplot as plt
out0, out1 = unary_spans
order = np.argsort(out0[:, 1])

plt.scatter(np.arange(order.size), np.log10(out1[order, 1]), label="extend-edges", s=4, alpha=0.1)
plt.plot(np.arange(order.size), np.log10(out0[order, 1]), label="simulated", c='black')
plt.xlabel("Rank")
plt.ylabel("Unary span (log10)")
plt.legend()
plt.savefig("ee_test_log.png")

plt.scatter(np.arange(order.size), out1[order, 1], label="extend-edges", s=4, alpha=0.1)
plt.plot(np.arange(order.size), out0[order, 1], label="simulated", c='black')
plt.xlabel("Rank")
plt.ylabel("Unary span (log10)")
plt.legend()
plt.savefig("ee_test.png")
