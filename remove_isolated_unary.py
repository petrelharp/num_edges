#!/usr/bin/env python

import sys
import tskit, msprime
import numpy as np

usage = f"""
Usage:
    python {sys.argv[0]} (input .trees filename) (output .trees filename)

"""

def get_intervals(ts):
    last_start = np.full(ts.num_nodes, 0.0)
    is_coalescent = np.full(ts.num_nodes, False)
    last_nodes = set()
    remove_intervals = { n: [] for n in np.arange(ts.num_nodes) }
    for t in ts.trees():
        pos = t.interval.left
        these_nodes = set(t.nodes())
        added = these_nodes.difference(last_nodes)
        removed = last_nodes.difference(these_nodes)
        for n in removed:
            if not is_coalescent[n]:
                remove_intervals[n].append([last_start[n], pos])
        for n in added:
            is_coalescent[n] = False
            last_start[n] = pos
        for n in these_nodes:
            if t.num_children(n) > 1:
                is_coalescent[n] = True
        last_nodes = these_nodes
    # deal with last tree
    for n in last_nodes:
        if not is_coalescent[n]:
            remove_intervals[n].append([last_start[n], pos])
    # whoops, except for samples
    for n in ts.samples():
        remove_intervals[n] = []
    return remove_intervals

def in_interval(x, ab):
    out = False
    for left, right in ab:
        out |= x >= left and x < right
    return out

def get_node_map(ts, remove_intervals):
    """
    Find for each node the distinct intervals and remapped parents on those intervals
    """
    last_nodes = set()
    last_start = np.full(ts.num_nodes, 0.0)
    node_map = { n: [] for n in np.arange(ts.num_nodes) }
    parent = np.full(ts.num_nodes, tskit.NULL)
    to_remap = set()
    for t in ts.trees():
        pos = t.interval.left
        these_nodes = set(t.nodes())
        added = these_nodes.difference(last_nodes)
        removed = last_nodes.difference(these_nodes)
        for n in added:
            last_start[n] = pos
            if in_interval(pos, remove_intervals[n]):
                to_remap.add(n)
                parent[n] = tskit.NULL
        for n in to_remap:
            u = n
            while u != tskit.NULL and u in to_remap:
                u = t.parent(u)
            if u != parent[n]:
                # parent changed!
                if last_start[n] < pos:
                    node_map[n].append(([float(last_start[n]), pos], int(parent[n])))
                    last_start[n] = pos
                parent[n] = u
        to_remap.difference_update(removed)
        last_nodes = these_nodes
    # last tree
    pos = ts.sequence_length
    for n in to_remap:
        if parent[n] != tskit.NULL and last_start[n] < pos:
            # parent changed!
            node_map[n].append(([float(last_start[n]), pos], int(parent[n])))
    return node_map

def check_node_map(ts, node_map):
    for n in node_map:
        for (a, b), x in node_map[n]:
            assert a < b, f"{n}: ({a}, {b}) -> {x}"
            assert ts.nodes_time[x] > ts.nodes_time[n], f"{n} (t={ts.nodes_time[n]}): ({a}, {b}) -> {x} (t={ts.nodes_time[x]})"

def overlaps(a, b):
    return b[1] > a[0] and a[1] > b[0]

def overlaps_any(a, bb):
    out = False
    for b in bb:
        out |= overlaps(a, b[0])
    return out

def overlaps_map(a, bb):
    for b in bb:
        if overlaps(a, b[0]):
            yield b


def remove_isolated_unary(ts, debug=False):
    remove_intervals = get_intervals(ts)
    node_map = get_node_map(ts, remove_intervals)
    if debug:
        check_node_map(ts, node_map)
    tables = ts.dump_tables()
    edges = tables.edges
    edges.clear()
    for e in ts.edges():
        if not overlaps_any((e.left, e.right), node_map[e.child]):
            remapped = False
            for (left, right), new_parent in overlaps_map((e.left, e.right), node_map[e.parent]):
                remapped = True
                assert max(e.left, left) < min(e.right, right), f"max({e.left}, {left}) < min({e.right}, {right}) "
                if new_parent != tskit.NULL:
                    edges.add_row(parent=new_parent, child=e.child, left=max(e.left, left), right=min(e.right, right))
            if not remapped:
                edges.append(e)
    tables.sort()
    return tables.tree_sequence()


if False: # for testing
    for seed in range(101, 120):
        print("seed", seed)
        ts = msprime.sim_ancestry(2, sequence_length=1e4, recombination_rate=1e-8, coalescing_segments_only=False, population_size=1e4, model='dtwf', random_seed=seed)
        new_ts = remove_isolated_unary(ts, debug=True)

        print(ts.draw_text())
        print(new_ts.draw_text())

        for _, t1, t2 in ts.coiterate(new_ts):
            for a in ts.samples():
                for b in ts.samples():
                    assert t1.mrca(a, b) == t2.mrca(a, b), f"{a}, {b}: {t1.tmrca(a, b)} != {t2.tmrca(a, b)}"

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(usage)
        exit(0)
    infile = sys.argv[1]
    outfile = sys.argv[2]
    ts = tskit.load(infile)
    rts = remove_isolated_unary(ts)
    rts.dump(outfile)
