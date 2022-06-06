import numpy as np
import tskit

# from https://github.com/tskit-dev/tskit/issues/2118

def mirror_coordinates(ts):
    """
    Returns a copy of the specified tables in which all
    coordinates x are transformed into L - x.
    """
    L = ts.sequence_length
    tables = ts.dump_tables()
    left = tables.edges.left
    right = tables.edges.right
    tables.edges.left = L - right
    tables.edges.right = L - left
    tables.sites.position = L - tables.sites.position - 1
    tables.sort()
    return tables.tree_sequence()

def forward_extend(ts, verbose=False, return_tables=False):
    edgediffs = ts.edge_diffs()
    # skip first tree
    num_edges = np.full(ts.num_nodes, 0)
    (_, edges_out, edges_in) = next(edgediffs)
    assert len(edges_out) == 0
    for e in edges_in:
        num_edges[e.parent] += 1
        num_edges[e.child] += 1


    t = ts.tables
    edges = t.edges
    new_left = edges.left
    new_right = edges.right

    pending_out = []
    pending_in = []
    for (interval, edges_out, edges_in) in edgediffs:
        edges_out.extend(pending_out)
        edges_in.extend(pending_in)
        # first update number of referring edges to what will happen after these are removed/added
        for e in edges_out:
            num_edges[e.parent] -= 1
            num_edges[e.child] -=1
        for e in edges_in:
            num_edges[e.parent] += 1
            num_edges[e.child] +=1
        assert np.all(num_edges >= 0)
        pending_out = []
        pending_in = []
        extended = [False for _ in edges_out]
        # remove those edges we've removed entirely from consideration
        for e1 in pending_out:
            for e2 in pending_in:
                if e1.parent == e2.parent and e1.child == e2.child:
                    pending_out.remove(e1)
                    pending_in.remove(e2)
        for j1, e1 in enumerate(edges_out):
            if not extended[j1]:
                for j2, e2 in enumerate(edges_out):
                    if not extended[j2]:
                        # need the intermediate node to not be present in the new tree
                        if (e1.parent == e2.child) and (num_edges[e2.child] == 0):
                            for e_in in edges_in:
                                if e_in.right > interval.left:
                                    if e1.child == e_in.child and e2.parent == e_in.parent:
                                        # extend e2->e1 and postpone e_in
                                        # print("     ", interval, e1.child, "-", e1.parent, "-", e2.parent, " -> ", e_in)
                                        if verbose:
                                            print("ping!!", e1.id, "+", e2.id, "=", e_in.id)
                                            print("      ", f"{e1.child}:{e1.parent} + "
                                                  f"{e2.child}:{e2.parent} -> {e_in.child}:{e_in.parent}")
                                        # extend e1 and e2, postpone e_in
                                        extended[j1] = True
                                        extended[j2] = True
                                        pending_out.extend([e1, e2])
                                        pending_in.append(e_in)
                                        new_right[e1.id] = interval.right
                                        new_right[e2.id] = interval.right
                                        new_left[e_in.id] = interval.right
                                        # amend num_edges: the intermediate node has 2 edges instead of 0
                                        num_edges[e1.parent] += 2
    keep = new_left < new_right
    edges.set_columns(
        left=new_left[keep],
        right=new_right[keep],
        parent=edges.parent[keep],
        child=edges.child[keep])
    t.build_index()
    if return_tables:
        return t
    else:
        return t.tree_sequence()


def extend_edges(ts, max_iter=100, verbose=False):
    num_edges = [ts.num_edges]
    for _ in range(max_iter):
        if verbose:
            print("Forwards")
        ts = forward_extend(ts, verbose=verbose)
        ts = mirror_coordinates(ts)
        if verbose:
            print("Reverse")
        ts = forward_extend(ts, verbose=verbose)
        ts = mirror_coordinates(ts)
        if ts.num_edges == num_edges[-1]:
            break
        else:
            num_edges.append(ts.num_edges)
    return ts, num_edges
