import numpy as np
import tskit

def node_tree_discrepancy(x, ts1, ts2):
    # here sim[y] is the span so far over which y has the same set of descendants as x
    sim = np.zeros(ts2.num_nodes)
   
    for interval, t1, t2 in ts1.coiterate(ts2):
        s = interval.span
        Dx = set(t1.samples(x))

        if len(Dx) > 0:
            if len(Dx) == 1:
                y = list(Dx)[0]
                MRCA = list(Dx)[0] 
            else:
                MRCA = t1.mrca(*list(Dx))
                y = t2.mrca(*list(Dx))
            if set(t2.samples(y)) == Dx:
                sim[y] += s 
                y = t2.parent(y)
                if x != MRCA:
                    while y != tskit.NULL and set(t2.samples(y)) == Dx:
                        sim[y] += s
                        y = t2.parent(y)
    return ts1.sequence_length - np.max(sim)

def discrepancy(ts1, ts2):
    dis = 0
    for n in range(ts1.num_nodes):
        dis += node_tree_discrepancy(n, ts1, ts2)
    return dis