import os, time
import tskit

orig = [x for x in os.listdir("trees2") if x[:3] == "ts_"]

def run(ts):
    before = time.time()
    for _ in range(40):
        _ = ts.Tajimas_D()
    after = time.time()
    return after - before

edges = {x : {} for x in orig}
runtime = {x : {} for x in orig}
colnames = ['', 'file', 'orig', 'I', 'IE', 'IS', 'ISE']
with open("table_s3.tsv", "w") as f:
    f.write("\t".join(colnames) + "\n")
    f.flush()
    for x in orig:
        ex = edges[x]
        rx = runtime[x]
        ex['file'] = rx['file'] = x
        print(x)
        # original
        ts = tskit.load("trees2/" + x)
        ex['orig'] = ts.num_edges
        rx['orig'] = run(ts)
        # I
        ts = tskit.load("trees2/infer_" + x)
        ex['I'] = ts.num_edges
        rx['I'] = run(ts)
        # IE
        ets = ts.extend_haplotypes()
        ex['IE'] = ets.num_edges
        rx['IE'] = run(ets)
        # IS
        sts = ts.simplify()
        ex['IS'] = sts.num_edges
        rx['IS'] = run(sts)
        # ISE
        ests = sts.extend_haplotypes()
        ex['ISE'] = ests.num_edges
        rx['ISE'] = run(ests)
        f.write("\t".join(['edges', x] + [str(ex[u]) for u in colnames[2:]]) + "\n")
        f.write("\t".join(['runtime', x] + [f"{rx[u]:.2f}" for u in colnames[2:]]) + "\n")
        f.flush()
