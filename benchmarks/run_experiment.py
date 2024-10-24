import argparse
import sys
import json

import tskit
import numpy as np

import time

def parse_args():
    dstring = "Run extend haplotypes, comparing numbers of edges and runtime for Tajima's D.",
    parser = argparse.ArgumentParser(description=dstring,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input', type=str,
                        help="Tree sequence filename")
    parser.add_argument('output', type=str,
                        help="Output file base name")
    parser.add_argument('--num_samples', '-N', type=int, required=True,
                        help="Number of samples")
    parser.add_argument('--length', '-L', type=float,
                        help="Genome length")
    parser.add_argument('--seed', '-S', type=int, required=True,
                        help="Random seed")
    return parser

if __name__ == "__main__":
    parser = parse_args()
    args = parser.parse_args(sys.argv[1:])

    num_samples = args.num_samples

    rng = np.random.default_rng(seed=args.seed)

    ts = tskit.load(args.input).trim()

    L = args.length
    if L is None:
        L = ts.sequence_length

    if L < ts.sequence_length:
        left = rng.choice(int(ts.sequence_length - L))
        ts = ts.keep_intervals([[left, left + L]]).trim()

    samples = rng.choice(list(ts.samples()), size=num_samples, replace=False)
    ts = ts.simplify(samples)
    num_edges_before = ts.num_edges

    start_time = time.time()
    _ = ts.Tajimas_D(mode="branch")
    runtime_before = time.time() - start_time

    start_time = time.time()
    ets = ts.extend_haplotypes()
    num_edges_after = ets.num_edges
    extend_time = time.time() - start_time
    print(f"Extending ts done in {extend_time:.2f} seconds, reducing {num_edges_before} to {num_edges_after} edges.")

    start_time = time.time()
    _ = ets.Tajimas_D(mode="branch")
    runtime_after = time.time() - start_time

    results = vars(args)
    results["L"] = L
    results.update({
        "num_edges_before": num_edges_before,
        "num_edges_after": num_edges_after,
        "extend_time": extend_time,
        "runtime_before": runtime_before,
        "runtime_after": runtime_after,
    })
    with open(f"{args.output}.json", "w") as f:
        json.dump(results, f)
