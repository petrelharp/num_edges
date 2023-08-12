import msprime, tskit
import numpy as np
import collections

ts = msprime.sim_ancestry(4, sequence_length=1e6, recombination_rate=1e-8, population_size=100, random_seed=123)
ts2 = msprime.sim_ancestry(4, sequence_length=1e6, recombination_rate=1e-8, population_size=100, random_seed=123).extend_edges()

class NodeMapper:

    def __init__(self, ts, ts2, *args, **kwargs):
        """
        Identify, for each node in ts, matching nodes in ts2 as follows:
            1. put down infinite-sites mutations at rate `mutation_rate` on ts,
            2. map these to ts2 with `tskit.Tree.map_mutations`,
            3. record "matches" between the nodes adjacent to the mutation in ts
                and any mutations in ts2
        The "matches" are defined as follows: if the mutation in ts is on an edge
        from `p` to `c`, and a mapped mutation in ts2 is on a series of unary edges
        from `p2` through more non-sample nodes `n` to `c`, then we record matches
        between `p` and `p2`, `c` and `c2`, and both `p` and `c` to all nodes in
        `n`.  (This reflects uncertainty about where in a series of edges separated
        by unary nodes a mutation should be located.) Furthermore, samples always match,
        and non-samples are not allowed to match to samples.

        This is recorded in the `NodeMapper.node_map`, a dictionary with non-sample
        nodes in `ts` as keys and a Counter of non-sample nodes in `ts2` as values,
        recording the number of matches.
        """
        assert ts.num_samples == ts2.num_samples, "Numbers of samples must match."
        self.ts = ts
        self.ts2 = ts2
        self.node_map = {}
        self.extra_matches = 0
        self.num_mutations = 0
        self.num_uniquely_mapped = 0
        self._map_nodes(*args, **kwargs)

    def _add_to_node_map(self, n, n2, extra=False):
        if not self.ts.node(n).is_sample() and not self.ts2.node(n2).is_sample():
            if n not in self.node_map:
                self.node_map[n] = collections.Counter()
            self.node_map[n].update((n2,))
            if extra:
                self.extra_matches += 1

    def _map_nodes(self, rate, random_seed=None):
        mts = msprime.sim_mutations(self.ts, rate=rate, discrete_genome=False, keep=False, model=msprime.BinaryMutationModel(), random_seed=random_seed)
        self.num_mutations += mts.num_mutations

        t = mts.first()
        t2 = self.ts2.first()
        for v in mts.variants():
            assert len(v.site.mutations) == 1
            assert v.site.ancestral_state == "0"
            mut = v.site.mutations[0]
            t.seek(v.position)
            t2.seek(v.position)
            ancestral_state, mutations = t2.map_mutations(v.genotypes, v.alleles)
            n = mut.node
            p = t.parent(n)
            self.num_uniquely_mapped += (len(mutations) == 1)
            for mut2 in mutations:
                if mut2.derived_state == "1":
                    n2 = mut2.node
                    self._add_to_node_map(n, n2)
                    self._add_to_node_map(p, t2.parent(n2))
                    while t2.num_children(n2) == 1 and not self.ts2.node(n2).is_sample():
                        self._add_to_node_map(p, n2, extra=True)
                        n2 = t2.left_child(n2)
                        self._add_to_node_map(n, n2, extra=True)

    def match_nodes(self):
        """
        Identify, for each node in ts, the best-matching node in ts2, by returning
        a (num nodes in ts) x 3 array with the j-th row recording for node `j`
        in ts`, the following things (in this order):
            1. best-matching node in ts2 (or -1 if none match, ties resolved arbitrarily)
            2. number of matches to that node
            3. number of other matches
        """
        # Identify best match:
        #  for each node in ts, find the node in ts2 that it matches most
        #  columns are: matching node (if any), number of supporting mutations
        match = np.full((self.ts.num_nodes, 3), [-1, 0, 0])
        for n in self.ts.samples():
            match[n] = [n, 0, 0]
        for n in self.node_map:
            top = self.node_map[n].most_common(1)[0]
            match[n] = (top[0], top[1], self.node_map[n].total() - top[1])
        return match

    def prop_uniquely_matching(self):
        """
        Returns the proportion of mutations that map uniquely to ts2 (except
        for uncertainty having to do with unary edges).
        """
        return self.num_uniquely_mapped / self.num_mutations

    def prop_matching(self):
        """
        Returns the proportion of matches that are between best-matched nodes.

        When mapping to a sequence of k unary edges, the mutation matches the
        input parent and child to k nodes each, of which at most one can be
        correct. So, the reported fraction for "proportion of matches" is
        the number of matches supporting the best matches divided by the total
        number of matches minus 2 * (number of 'extra' edges), where each
        mutation mapped to a sequence of k unary edges contributes k-1 to the
        count of 'extra' edges.
        """
        match = self.match_nodes()
        return np.sum(match[:,1]) / (np.sum(match[:,1:3]) - self.extra_matches)

match_args = {'rate': 2e-8, 'random_seed': 345}
nm = NodeMapper(ts, ts2, **match_args)
match = nm.match_nodes()
correct = (match[:,0] == np.arange(ts.num_nodes))

print(f"Proportion matching: {nm.prop_matching()}")
print(f"Proportion uniquely mapping: {nm.prop_uniquely_matching()}")
print(f"Number of extra matches: {nm.extra_matches}")
print("node match matches non_matches correct?")
for x in np.column_stack([np.arange(ts.num_nodes), match, correct]):
    if not ts.node(x[0]).is_sample():
        print("\t".join(map(str, x)))

# correctly matched nodes are black, incorrect are green,
# with the node they are incorrectly matched to listed in parentheses
styles = [
    f".node.n{n.id} > .sym, .node.n{n.id} > .lab" + "{" + f"fill: {'black' if correct[n.id] else 'green'}" + "}"
    for n in ts.nodes()
]
node_labels = { a : f"{a}{'' if correct[a] else ' ('+str(b)+')'}" for a, b in enumerate(match[:,0]) }
mts = msprime.sim_mutations(ts, **match_args, discrete_genome=False, keep=False, model=msprime.BinaryMutationModel())
_ = mts.draw_svg("orig.svg", size=(1200, 700), style=" ".join(styles), node_labels=node_labels)
