# MIT License
#
# Copyright (c) 2024 Tskit Developers
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
"""
Tools for comparing node times between tree sequences with different node sets
"""
import copy
from collections import defaultdict
from dataclasses import dataclass
from itertools import product

import numpy as np
import scipy.sparse
import tskit


def node_spans(ts):
    """
    Returns the array of "node spans", i.e., the `j`th entry gives
    the total span over which node `j` is in the tree sequence
    (i.e., does not have 'missing data' there).

    """
    child_spans = np.bincount(
        ts.edges_child,
        weights=ts.edges_right - ts.edges_left,
        minlength=ts.num_nodes,
    )
    for t in ts.trees():
        span = t.span
        for r in t.roots:
            # do this check to exempt 'missing data'
            if t.num_children(r) > 0:
                child_spans[r] += span
    return child_spans


class CladeMap:
    """
    An iterator across trees that maintains a mapping from a clade (a `frozenset` of
    sample IDs) to a `set` of nodes. When there are unary nodes, there may be multiple
    nodes associated with each clade.
    """

    def __init__(self, ts):
        self._nil = frozenset()
        self._nodes = defaultdict(set)  # nodes[clade] = {node ids}
        self._clades = defaultdict(frozenset)  # clades[node] = {sample ids}
        self.tree_sequence = ts
        self.tree = ts.first(sample_lists=True)
        for node in self.tree.nodes():
            clade = frozenset(self.tree.samples(node))
            self._nodes[clade].add(node)
            self._clades[node] = clade
        self._prev = copy.deepcopy(self._clades)
        self._diff = ts.edge_diffs()
        next(self._diff)

    def _propagate(self, edge, downdate=False):
        """
        Traverse path from `edge.parent` to root, either adding or removing the
        state (clade) associated with `edge.child` from the state of each
        visited node. Return a set with the node ids encountered during
        traversal.
        """
        nodes = set()
        node = edge.parent
        clade = self._clades[edge.child]
        while node != tskit.NULL:
            last = self._clades[node]
            self._clades[node] = last - clade if downdate else last | clade
            if len(last):
                self._nodes[last].remove(node)
                if len(self._nodes[last]) == 0:
                    del self._nodes[last]
            self._nodes[self._clades[node]].add(node)
            nodes.add(node)
            node = self.tree.parent(node)
        return nodes

    def next(self):  # noqa: A003
        """
        Advance to the next tree, returning the difference between trees as a
        dictionary of the form `node : (last_clade, next_clade)`
        """
        nodes = set()  # nodes with potentially altered clades
        diff = {}  # diff[node] = (prev_clade, curr_clade)

        if self.tree.index + 1 == self.tree_sequence.num_trees:
            return None

        # Subtract clades subtended by outgoing edges
        edge_diff = next(self._diff)
        for eo in edge_diff.edges_out:
            nodes |= self._propagate(eo, downdate=True)

        # Prune nodes that are no longer in tree
        for node in self._nodes[self._nil]:
            diff[node] = (self._prev[node], self._nil)
            del self._clades[node]
        nodes -= self._nodes[self._nil]
        self._nodes[self._nil].clear()

        # Add clades subtended by incoming edges
        self.tree.next()
        for ei in edge_diff.edges_in:
            nodes |= self._propagate(ei, downdate=False)

        # Find difference in clades between adjacent trees
        for node in nodes:
            diff[node] = (self._prev[node], self._clades[node])
            if self._prev[node] == self._clades[node]:
                del diff[node]

        # Sync previous and current states
        for node, (_, curr) in diff.items():
            if curr == self._nil:
                del self._prev[node]
            else:
                self._prev[node] = curr

        return diff

    @property
    def interval(self):
        """
        Return interval spanned by tree
        """
        return self.tree.interval

    def clades(self):
        """
        Return set of clades in tree
        """
        return self._nodes.keys() - self._nil

    def __getitem__(self, clade):
        """
        Return set of nodes associated with a given clade.
        """
        return frozenset(self._nodes[clade]) if frozenset(clade) in self else self._nil

    def __contains__(self, clade):
        """
        Check if a clade is present in the tree
        """
        return clade in self._nodes


def shared_node_spans(ts, other):
    """
    Calculate the spans over which pairs of nodes in two tree sequences are
    ancestral to identical sets of samples.


    Returns a sparse matrix where rows correspond to nodes in `ts` and columns
    correspond to nodes in `other`, and whose value is the total amount of span
    over which the set of samples inheriting from the two nodes is identical.

    :return: A sparse matrix of class `scipy.sparse.csr_matrix`.
    """

    if ts.sequence_length != other.sequence_length:
        raise ValueError("Tree sequences must be of equal sequence length.")

    if ts.num_samples != other.num_samples:
        raise ValueError("Tree sequences must have the same numbers of samples.")

    nil = frozenset()

    # Initialize clade iterators
    query = CladeMap(ts)
    target = CladeMap(other)

    # Initialize buffer[clade] = (query_nodes, target_nodes, left_coord)
    modified = query.clades() | target.clades()
    buffer = {}

    # Build sparse matrix of matches in triplet format
    query_node = []
    target_node = []
    shared_span = []
    right = 0
    while True:
        left = right
        right = min(query.interval[1], target.interval[1])

        # Flush pairs of nodes that no longer have matching clades
        for clade in modified:  # flush:
            if clade in buffer:
                n_i, n_j, start = buffer.pop(clade)
                span = left - start
                for i, j in product(n_i, n_j):
                    query_node.append(i)
                    target_node.append(j)
                    shared_span.append(span)

        # Add new pairs of nodes with matching clades
        for clade in modified:
            assert clade not in buffer
            if clade in query and clade in target:
                n_i, n_j = query[clade], target[clade]
                buffer[clade] = (n_i, n_j, left)

        if right == ts.sequence_length:
            break

        # Find difference in clades with advance to next tree
        modified.clear()
        for clade_map in (query, target):
            if clade_map.interval[1] == right:
                clade_diff = clade_map.next()
                for prev, curr in clade_diff.values():
                    if prev != nil:
                        modified.add(prev)
                    if curr != nil:
                        modified.add(curr)

    # Flush final tree
    for clade in buffer:
        n_i, n_j, start = buffer[clade]
        span = right - start
        for i, j in product(n_i, n_j):
            query_node.append(i)
            target_node.append(j)
            shared_span.append(span)

    numer = scipy.sparse.coo_matrix(
        (shared_span, (query_node, target_node)),
        shape=(ts.num_nodes, other.num_nodes),
    ).tocsr()

    return numer


def match_node_ages(ts, other):
    """
        For each node in `ts`, return the age of a matched node from `other`.  Node
        matching is accomplished as described in :func:`.compare`.


        Returns a tuple of three vectors of length `ts.num_nodes`, in this order:
        the age of the best matching node in `other`;
        the proportion of the node span in `ts` that is covered by the best match;
        and the node id of the best match in `other`.


    :return: A tuple of arrays of length `ts.num_nodes` containing
        (time of matching node, proportion overlap, and node ID of match).
    """

    shared_spans = shared_node_spans(ts, other)
    matched_span = shared_spans.max(axis=1).todense().A1
    best_match = shared_spans.argmax(axis=1).A1
    # NB: if there are multiple nodes with the largest span in a row,
    # argmax returns the node with the smallest integer id
    matched_time = other.nodes_time[best_match]

    best_match[matched_span == 0] = tskit.NULL
    matched_time[matched_span == 0] = np.nan

    return matched_time, matched_span, best_match


@dataclass
class ARFResult:
    """
    The result of a call to tscompare.compare(ts, other),
    returning metrics associated with the ARG Robinson-Foulds
    measures of similarity and dissimilarity.
    This contains:

    `arf`:
        The ARG Robinson-Foulds relative dissimilarity:
        the proportion of the total span of `ts` that is *not* represented in `other`.
        This is: `dissimilarity / total_span[0]`

    `tpr`:
        The "true proportion represented":
        the proportion of the total span of `other` that is represented in `ts`.
        This is: `(total_span[0] - dissimilarity) / total_span[1]`

    `dissimilarity`:
        The total span of `ts` that is not represented in `other`.

    `inverse_dissimilarity`:
        The total span of `other` that is not represented in `ts`.

    `total_span`:
        The total of all node spans of the two tree sequences, in order (`ts`, `other`).

    `rmse`:
        The root-mean-squared error between the transformed times of the nodes in
        `ts` and the transformed times of their best-matching nodes in `other`, with
        the average taken weighting by span in `ts`.

    `transform`:
        The transformation function used to transform times for computing `rmse`.
    """

    arf: float
    tpr: float
    dissimilarity: float
    inverse_dissimilarity: float
    total_span: tuple
    rmse: float
    transform: callable

    def __str__(self):
        """
        Return a plain text summary of the ARF result.
        """
        out = "Tree sequence comparison:\n"
        out += f"    ARF: {100 * self.arf:.2f}%\n"
        out += f"    TPR: {100 * self.tpr:.2f}%\n"
        out += f"    dissimilarity: {self.dissimilarity}\n"
        out += f"    inverse_dissimilarity: {self.inverse_dissimilarity}\n"
        out += (
            f"    total span (ts, other): {self.total_span[0]}, {self.total_span[1]}\n"
        )
        out += f"    time RMSE: {self.rmse}\n"
        return out


def compare(ts, other, transform=None):
    """
    For two tree sequences `ts` and `other`,
    this method returns an object of type :class:`.ARFResult`.
    The values reported summarize the degree to which nodes in `ts`
    "match" corresponding nodes in `other`.

    To match nodes,
    for each node in `ts`, the best matching node(s) from `other`
    has the longest matching span using :func:`.shared_node_spans`.
    If there are multiple matches with the same longest shared span
    for a single node, the best match is the match that is closest in time.

    For each node in `other` we compute the best matched span
    as the maximum shared span amongst all nodes in `ts` which are its match.
    The similarity will then not exceed the total node span of `other`,
    bounding `tpr` to a proportion between 0 and 1.

    Then, :class:`.ARFResult` contains:

    - (`dissimilarity`)
        The total "matching span", which is the total span of
        all nodes in `ts` over which each node is ancestral to the same set of
        samples as its best match in `other`.

    - (`inverse_dissimiliarity`)
        The total "inverse matching span", which is the total
        span of all nodes in `other` over which each node is ancestral
        to the same set of samples as its best match in `ts`.

    - (`arf`)
        The fraction of the total span of `ts` over which each nodes'
        descendant sample set does not match its' best match's descendant
        sample set (i.e., the total *un*-matched span divided by the total
        span of `ts`).

    - (`tpr`)
        The proportion of the span in `other` that is correctly
        represented in `ts` (i.e., the total inverse matching span divided
        by the total span of `other`).

    - (`rmse`)
        The root mean squared difference
        between the transformed times of the nodes in `ts`
        and transformed times of their best matching nodes in `other`,
        with the average weighted by the nodes' spans in `ts`.

    - (`total_spans`)
        The total node spans of `ts` and `other`.

    The callable `transform` is used to transform times before computing
    root-mean-squared error (see :class:`.ARFResult`); the default
    is `log(1 + t)`.

    :param ts: The focal tree sequence.
    :param other: The tree sequence we compare to.
    :param transform: A callable that can take an array of times and
        return another array of numbers.
    :return: The three quantities above.
    :rtype: ARFResult
    """

    def f(t):
        return np.log(1 + t)

    if transform is None:
        transform = f

    shared_spans = shared_node_spans(ts, other)
    # Find all potential matches for a node based on max shared span length
    max_span = shared_spans.max(axis=1).toarray().flatten()
    col_ind = shared_spans.indices
    row_ind = np.repeat(
        np.arange(shared_spans.shape[0]), repeats=np.diff(shared_spans.indptr)
    )
    # mask to find all potential node matches
    match = shared_spans.data == max_span[row_ind]
    # scale with difference in node times
    # determine best matches with the best_match_matrix
    ts_times = ts.nodes_time[row_ind[match]]
    other_times = other.nodes_time[col_ind[match]]
    time_difference = np.absolute(
        np.asarray(transform(ts_times) - transform(other_times))
    )
    # If a node x in `ts` has no match then we set time_difference to zero
    # This node then does not effect the rmse
    for j in range(len(shared_spans.data[match])):
        if shared_spans.data[match][j] == 0:
            time_difference[j] = 0.0
    # If two nodes have the same time, then
    # time_difference is zero, which causes problems with argmin
    # Instead we store data as 1/(1+x) and find argmax
    best_match_matrix = scipy.sparse.coo_matrix(
        (
            1 / (1 + time_difference),
            (row_ind[match], col_ind[match]),
        ),
        shape=(ts.num_nodes, other.num_nodes),
    )
    # Between each pair of nodes, find the maximum shared span
    # n1_match is the matching N1 -> N2 (for arf, dissimilarity)
    # n2_match is finds the max match between nodes in N2 and their
    # best match in N1 based on max-span (for tpr, inverse_dissimilarity)
    best_n1_match = best_match_matrix.argmax(axis=1).A1
    n2_match_matrix = scipy.sparse.lil_matrix((ts.num_nodes, other.num_nodes))
    n2_match_matrix[np.arange(ts.num_nodes), best_n1_match] = best_match_matrix.tocsr()[
        np.arange(ts.num_nodes), best_n1_match
    ]
    n2_match_matrix = n2_match_matrix.tocsr()
    best_n2_match = n2_match_matrix.argmax(axis=0).A1
    n2_match_mask = best_n1_match[best_n2_match] == np.arange(other.num_nodes)
    best_match_n1_spans = shared_spans[np.arange(ts.num_nodes), best_n1_match].reshape(
        -1
    )
    best_match_n2_spans = shared_spans[
        best_n2_match, np.arange(other.num_nodes)
    ].reshape(-1)[0, n2_match_mask]
    total_match_n1_span = np.sum(best_match_n1_spans)
    total_match_n2_span = np.sum(best_match_n2_spans)
    ts_node_spans = node_spans(ts)
    total_span_ts = np.sum(ts_node_spans)
    total_span_other = np.sum(node_spans(other))
    # Compute the root-mean-square difference in transformed time
    # with the average weighted by span in ts
    time_matrix = scipy.sparse.csr_matrix(
        (time_difference, (row_ind[match], col_ind[match])),
        shape=(ts.num_nodes, other.num_nodes),
    )
    time_discrepancies = np.asarray(
        time_matrix[np.arange(len(best_n1_match)), best_n1_match].reshape(-1)
    )
    product = np.multiply((time_discrepancies**2), ts_node_spans)
    rmse = np.sqrt(np.sum(product) / total_span_ts)
    return ARFResult(
        arf=1.0 - total_match_n1_span / total_span_ts,
        tpr=total_match_n2_span / total_span_other,
        dissimilarity=total_span_ts - total_match_n1_span,
        inverse_dissimilarity=total_span_other - total_match_n2_span,
        total_span=(total_span_ts, total_span_other),
        rmse=rmse,
        transform=transform,
    )
