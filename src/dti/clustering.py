"""Grouping compounds and targets so that a split can keep like with like.

Two different algorithms, for two different reasons -- and the difference is a
result the thesis has to state, not an implementation detail.

**Targets: single linkage.** R3 asks the model about a protein it has never
seen, which means no test target may sit above the identity cutoff from any
training target. The only grouping that *guarantees* that is connected
components of the "above cutoff" graph: if two targets are linked, they are in
the same group, so they land on the same side. Chaining (a-b-c where a and c
are dissimilar) is the usual objection to single linkage, but with 565 targets
the groups stay inspectable, and a guarantee is worth more here than compact
clusters.

**Compounds: Butina sphere exclusion.** The same single-linkage guarantee is
not available for 190,339 drug-like compounds: at a 0.4 Tanimoto cutoff the
neighbourhood graph is dense enough that connected components collapse into
one giant component covering most of chemical space, and no split exists at
all. Sphere exclusion instead picks cluster centres greedily from the densest
compounds outwards, which keeps clusters compact but, crucially, does NOT
bound the similarity between compounds in *different* clusters.

So R2 and R4a guarantee only "a test compound comes from a cluster no training
compound came from". Whether that is strict enough is an empirical question,
which is why the leakage audit measures the maximum train-test Tanimoto for
every regime and why R4b exists as an explicit post-split filter. Reporting
the audit number is what makes the weaker guarantee honest.
"""

from __future__ import annotations

import numpy as np

from .similarity import (
    neighbour_counts,
    neighbours_above,
    pick_device,
    to_device,
)

# Candidates whose neighbourhoods are fetched in one GPU pass. Bounds memory
# while keeping the greedy order exact: similarities do not depend on what has
# already been assigned, only the "is it still free?" test does, so computing a
# block ahead of time gives the same clusters as computing one at a time.
CANDIDATE_BLOCK = 256


def sphere_exclusion(
    fp: np.ndarray,
    cutoff: float,
    *,
    device=None,
    progress: bool = False,
) -> np.ndarray:
    """Butina/Taylor sphere-exclusion clusters. Returns a label per row.

    The order is deterministic: candidates are taken by descending neighbour
    count, ties broken by index. Without the tie-break the clustering would
    depend on the sort's internal order and the dataset would not be
    reproducible.
    """
    device = device or pick_device()
    # Uploaded once and reused for the counting pass and every candidate
    # block: re-sending 780 MB per block would cost more than the arithmetic.
    reference = to_device(fp, device)
    counts = neighbour_counts(reference, cutoff, device=device)
    n = len(counts)

    # Descending count, ascending index on ties: lexsort reads keys last-first.
    order = np.lexsort((np.arange(n), -counts))

    labels = np.full(n, -1, dtype=np.int64)
    next_label = 0

    for start in range(0, n, CANDIDATE_BLOCK):
        block = order[start : start + CANDIDATE_BLOCK]
        free = block[labels[block] < 0]
        if len(free) == 0:
            continue
        neighbourhoods = neighbours_above(
            fp[free], reference, cutoff, device=device
        )
        for candidate, neighbours in zip(free, neighbourhoods):
            if labels[candidate] >= 0:
                # Claimed by an earlier candidate in this same block.
                continue
            labels[candidate] = next_label
            unclaimed = neighbours[labels[neighbours] < 0]
            labels[unclaimed] = next_label
            next_label += 1
        if progress:
            done = int((labels >= 0).sum())
            print(
                f"\r  clustering {done:,} / {n:,} ({100 * done / n:.1f}%)  "
                f"{next_label:,} clusters",
                end="",
                flush=True,
            )
    if progress:
        print()

    if (labels < 0).any():
        raise AssertionError("sphere exclusion left compounds unassigned")
    return labels


def connected_components(n: int, pairs: np.ndarray) -> np.ndarray:
    """Single-linkage groups from a list of linked index pairs.

    Union-find with path compression. `pairs` is an (m, 2) array; isolated
    items become singleton groups. Labels are renumbered in order of first
    appearance so they do not depend on the union order.
    """
    parent = np.arange(n, dtype=np.int64)

    def find(x: int) -> int:
        root = x
        while parent[root] != root:
            root = int(parent[root])
        while parent[x] != root:   # path compression
            parent[x], x = root, int(parent[x])
        return root

    for a, b in np.asarray(pairs, dtype=np.int64).reshape(-1, 2):
        ra, rb = find(int(a)), find(int(b))
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)

    roots = np.array([find(i) for i in range(n)], dtype=np.int64)
    _, labels = np.unique(roots, return_inverse=True)
    return labels.astype(np.int64)


def identity_clusters(
    identity: np.ndarray,
    coverage: np.ndarray,
    min_identity: float,
    min_coverage: float,
) -> np.ndarray:
    """Single-linkage target groups from an identity/coverage matrix.

    A pair counts as linked only if it clears BOTH thresholds: a short motif
    matching at 90% identity over 10% of the protein is not the same protein,
    and without the coverage test such hits would chain unrelated families
    together.

    The pipeline now passes `identity_over_shorter` here with `min_coverage=0`,
    because that metric already carries the protection the coverage threshold
    was added for -- a ten-residue motif cannot account for 40% of a 250-residue
    domain. The two-argument form is kept because it is what makes the coverage
    requirement testable, and because the earlier definition is an ablation.

    Using one definition for clustering and another for the audit was measured
    and rejected: with identity+coverage clustering and identity-over-shorter
    auditing, 10.2% of R3 test targets came back above the cutoff -- not a leak
    in the split so much as two definitions disagreeing about what "the same
    protein" means. One definition, used for both, makes the R3 guarantee hold
    by construction and the audit confirm it at 0%.
    """
    if identity.shape != coverage.shape or identity.ndim != 2:
        raise ValueError("identity and coverage must be the same square matrix shape")
    n = identity.shape[0]
    linked = (identity >= min_identity) & (coverage >= min_coverage)
    linked = np.triu(linked, k=1)
    pairs = np.argwhere(linked)
    return connected_components(n, pairs)
