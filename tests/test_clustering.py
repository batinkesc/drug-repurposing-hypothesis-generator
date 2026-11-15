"""Clustering properties the splits depend on.

The claims being pinned here are the ones the thesis will make in words:
sphere exclusion is deterministic and covers everything, single linkage
*guarantees* separation while sphere exclusion does not, and the coverage
threshold stops a short motif from chaining unrelated proteins together.
"""

from __future__ import annotations

import numpy as np
import pytest

from dti.clustering import (
    connected_components,
    identity_clusters,
    sphere_exclusion,
)
from dti.molecules import fingerprint_matrix, make_fingerprint_generator, tanimoto_matrix

# Three obvious chemical groups plus one outlier: short alcohols, benzene
# derivatives, and a sugar.
GROUPED_SMILES = [
    "CCO", "CCCO", "CCCCO", "CCCCCO",
    "c1ccccc1C", "c1ccccc1CC", "c1ccccc1CCC",
    "OCC1OC(O)C(O)C(O)C1O",
]


@pytest.fixture(scope="module")
def fps() -> np.ndarray:
    return fingerprint_matrix(GROUPED_SMILES, make_fingerprint_generator())


def test_every_compound_lands_in_a_cluster(fps):
    labels = sphere_exclusion(fps, 0.4)
    assert (labels >= 0).all()
    assert len(labels) == len(GROUPED_SMILES)


def test_sphere_exclusion_is_deterministic(fps):
    """Reproducibility is not optional: the cutoff is an ablation parameter, so
    the clustering is re-run and the two runs have to be comparable."""
    first = sphere_exclusion(fps, 0.4)
    second = sphere_exclusion(fps, 0.4)
    assert np.array_equal(first, second)


def test_a_high_cutoff_isolates_and_a_low_cutoff_merges(fps):
    strict = sphere_exclusion(fps, 0.95)
    loose = sphere_exclusion(fps, 0.05)
    assert len(np.unique(strict)) == len(GROUPED_SMILES)
    assert len(np.unique(loose)) < len(np.unique(strict))


def test_cluster_members_are_within_the_cutoff_of_their_centre(fps):
    """The property sphere exclusion actually provides -- and the one it does
    not. Members are close to their centre; two members of DIFFERENT clusters
    are not guaranteed to be far apart, which is exactly why the leakage audit
    measures cross-cluster similarity instead of assuming it."""
    cutoff = 0.4
    labels = sphere_exclusion(fps, cutoff)
    sim = tanimoto_matrix(fps, fps)
    for label in np.unique(labels):
        members = np.flatnonzero(labels == label)
        # The centre is the member that is within the cutoff of all others.
        reachable = [
            m for m in members if (sim[m, members] >= cutoff).all()
        ]
        assert reachable, f"cluster {label} has no centre within the cutoff"


def test_connected_components_merges_chains():
    """Single linkage links a to c through b, which is the behaviour R3 needs:
    a test target must not be similar to ANY training target, even indirectly
    via a third protein that would otherwise bridge the two sides."""
    labels = connected_components(5, np.array([[0, 1], [1, 2], [3, 4]]))
    assert labels[0] == labels[1] == labels[2]
    assert labels[3] == labels[4]
    assert labels[0] != labels[3]


def test_isolated_items_become_singletons():
    labels = connected_components(4, np.empty((0, 2), dtype=np.int64))
    assert len(np.unique(labels)) == 4


def test_component_labels_do_not_depend_on_pair_order():
    a = connected_components(6, np.array([[0, 1], [2, 3], [1, 4]]))
    b = connected_components(6, np.array([[1, 4], [2, 3], [1, 0]]))
    assert np.array_equal(a, b)


def test_identity_clusters_require_both_thresholds():
    """A short motif matching at 95% identity is not the same protein.

    Without the coverage test, targets 0 and 2 would be linked through a hit
    covering a tenth of the sequence, and R3 would hold out a target group that
    merged two unrelated families.
    """
    identity = np.array(
        [
            [1.00, 0.95, 0.95],
            [0.95, 1.00, 0.10],
            [0.95, 0.10, 1.00],
        ],
        dtype=np.float32,
    )
    coverage = np.array(
        [
            [1.00, 0.90, 0.10],   # 0-1 well covered, 0-2 barely
            [0.90, 1.00, 0.10],
            [0.10, 0.10, 1.00],
        ],
        dtype=np.float32,
    )
    labels = identity_clusters(identity, coverage, 0.40, 0.80)
    assert labels[0] == labels[1]
    assert labels[2] != labels[0]


def test_identity_clusters_separate_everything_below_the_cutoff():
    identity = np.full((4, 4), 0.2, dtype=np.float32)
    np.fill_diagonal(identity, 1.0)
    coverage = np.ones((4, 4), dtype=np.float32)
    labels = identity_clusters(identity, coverage, 0.40, 0.80)
    assert len(np.unique(labels)) == 4
