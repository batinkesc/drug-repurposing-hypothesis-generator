"""The identity metric, and the trap it was chosen to avoid.

The first leakage audit reported a maximum train-test identity of 1.000 for
every regime, R3 included -- where no target sits on both sides at all. The
cause was the metric: local alignment finds the best-scoring segment, so two
unrelated proteins sharing a short motif come back at high "identity" over a
window covering a few percent of the sequence.

These tests pin the distinction, because it is the difference between an audit
that measures protein similarity and one that measures nothing.
"""

from __future__ import annotations

import numpy as np
import pytest

from dti.sequences import identity_matrix

# Two unrelated 120-residue sequences that share one identical 12-residue
# motif. Built from repeating blocks so the only alignable region is the motif.
MOTIF = "WYWYWYWYWYWY"
LEFT = "AAAAAAAAAAAAAAAAAAAAAAAA" * 2 + MOTIF + "AAAAAAAAAAAAAAAAAAAAAAAA" * 2
RIGHT = "GGGGGGGGGGGGGGGGGGGGGGGG" * 2 + MOTIF + "GGGGGGGGGGGGGGGGGGGGGGGG" * 2


@pytest.fixture(scope="module")
def motif_pair():
    identity, coverage, over_shorter = identity_matrix(
        [LEFT, RIGHT], workers=1, progress=False
    )
    return identity[0, 1], coverage[0, 1], over_shorter[0, 1]


def test_local_alignment_identity_is_inflated_by_a_shared_motif(motif_pair):
    """The trap, demonstrated: raw identity calls these two the same protein."""
    identity, coverage, _ = motif_pair
    assert identity > 0.8, "expected the motif to dominate the local alignment"
    assert coverage < 0.3, "expected the alignment to cover little of the protein"


def test_identity_over_shorter_is_not_inflated(motif_pair):
    """The fix: a 12-residue motif cannot account for 40% of a 120-mer."""
    _, _, over_shorter = motif_pair
    assert over_shorter < 0.2
    assert over_shorter == pytest.approx(len(MOTIF) / len(LEFT), abs=0.05)


def test_identical_sequences_score_one_on_every_measure():
    identity, coverage, over_shorter = identity_matrix(
        [LEFT, LEFT], workers=1, progress=False
    )
    assert identity[0, 1] == pytest.approx(1.0)
    assert coverage[0, 1] == pytest.approx(1.0)
    assert over_shorter[0, 1] == pytest.approx(1.0)


def test_a_contained_subsequence_scores_high_over_the_shorter_sequence():
    """A short domain wholly inside a longer one IS the same protein region,
    and the metric has to say so -- this is the case the coverage threshold got
    wrong, since coverage over the shorter sequence is 1.0 but the pair would
    be missed if coverage were measured over the longer one."""
    short = LEFT[:60]
    _, coverage, over_shorter = identity_matrix([LEFT, short], workers=1, progress=False)
    assert over_shorter[0, 1] > 0.9
    assert coverage[0, 1] > 0.9


def test_matrices_are_symmetric_with_unit_diagonals():
    identity, coverage, over_shorter = identity_matrix(
        [LEFT, RIGHT, LEFT[:60]], workers=1, progress=False
    )
    for matrix in (identity, coverage, over_shorter):
        assert np.allclose(matrix, matrix.T)
        assert np.allclose(np.diag(matrix), 1.0)


def test_a_single_sequence_needs_no_alignment():
    identity, coverage, over_shorter = identity_matrix([LEFT], workers=1, progress=False)
    assert identity.shape == (1, 1)
    assert over_shorter[0, 0] == pytest.approx(1.0)
