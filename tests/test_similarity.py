"""The tiled GPU path must agree with the dense reference implementation.

This is the test that matters most in the whole suite. `molecules.tanimoto_matrix`
is obviously correct and obviously too slow; `similarity` is fast and full of
places to be subtly wrong -- fp16 rounding, tile boundaries, the division-free
threshold rearrangement, the self-exclusion mask. If the two ever disagree,
every split and every leakage number is suspect.
"""

from __future__ import annotations

import numpy as np
import pytest

from dti import similarity
from dti.molecules import (
    fingerprint_matrix,
    make_fingerprint_generator,
    tanimoto_matrix,
)

SMILES = [
    "CCO",
    "CCOC",
    "c1ccccc1",
    "c1ccccc1O",
    "c1ccccc1OC",
    "CC(=O)Oc1ccccc1C(=O)O",         # aspirin
    "CC(=O)Nc1ccc(O)cc1",            # paracetamol
    "CC(C)Cc1ccc(C(C)C(=O)O)cc1",    # ibuprofen
    "CN1C=NC2=C1C(=O)N(C)C(=O)N2C",  # caffeine
    "CCN(CC)CCOC(=O)c1ccc(N)cc1",    # procaine
    "OCC1OC(O)C(O)C(O)C1O",
    "NCCc1ccc(O)c(O)c1",             # dopamine
]


@pytest.fixture(scope="module")
def fps() -> np.ndarray:
    generator = make_fingerprint_generator(radius=2, n_bits=2048)
    return fingerprint_matrix(SMILES, generator)


@pytest.fixture(scope="module")
def reference(fps: np.ndarray) -> np.ndarray:
    return tanimoto_matrix(fps, fps)


def test_max_similarity_matches_dense_reference(fps, reference):
    best, idx = similarity.max_similarity(fps, fps)
    expected = reference.max(axis=1)
    assert np.allclose(best, expected, atol=1e-3)
    # Every molecule's closest match in the full set is itself.
    assert np.array_equal(idx, np.arange(len(SMILES)))


def test_self_exclusion_finds_the_second_best(fps, reference):
    best, idx = similarity.max_similarity(fps, fps, same_set=True)
    masked = reference.copy()
    np.fill_diagonal(masked, -1.0)
    assert np.allclose(best, masked.max(axis=1), atol=1e-3)
    assert (idx != np.arange(len(SMILES))).all()


def test_max_similarity_against_a_disjoint_set(fps, reference):
    query, ref = [0, 2, 4], [1, 3, 5, 7]
    best, idx = similarity.max_similarity(fps[query], fps[ref])
    expected_block = reference[np.ix_(query, ref)]
    assert np.allclose(best, expected_block.max(axis=1), atol=1e-3)
    assert np.array_equal(
        np.asarray(ref)[idx], np.asarray(ref)[expected_block.argmax(axis=1)]
    )


@pytest.mark.parametrize("cutoff", [0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 0.9])
def test_neighbour_counts_match_the_reference_at_every_cutoff(fps, reference, cutoff):
    counts = similarity.neighbour_counts(fps, cutoff)
    masked = reference.copy()
    np.fill_diagonal(masked, -1.0)
    expected = (masked >= cutoff).sum(axis=1)
    assert np.array_equal(counts, expected)


@pytest.mark.parametrize("cutoff", [0.2, 0.4, 0.6])
def test_neighbours_above_match_the_reference(fps, reference, cutoff):
    found = similarity.neighbours_above(fps, fps, cutoff)
    for i, hits in enumerate(found):
        expected = np.flatnonzero(reference[i] >= cutoff)
        # Self is included here (similarity 1.0 clears every cutoff); sphere
        # exclusion relies on that, since a cluster contains its own centre.
        assert np.array_equal(hits, expected), f"row {i} at cutoff {cutoff}"


def test_tiling_does_not_change_the_answer(fps, reference, monkeypatch):
    """Tile boundaries are where an off-by-one would hide."""
    monkeypatch.setattr(similarity, "QUERY_TILE", 1)
    monkeypatch.setattr(similarity, "REF_CHUNK", 2)
    best, _ = similarity.max_similarity(fps, fps, same_set=True)
    masked = reference.copy()
    np.fill_diagonal(masked, -1.0)
    assert np.allclose(best, masked.max(axis=1), atol=1e-3)
    counts = similarity.neighbour_counts(fps, 0.3)
    assert np.array_equal(counts, (masked >= 0.3).sum(axis=1))


def test_threshold_is_exact_on_the_cutoff():
    """A pair sitting exactly on the cutoff must resolve deterministically.

    Two fingerprints with 2 shared bits out of a 4-bit union have similarity
    exactly 0.5. The division-free comparison must call that `>= 0.5` true --
    if it were done in fp16 the value could land just below and the pair would
    silently leave the cluster.
    """
    a = np.zeros((1, 2048), dtype=np.uint8)
    b = np.zeros((1, 2048), dtype=np.uint8)
    a[0, [0, 1, 2]] = 1
    b[0, [1, 2, 3]] = 1          # intersection 2, union 4 -> exactly 0.5
    assert similarity.neighbour_counts(np.vstack([a, b]), 0.5).tolist() == [1, 1]
    assert similarity.neighbour_counts(np.vstack([a, b]), 0.51).tolist() == [0, 0]


def test_all_zero_fingerprint_has_no_neighbours(fps):
    """An unreadable SMILES becomes an all-zero row; it must not match
    everything through a zero/zero union."""
    padded = np.vstack([fps, np.zeros((1, fps.shape[1]), dtype=np.uint8)])
    best, _ = similarity.max_similarity(padded, padded, same_set=True)
    assert best[-1] == 0.0
    counts = similarity.neighbour_counts(padded, 0.0)
    assert counts[-1] == 0


def test_cpu_and_gpu_agree(fps):
    """If the GPU path drifted from the CPU path the pipeline would not be
    reproducible on a machine without a card."""
    import torch

    if not torch.cuda.is_available():
        pytest.skip("no GPU in this environment")
    cpu, _ = similarity.max_similarity(
        fps, fps, same_set=True, device=torch.device("cpu")
    )
    gpu, _ = similarity.max_similarity(
        fps, fps, same_set=True, device=torch.device("cuda")
    )
    assert np.allclose(cpu, gpu, atol=1e-3)
