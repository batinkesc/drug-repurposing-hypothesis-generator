"""Tests for molecule standardisation.

The load-bearing claim being tested is the one the R4b split rests on: salt,
charge and tautomer variants of the same compound must collapse to one
canonical SMILES. If they do not, the "chemically new ligand" definition leaks
-- the same molecule can sit on both sides of the split under two spellings.
"""

from __future__ import annotations

import numpy as np
import pytest

from dti.molecules import (
    Standardiser,
    fingerprint_matrix,
    make_fingerprint_generator,
    tanimoto_matrix,
)


@pytest.fixture(scope="module")
def std() -> Standardiser:
    return Standardiser()


def test_plain_molecule_round_trips(std):
    result = std.standardise("c1ccccc1C(=O)Nc1ccccc1")
    assert result.ok
    assert result.inchikey
    assert result.n_heavy_atoms == 15


def test_salt_is_stripped(std):
    """A sodium salt and its free acid must become the same molecule.

    Uses ibuprofen rather than acetate: the free acid has to clear
    `min_heavy_atoms`, which is applied after the fragment is chosen.
    """
    salt = std.standardise("CC(C)Cc1ccc(cc1)C(C)C(=O)[O-].[Na+]")
    acid = std.standardise("CC(C)Cc1ccc(cc1)C(C)C(=O)O")
    assert salt.ok and acid.ok
    assert salt.canonical_smiles == acid.canonical_smiles


def test_salt_whose_free_acid_is_below_the_size_floor_is_dropped(std):
    """The size floor applies to the stripped molecule, not the input.

    Sodium acetate becomes 4 heavy atoms, below the floor of 6, so it is
    dropped rather than kept as a fragment that ChemBERTa-2 would embed
    meaninglessly.
    """
    result = std.standardise("CC(=O)[O-].[Na+]")
    assert not result.ok
    assert result.reason == "too_small"
    assert result.n_heavy_atoms == 4


def test_counter_ion_is_stripped_keeping_the_drug(std):
    hydrochloride = std.standardise("CN(C)CCc1ccccc1.Cl")
    freebase = std.standardise("CN(C)CCc1ccccc1")
    assert hydrochloride.canonical_smiles == freebase.canonical_smiles


def test_charge_is_neutralised(std):
    charged = std.standardise("CCCC(=O)[O-]")
    neutral = std.standardise("CCCC(=O)O")
    assert charged.canonical_smiles == neutral.canonical_smiles


def test_tautomers_collapse_to_one_form(std):
    """Keto and enol spellings of the same compound must not survive as two."""
    keto = std.standardise("CC(=O)CC(=O)C")
    enol = std.standardise("CC(O)=CC(=O)C")
    assert keto.ok and enol.ok
    assert keto.canonical_smiles == enol.canonical_smiles


def test_different_input_spellings_agree(std):
    """Kekulised, aromatic and atom-reordered spellings of aspirin."""
    spellings = [
        "CC(=O)Oc1ccccc1C(=O)O",
        "CC(=O)OC1=CC=CC=C1C(=O)O",
        "O=C(O)c1ccccc1OC(C)=O",
    ]
    canonical = {std.standardise(s).canonical_smiles for s in spellings}
    assert len(canonical) == 1


def test_standardisation_is_idempotent(std):
    """Re-standardising a canonical SMILES must be a no-op."""
    first = std.standardise("CC(=O)[O-].[Na+]")
    second = std.standardise(first.canonical_smiles)
    assert second.canonical_smiles == first.canonical_smiles


def test_unparseable_smiles_is_reported_not_raised(std):
    result = std.standardise("this is not a molecule")
    assert not result.ok
    assert result.reason == "unparseable"


def test_empty_input_is_reported(std):
    assert std.standardise("").reason == "empty_smiles"
    assert std.standardise(None).reason == "empty_smiles"


def test_size_bounds_are_enforced(std):
    too_small = std.standardise("CCO")  # 3 heavy atoms, below min 6
    assert not too_small.ok
    assert too_small.reason == "too_small"

    peptide = std.standardise("C" * 150)  # long chain, above max 100
    assert not peptide.ok
    assert peptide.reason == "too_large"


def test_drop_reasons_are_distinguishable(std):
    """The drop table needs one reason per stage, not a generic failure."""
    reasons = {
        std.standardise("").reason,
        std.standardise("not a molecule").reason,
        std.standardise("CCO").reason,
    }
    assert reasons == {"empty_smiles", "unparseable", "too_small"}


# --- fingerprints and similarity ------------------------------------------


def test_fingerprint_matrix_shape_and_content():
    generator = make_fingerprint_generator(radius=2, n_bits=1024)
    matrix = fingerprint_matrix(["CCO", "c1ccccc1"], generator)
    assert matrix.shape == (2, 1024)
    assert matrix.dtype == np.uint8
    assert set(np.unique(matrix)) <= {0, 1}
    assert matrix.sum() > 0


def test_tanimoto_of_identical_molecules_is_one():
    generator = make_fingerprint_generator()
    matrix = fingerprint_matrix(["CC(=O)Oc1ccccc1C(=O)O"], generator)
    similarity = tanimoto_matrix(matrix, matrix)
    assert similarity[0, 0] == pytest.approx(1.0)


def test_tanimoto_of_dissimilar_molecules_is_low():
    generator = make_fingerprint_generator()
    a = fingerprint_matrix(["CCCCCCCC"], generator)
    b = fingerprint_matrix(["c1ccc2[nH]ccc2c1"], generator)
    assert tanimoto_matrix(a, b)[0, 0] < 0.2


def test_tanimoto_of_analogues_is_high():
    """The R4b filter exists because analogues score high: ethyl vs propyl."""
    generator = make_fingerprint_generator()
    a = fingerprint_matrix(["CCc1ccccc1C(=O)NC"], generator)
    b = fingerprint_matrix(["CCCc1ccccc1C(=O)NC"], generator)
    assert tanimoto_matrix(a, b)[0, 0] > 0.4


def test_tanimoto_matrix_is_rectangular_and_symmetric_where_it_should_be():
    generator = make_fingerprint_generator()
    train = fingerprint_matrix(["CCO", "c1ccccc1", "CC(=O)O"], generator)
    test = fingerprint_matrix(["CCCO", "c1ccncc1"], generator)

    similarity = tanimoto_matrix(test, train)
    assert similarity.shape == (2, 3)
    assert np.all(similarity >= 0) and np.all(similarity <= 1)

    # What the R4b filter actually asks: the nearest training neighbour.
    max_per_test = similarity.max(axis=1)
    assert max_per_test.shape == (2,)


def test_tanimoto_handles_an_all_zero_row():
    """An empty fingerprint must give 0, not a nan that silently passes a
    similarity cutoff."""
    a = np.zeros((1, 64), dtype=np.uint8)
    b = np.zeros((1, 64), dtype=np.uint8)
    b[0, :5] = 1
    assert tanimoto_matrix(a, b)[0, 0] == 0.0
    assert tanimoto_matrix(a, a)[0, 0] == 0.0
