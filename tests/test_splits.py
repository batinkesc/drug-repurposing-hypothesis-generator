"""Split invariants.

Each test below corresponds to a sentence the thesis will assert about a
regime. The grouped regimes are checked for the property they claim -- no
group straddles train and test, *and* no group straddles train and validation,
because early stopping on a leaky validation set selects the model that best
exploits the leak.
"""

from __future__ import annotations

import numpy as np
import pytest

from dti.splits import (
    Split,
    assign_groups,
    filter_r4b,
    split_grouped,
    split_r1,
    split_r4a,
)


def test_overlapping_folds_are_rejected_on_construction():
    with pytest.raises(AssertionError, match="both train and test"):
        Split(
            regime="broken",
            seed=0,
            train=np.array([0, 1, 2]),
            val=np.array([3]),
            test=np.array([2, 4]),
        )


def test_r1_covers_every_row_exactly_once():
    split = split_r1(1000, seed=0, test_fraction=0.15, val_fraction=0.15)
    everything = np.concatenate([split.train, split.val, split.test])
    assert np.array_equal(np.sort(everything), np.arange(1000))


def test_r1_fractions_are_hit_closely():
    split = split_r1(10_000, seed=3, test_fraction=0.15, val_fraction=0.15)
    assert split.fractions["test"] == pytest.approx(0.15, abs=0.005)
    assert split.fractions["val"] == pytest.approx(0.15, abs=0.005)


def test_r1_is_reproducible_and_seed_dependent():
    a = split_r1(5000, seed=1, test_fraction=0.15, val_fraction=0.15)
    b = split_r1(5000, seed=1, test_fraction=0.15, val_fraction=0.15)
    c = split_r1(5000, seed=2, test_fraction=0.15, val_fraction=0.15)
    assert np.array_equal(a.test, b.test)
    assert not np.array_equal(a.test, c.test)


def _grouped_fixture(n_groups: int = 60, per_group: int = 20):
    groups = np.repeat(np.arange(n_groups), per_group)
    return groups


@pytest.mark.parametrize("seed", [0, 1, 2, 3, 4])
def test_no_group_straddles_any_two_folds(seed):
    groups = _grouped_fixture()
    split = split_grouped("R2", groups, seed, 0.15, 0.15, unit="compound_cluster")
    train_groups = set(groups[split.train].tolist())
    val_groups = set(groups[split.val].tolist())
    test_groups = set(groups[split.test].tolist())
    assert not (train_groups & test_groups)
    assert not (val_groups & test_groups)
    assert not (train_groups & val_groups)


def test_grouped_split_covers_every_row():
    groups = _grouped_fixture()
    split = split_grouped("R3", groups, 0, 0.15, 0.15, unit="target_group")
    everything = np.concatenate([split.train, split.val, split.test])
    assert np.array_equal(np.sort(everything), np.arange(len(groups)))


def test_grouped_fractions_are_approximate_not_exact():
    """Documented behaviour, not a defect: whole groups move, so the achieved
    fraction overshoots. The alternative -- stopping before a group that would
    overshoot -- would systematically prefer small groups for the test set."""
    groups = np.repeat(np.arange(10), 100)   # coarse groups, 10% each
    split = split_grouped("R2", groups, 0, 0.15, 0.15, unit="compound_cluster")
    assert split.fractions["test"] >= 0.15
    assert split.fractions["test"] <= 0.30


def test_one_dominant_group_does_not_break_the_split():
    """A single cluster holding most of the rows is the realistic failure mode
    for Butina on a large compound set."""
    groups = np.concatenate([np.zeros(900, dtype=int), np.arange(1, 101)])
    split = split_grouped("R2", groups, 0, 0.15, 0.15, unit="compound_cluster")
    everything = np.concatenate([split.train, split.val, split.test])
    assert np.array_equal(np.sort(everything), np.arange(len(groups)))
    assert len(split.test) > 0


def test_assign_groups_returns_a_fold_for_every_row():
    groups = _grouped_fixture()
    folds = assign_groups(groups, 0.15, 0.15, np.random.default_rng(0))
    assert set(np.unique(folds)) <= {"train", "val", "test"}
    assert not (folds == None).any()  # noqa: E711 -- object array


def test_r4a_holds_out_exactly_the_named_family():
    families = np.array(["kinase"] * 100 + ["gpcr"] * 80 + ["protease"] * 60)
    split = split_r4a(families, "gpcr", seed=0, val_fraction=0.15)
    assert set(families[split.test].tolist()) == {"gpcr"}
    assert "gpcr" not in set(families[split.train].tolist())
    assert "gpcr" not in set(families[split.val].tolist())


def test_r4a_validation_comes_from_training_families_only():
    """The rule that keeps R4 honest: a held-out-family validation set would
    leak that family into model selection."""
    families = np.array(["kinase"] * 100 + ["gpcr"] * 80)
    split = split_r4a(families, "gpcr", seed=1, val_fraction=0.20)
    assert len(split.val) > 0
    assert set(families[split.val].tolist()) == {"kinase"}


def test_r4a_rejects_an_absent_family():
    families = np.array(["kinase"] * 10)
    with pytest.raises(ValueError, match="no rows"):
        split_r4a(families, "gpcr", seed=0, val_fraction=0.15)


def _r4a_for_filter():
    families = np.array(["kinase"] * 60 + ["gpcr"] * 40)
    split = split_r4a(families, "gpcr", seed=0, val_fraction=0.10)
    # One compound per row, so compound index == row index.
    compound_index = np.arange(len(families))
    return split, compound_index


def test_r4b_removes_test_rows_above_the_cutoff():
    split, compound_index = _r4a_for_filter()
    similarity = np.zeros(100, dtype=np.float32)
    similarity[split.test[:10]] = 0.9        # too close to training
    filtered = filter_r4b(split, compound_index, similarity, 0.4, min_test_records=1)
    assert len(filtered.test) == len(split.test) - 10
    assert filtered.meta["n_removed"] == 10
    assert filtered.meta["max_train_test_similarity"] < 0.4


def test_r4b_keeps_training_and_validation_untouched():
    """Only the test side is filtered: dropping training data would change the
    model as well as the test, and the two effects could not be separated."""
    split, compound_index = _r4a_for_filter()
    similarity = np.full(100, 0.9, dtype=np.float32)
    filtered = filter_r4b(split, compound_index, similarity, 0.4, min_test_records=1)
    assert np.array_equal(filtered.train, split.train)
    assert np.array_equal(filtered.val, split.val)
    assert len(filtered.test) == 0


def test_r4b_marks_itself_unscorable_below_the_minimum():
    """The threshold is fixed in advance so that the decision not to score a
    family cannot be taken after seeing the model's result on it."""
    split, compound_index = _r4a_for_filter()
    similarity = np.full(100, 0.9, dtype=np.float32)
    similarity[split.test[:5]] = 0.0
    filtered = filter_r4b(split, compound_index, similarity, 0.4, min_test_records=200)
    assert len(filtered.test) == 5
    assert filtered.meta["scorable"] is False

    generous = filter_r4b(split, compound_index, similarity, 0.4, min_test_records=5)
    assert generous.meta["scorable"] is True


def test_r4b_renames_the_regime():
    split, compound_index = _r4a_for_filter()
    filtered = filter_r4b(
        split, compound_index, np.zeros(100, dtype=np.float32), 0.4, 1
    )
    assert filtered.regime == "R4b_gpcr"
    assert split.regime == "R4a_gpcr"
