"""The four split regimes.

Each regime answers a different generalisation question, and the only thing
that changes between them is what is held out:

* **R1** random rows -- the optimistic baseline every DTI paper reports.
* **R2** unseen compound -- test compounds come from Butina clusters no
  training compound came from.
* **R3** unseen target -- test targets come from identity groups no training
  target came from.
* **R4a** out-of-family -- one whole protein family is held out.
* **R4b** out-of-family *and* chemically novel -- R4a plus a post-split filter
  removing test compounds too similar to anything in training.

Two rules are applied uniformly, because breaking either would quietly make a
regime easier than it claims to be:

1. **Validation obeys the same constraint as test.** In R2 the validation
   compounds are from their own clusters, in R3 their own target groups. Early
   stopping on a validation set that violates the regime would select the model
   that best exploits the leak. In R4 validation comes from the training
   families only -- a held-out-family validation set would leak that family
   into model selection, which is the one thing R4 must not do.
2. **Groups, not rows, are assigned.** Fractions are therefore approximate;
   the achieved fractions are recorded in `meta` and reported, rather than
   assumed to equal the configured ones.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass(frozen=True)
class Split:
    """One train/val/test partition of the dataset's row indices."""

    regime: str
    seed: int
    train: np.ndarray
    val: np.ndarray
    test: np.ndarray
    meta: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        # Checked on construction rather than in a test: a row in two folds is
        # the one bug that would invalidate every number downstream while
        # leaving the pipeline looking healthy.
        for a, b in (("train", "test"), ("val", "test"), ("train", "val")):
            shared = np.intersect1d(
                getattr(self, a), getattr(self, b), assume_unique=False
            )
            if len(shared):
                raise AssertionError(
                    f"{self.regime}: {len(shared)} rows in both {a} and {b}"
                )
        self.meta["n_rows"] = len(self.train) + len(self.val) + len(self.test)

    @property
    def fractions(self) -> dict[str, float]:
        total = len(self.train) + len(self.val) + len(self.test)
        if total == 0:
            return {"train": 0.0, "val": 0.0, "test": 0.0}
        return {
            "train": len(self.train) / total,
            "val": len(self.val) / total,
            "test": len(self.test) / total,
        }


def assign_groups(
    groups: np.ndarray,
    test_fraction: float,
    val_fraction: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """Assign whole groups to folds, targeting row fractions.

    Returns an array of "train" / "val" / "test" per row. Groups are shuffled
    and then poured into test until its row quota is met, then into validation,
    with the remainder going to training. The group that crosses a quota is not
    split, so the achieved fraction overshoots slightly; this is preferred to
    the alternative of favouring small groups, which would systematically put
    the sparsest targets or the rarest chemotypes in the test set.
    """
    unique, counts = np.unique(groups, return_counts=True)
    order = rng.permutation(len(unique))
    total = int(counts.sum())
    test_quota = test_fraction * total
    val_quota = val_fraction * total

    fold_of_group: dict = {}
    test_rows = val_rows = 0
    for idx in order:
        group, size = unique[idx], int(counts[idx])
        if test_rows < test_quota:
            fold_of_group[group] = "test"
            test_rows += size
        elif val_rows < val_quota:
            fold_of_group[group] = "val"
            val_rows += size
        else:
            fold_of_group[group] = "train"

    folds = np.empty(len(groups), dtype=object)
    for group, fold in fold_of_group.items():
        folds[groups == group] = fold
    return folds


def _indices(folds: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    return (
        np.flatnonzero(folds == "train"),
        np.flatnonzero(folds == "val"),
        np.flatnonzero(folds == "test"),
    )


def split_r1(n_rows: int, seed: int, test_fraction: float, val_fraction: float) -> Split:
    """Random rows. The same compound and the same target appear on both sides."""
    rng = np.random.default_rng(seed)
    order = rng.permutation(n_rows)
    n_test = int(round(test_fraction * n_rows))
    n_val = int(round(val_fraction * n_rows))
    test, val, train = order[:n_test], order[n_test : n_test + n_val], order[n_test + n_val :]
    return Split(
        regime="R1",
        seed=seed,
        train=np.sort(train),
        val=np.sort(val),
        test=np.sort(test),
        meta={"unit": "row"},
    )


def split_grouped(
    regime: str,
    groups: np.ndarray,
    seed: int,
    test_fraction: float,
    val_fraction: float,
    *,
    unit: str,
) -> Split:
    """R2 and R3: identical mechanics, different grouping."""
    rng = np.random.default_rng(seed)
    folds = assign_groups(groups, test_fraction, val_fraction, rng)
    train, val, test = _indices(folds)
    n_groups = len(np.unique(groups))
    return Split(
        regime=regime,
        seed=seed,
        train=train,
        val=val,
        test=test,
        meta={
            "unit": unit,
            "n_groups": n_groups,
            "n_groups_test": len(np.unique(groups[test])) if len(test) else 0,
            "largest_group_rows": int(np.bincount(
                np.unique(groups, return_inverse=True)[1]
            ).max()),
        },
    )


def split_r4a(
    families: np.ndarray,
    held_out: str,
    seed: int,
    val_fraction: float,
) -> Split:
    """Leave one family out. Validation is drawn from the training families."""
    rng = np.random.default_rng(seed)
    test = np.flatnonzero(families == held_out)
    remaining = np.flatnonzero(families != held_out)
    if len(test) == 0:
        raise ValueError(f"no rows for held-out family {held_out!r}")

    shuffled = rng.permutation(remaining)
    n_val = int(round(val_fraction * len(remaining)))
    val, train = shuffled[:n_val], shuffled[n_val:]
    return Split(
        regime=f"R4a_{held_out}",
        seed=seed,
        train=np.sort(train),
        val=np.sort(val),
        test=np.sort(test),
        meta={"unit": "family", "held_out_family": held_out},
    )


def filter_r4b(
    split: Split,
    compound_index: np.ndarray,
    max_similarity_to_train: np.ndarray,
    cutoff: float,
    min_test_records: int,
) -> Split:
    """R4b: drop test rows whose compound is too close to the training set.

    `max_similarity_to_train` is indexed by compound, not by row, because the
    similarity is a property of the molecule -- computing it per row would
    repeat the same query for every target the compound was measured against.

    If too few records survive, the split is still returned but marked
    `scorable: False`. The threshold is fixed in the config in advance so that
    the decision not to score a family cannot be made after seeing how the
    model did on it.
    """
    test_similarity = max_similarity_to_train[compound_index[split.test]]
    keep = test_similarity < cutoff
    kept_test = split.test[keep]
    removed = int((~keep).sum())

    surviving = test_similarity[keep]
    return Split(
        regime=split.regime.replace("R4a", "R4b"),
        seed=split.seed,
        train=split.train,
        val=split.val,
        test=kept_test,
        meta={
            **split.meta,
            "r4b_cutoff": cutoff,
            "n_removed": removed,
            "n_test_before": len(split.test),
            "max_train_test_similarity": (
                float(surviving.max()) if len(surviving) else 0.0
            ),
            "scorable": len(kept_test) >= min_test_records,
            "min_test_records": min_test_records,
        },
    )
