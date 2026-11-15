"""The leakage audit.

No score from any regime is interpretable without these numbers. A split can
satisfy its own definition and still be easy: R2 guarantees that a test
compound comes from a cluster no training compound came from, but Butina
clusters do not bound cross-cluster similarity (see `clustering`), so a test
molecule may still be a close analogue of a training molecule. The audit
measures that directly instead of trusting the construction.

Three quantities per regime:

* **identical molecules** shared between training and test. Should be zero in
  R2 and R4b by construction; a non-zero value is a bug, not a finding.
* **maximum train-test Tanimoto**, with the distribution behind it. The
  maximum alone can be one unlucky pair, so the share of test compounds above
  the cutoff is reported next to it -- that is the number that says whether
  the regime is testing what it claims.
* **maximum train-test domain identity**, the same question for the protein
  tower. Relevant in every regime, not only R3: a held-out family in R4 could
  still contain a domain similar to a training one.

Which identity, exactly, is not a detail. The first version of this audit used
identity over the local alignment and reported a maximum of 1.000 for *every*
regime, R3 included -- where no target is on both sides at all. The cause was
the metric: a local alignment covering 3% of two unrelated proteins at 70%
identity is a motif match, not protein similarity. The audit therefore reports
`identity_over_shorter` (identical residues / length of the shorter domain),
which a short motif cannot inflate. See `sequences`.

A note on what "zero" means per regime, since getting this backwards would
label correct behaviour as a bug:

* **R2 and R4b** must share no identical molecule: both are built to prevent it.
* **R1** shares molecules and targets by design; it is the optimistic baseline.
* **R3** holds out *targets*, so the same compound appearing on both sides,
  measured against different proteins, is correct and expected.
* **R4a** holds out a *family*, so a promiscuous compound measured against both
  a kinase and a GPCR legitimately appears on both sides. That this happens to
  roughly half the out-of-family test data is precisely why R4b exists.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .similarity import FingerprintSet, max_similarity, pick_device


def _percentiles(values: np.ndarray) -> dict[str, float]:
    if len(values) == 0:
        return {"p50": float("nan"), "p95": float("nan"), "max": float("nan")}
    return {
        "p50": float(np.percentile(values, 50)),
        "p95": float(np.percentile(values, 95)),
        "max": float(values.max()),
    }


def compound_similarity_audit(
    test_compounds: np.ndarray,
    train_compounds: np.ndarray,
    fingerprints: np.ndarray | FingerprintSet,
    cutoff: float,
    *,
    device=None,
    precomputed: np.ndarray | None = None,
) -> dict:
    """Train-test similarity for the compound tower.

    `test_compounds` and `train_compounds` are compound indices (positions in
    the fingerprint matrix), deduplicated by the caller: the similarity is a
    property of the molecule, so computing it per row would repeat identical
    queries thousands of times.

    `precomputed` is a per-compound array of maximum similarity to the training
    set, supplied when the caller has already computed it -- R4b needs exactly
    this array to do its filtering, and R4a and R4b share one training set, so
    the N^2 pass is done once and read three times.
    """
    device = device or pick_device()
    shared = np.intersect1d(test_compounds, train_compounds)
    if len(test_compounds) == 0 or len(train_compounds) == 0:
        return {
            "n_identical_molecules": int(len(shared)),
            "share_test_compounds_above_cutoff": float("nan"),
            **_percentiles(np.empty(0, dtype=np.float32)),
        }

    if precomputed is not None:
        best = np.asarray(precomputed)[test_compounds]
    else:
        if not isinstance(fingerprints, np.ndarray):
            raise TypeError(
                "pass the dense fingerprint matrix: each split needs a "
                "different subset, so a device-resident set would have to be "
                "re-uploaded"
            )
        best, _ = max_similarity(
            fingerprints[test_compounds],
            fingerprints[train_compounds],
            device=device,
        )
    return {
        "n_identical_molecules": int(len(shared)),
        "share_test_compounds_above_cutoff": float((best >= cutoff).mean()),
        **_percentiles(best),
    }


def sequence_identity_audit(
    test_targets: np.ndarray,
    train_targets: np.ndarray,
    identity: np.ndarray,
    cutoff: float,
) -> dict:
    """Train-test identity for the protein tower, from the precomputed matrix.

    `identity` must be the identity-over-shorter matrix, not the raw local
    alignment identity; see the module docstring for why the distinction
    changes the answer.
    """
    if len(test_targets) == 0 or len(train_targets) == 0:
        return {
            "n_shared_targets": 0,
            "share_test_targets_above_cutoff": float("nan"),
            **_percentiles(np.empty(0, dtype=np.float32)),
        }
    block = identity[np.ix_(test_targets, train_targets)]
    best = block.max(axis=1)
    return {
        "n_shared_targets": int(len(np.intersect1d(test_targets, train_targets))),
        "share_test_targets_above_cutoff": float((best >= cutoff).mean()),
        **_percentiles(best),
    }


def audit_split(
    split,
    compound_index: np.ndarray,
    target_index: np.ndarray,
    fingerprints: np.ndarray,
    identity: np.ndarray,
    *,
    tanimoto_cutoff: float,
    identity_cutoff: float,
    device=None,
    compound_max_similarity: np.ndarray | None = None,
) -> dict:
    """Every audit number for one split, as a flat record."""
    test_compounds = np.unique(compound_index[split.test])
    train_compounds = np.unique(compound_index[split.train])
    test_targets = np.unique(target_index[split.test])
    train_targets = np.unique(target_index[split.train])

    compounds = compound_similarity_audit(
        test_compounds,
        train_compounds,
        fingerprints,
        tanimoto_cutoff,
        device=device,
        precomputed=compound_max_similarity,
    )
    sequences = sequence_identity_audit(
        test_targets, train_targets, identity, identity_cutoff
    )

    return {
        "regime": split.regime,
        "seed": split.seed,
        "n_train": len(split.train),
        "n_val": len(split.val),
        "n_test": len(split.test),
        "test_fraction": split.fractions["test"],
        "n_test_compounds": len(test_compounds),
        "n_test_targets": len(test_targets),
        "identical_molecules": compounds["n_identical_molecules"],
        "tanimoto_p50": compounds["p50"],
        "tanimoto_p95": compounds["p95"],
        "tanimoto_max": compounds["max"],
        "share_compounds_above_tanimoto_cutoff": compounds[
            "share_test_compounds_above_cutoff"
        ],
        "shared_targets": sequences["n_shared_targets"],
        "identity_p50": sequences["p50"],
        "identity_p95": sequences["p95"],
        "identity_max": sequences["max"],
        "share_targets_above_identity_cutoff": sequences[
            "share_test_targets_above_cutoff"
        ],
    }


def family_composition(split, families: np.ndarray) -> pd.DataFrame:
    """Rows per family on each side.

    Reported for every regime, not just R4: if a regime's test set happens to
    be dominated by one family, an aggregate score says more about that family
    than about the regime.
    """
    records = []
    for fold in ("train", "val", "test"):
        idx = getattr(split, fold)
        counts = pd.Series(families[idx]).value_counts()
        for family, n in counts.items():
            records.append({"fold": fold, "family": family, "rows": int(n)})
    return pd.DataFrame(records)
