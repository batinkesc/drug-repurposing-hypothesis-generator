"""Evaluation metrics.

Primary: RMSE and Pearson on the pChEMBL scale. Secondary: Spearman and the
concordance index, which measure ranking quality -- what actually matters in the
repurposing scenario.

The one non-obvious part is `rmse_variants`. Families differ in their pChEMBL
distribution, so a model transferred to a held-out family tends to make a
systematic offset error: raw RMSE explodes while the ranking survives. Reporting
raw RMSE alone would call that "transfer failed", which is a different finding
from "transfer shifted". The three variants separate the two.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import stats

# Above this many test records, the concordance index is estimated from a
# random sample of pairs instead of all of them: the exact computation is
# O(n^2) and 20k records already means 2e8 pairs.
CI_EXACT_MAX_N = 5000
CI_SAMPLE_PAIRS = 2_000_000


def _clean(y_true: np.ndarray, y_pred: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Drop non-finite entries, keeping the two arrays aligned."""
    y_true = np.asarray(y_true, dtype=float).ravel()
    y_pred = np.asarray(y_pred, dtype=float).ravel()
    if y_true.shape != y_pred.shape:
        raise ValueError(f"shape mismatch: {y_true.shape} vs {y_pred.shape}")
    keep = np.isfinite(y_true) & np.isfinite(y_pred)
    return y_true[keep], y_pred[keep]


def rmse(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    y_true, y_pred = _clean(y_true, y_pred)
    if y_true.size == 0:
        return float("nan")
    return float(np.sqrt(np.mean((y_true - y_pred) ** 2)))


def mae(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    y_true, y_pred = _clean(y_true, y_pred)
    if y_true.size == 0:
        return float("nan")
    return float(np.mean(np.abs(y_true - y_pred)))


def pearson(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    y_true, y_pred = _clean(y_true, y_pred)
    # A constant prediction has no correlation to report. This is not a corner
    # case here: B1 in the out-of-family regime predicts one value for
    # everything, so it would otherwise emit a divide-by-zero warning and nan.
    if y_true.size < 2 or np.std(y_pred) == 0 or np.std(y_true) == 0:
        return float("nan")
    # Tuple indexing rather than `.statistic`: works across scipy versions.
    return float(stats.pearsonr(y_true, y_pred)[0])


def spearman(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    y_true, y_pred = _clean(y_true, y_pred)
    if y_true.size < 2 or np.std(y_pred) == 0 or np.std(y_true) == 0:
        return float("nan")
    return float(stats.spearmanr(y_true, y_pred)[0])


def concordance_index(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    *,
    seed: int = 0,
    exact_max_n: int = CI_EXACT_MAX_N,
    sample_pairs: int = CI_SAMPLE_PAIRS,
) -> float:
    """Fraction of comparable pairs ranked correctly; tied predictions score 0.5.

    This is the CI used in the DTI regression literature. Pairs with equal true
    values are not comparable and are excluded. For large inputs the value is
    estimated from a random sample of pairs with a fixed seed, so repeated runs
    on the same data agree.
    """
    y_true, y_pred = _clean(y_true, y_pred)
    n = y_true.size
    if n < 2:
        return float("nan")

    if n <= exact_max_n:
        true_diff = y_true[:, None] - y_true[None, :]
        pred_diff = y_pred[:, None] - y_pred[None, :]
        # Upper triangle only: each unordered pair once.
        upper = np.triu(np.ones((n, n), dtype=bool), k=1)
        comparable = upper & (true_diff != 0)
        n_comparable = int(comparable.sum())
        if n_comparable == 0:
            return float("nan")
        concordant = (np.sign(true_diff) == np.sign(pred_diff)) & comparable
        ties = (pred_diff == 0) & comparable
        return float((concordant.sum() + 0.5 * ties.sum()) / n_comparable)

    rng = np.random.default_rng(seed)
    i = rng.integers(0, n, size=sample_pairs)
    j = rng.integers(0, n, size=sample_pairs)
    keep = i != j
    i, j = i[keep], j[keep]
    true_diff = y_true[i] - y_true[j]
    pred_diff = y_pred[i] - y_pred[j]
    comparable = true_diff != 0
    if not comparable.any():
        return float("nan")
    true_diff, pred_diff = true_diff[comparable], pred_diff[comparable]
    concordant = np.sign(true_diff) == np.sign(pred_diff)
    ties = pred_diff == 0
    return float((concordant.sum() + 0.5 * ties.sum()) / true_diff.size)


@dataclass
class RmseVariants:
    """RMSE before and after correcting for a systematic shift or scaling.

    - `raw`: the model as it is. The only variant that reflects usability.
    - `mean_shifted`: predictions shifted so their mean matches the test set's.
      Answers "how much of the error is a constant offset?"
    - `linear_calibrated`: the best affine rescaling of the predictions.
      Answers "how much survives if scale is given away too?"

    An oracle supplies the correction in both corrected variants -- they use
    test-set statistics and are therefore NOT usable as a model. They are
    diagnostics that separate ranking transfer from scale transfer, and must be
    labelled as such wherever they are reported.
    """

    raw: float
    mean_shifted: float
    linear_calibrated: float
    shift: float = field(default=float("nan"))
    slope: float = field(default=float("nan"))

    def as_dict(self) -> dict[str, float]:
        return {
            "rmse_raw": self.raw,
            "rmse_mean_shifted": self.mean_shifted,
            "rmse_linear_calibrated": self.linear_calibrated,
            "calibration_shift": self.shift,
            "calibration_slope": self.slope,
        }


def rmse_variants(y_true: np.ndarray, y_pred: np.ndarray) -> RmseVariants:
    """Compute the three RMSE variants plus the corrections they applied.

    The reported `shift` is itself a result: in the out-of-family regime it is
    the size of the systematic offset between the training families and the
    held-out one.
    """
    y_true, y_pred = _clean(y_true, y_pred)
    if y_true.size == 0:
        nan = float("nan")
        return RmseVariants(nan, nan, nan)

    raw = float(np.sqrt(np.mean((y_true - y_pred) ** 2)))

    shift = float(np.mean(y_true) - np.mean(y_pred))
    shifted = float(np.sqrt(np.mean((y_true - (y_pred + shift)) ** 2)))

    # Least-squares affine fit of predictions onto truth. A constant prediction
    # carries no slope information, so the calibrated value falls back to
    # predicting the test mean -- which is exactly what it is worth.
    if y_true.size >= 2 and np.std(y_pred) > 0:
        slope, intercept = np.polyfit(y_pred, y_true, 1)
        calibrated = float(
            np.sqrt(np.mean((y_true - (slope * y_pred + intercept)) ** 2))
        )
    else:
        slope = float("nan")
        calibrated = float(np.sqrt(np.mean((y_true - np.mean(y_true)) ** 2)))

    return RmseVariants(
        raw=raw,
        mean_shifted=shifted,
        linear_calibrated=calibrated,
        shift=shift,
        slope=float(slope),
    )


def evaluate(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    *,
    include_rmse_variants: bool = False,
    seed: int = 0,
) -> dict[str, float]:
    """All metrics for one prediction set.

    `include_rmse_variants` is switched on for the R4 regimes, where the
    shift/scale decomposition is needed to interpret the result at all.
    """
    y_true, y_pred = _clean(y_true, y_pred)
    result: dict[str, float] = {
        "n": float(y_true.size),
        "rmse": rmse(y_true, y_pred),
        "mae": mae(y_true, y_pred),
        "pearson": pearson(y_true, y_pred),
        "spearman": spearman(y_true, y_pred),
        "ci": concordance_index(y_true, y_pred, seed=seed),
    }
    if include_rmse_variants:
        result.update(rmse_variants(y_true, y_pred).as_dict())
    return result


def aggregate_seeds(runs: list[dict[str, float]]) -> dict[str, float]:
    """Mean and standard deviation across repeated runs.

    Every configuration is repeated with several seeds and reported as mean plus
    spread, never a single run -- so this is what goes into the tables.
    """
    if not runs:
        return {}
    keys = runs[0].keys()
    out: dict[str, float] = {}
    for key in keys:
        values = np.array([run[key] for run in runs], dtype=float)
        finite = values[np.isfinite(values)]
        if finite.size == 0:
            out[f"{key}_mean"] = float("nan")
            out[f"{key}_std"] = float("nan")
        else:
            out[f"{key}_mean"] = float(finite.mean())
            # Population std: we report the spread of the runs we made, not an
            # estimate of a wider population.
            out[f"{key}_std"] = float(finite.std())
        out[f"{key}_n_runs"] = float(values.size)
    return out
