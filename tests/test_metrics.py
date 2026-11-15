"""Tests for the evaluation metrics.

The cases that matter here are the degenerate ones, because they are not corner
cases in this project: B1 in the out-of-family regime predicts a single constant
for everything, and the R4 shift decomposition is the whole point of E4.
"""

from __future__ import annotations

import numpy as np
import pytest

from dti import metrics


def test_rmse_and_mae_on_known_values():
    y_true = np.array([1.0, 2.0, 3.0])
    y_pred = np.array([2.0, 2.0, 2.0])
    assert metrics.rmse(y_true, y_pred) == pytest.approx(np.sqrt(2 / 3))
    assert metrics.mae(y_true, y_pred) == pytest.approx(2 / 3)


def test_perfect_prediction():
    y = np.array([5.0, 6.0, 7.0, 8.0])
    assert metrics.rmse(y, y) == pytest.approx(0.0)
    assert metrics.pearson(y, y) == pytest.approx(1.0)
    assert metrics.concordance_index(y, y) == pytest.approx(1.0)


def test_non_finite_entries_are_dropped_pairwise():
    y_true = np.array([1.0, 2.0, np.nan, 4.0])
    y_pred = np.array([1.0, np.inf, 3.0, 4.0])
    # Only indices 0 and 3 survive, and both are exact.
    assert metrics.rmse(y_true, y_pred) == pytest.approx(0.0)


def test_constant_prediction_yields_nan_correlation_not_a_crash():
    """B1 in the out-of-family regime predicts one value for everything."""
    y_true = np.array([4.0, 5.0, 6.0, 7.0])
    y_pred = np.full(4, 5.5)
    assert np.isnan(metrics.pearson(y_true, y_pred))
    assert np.isnan(metrics.spearman(y_true, y_pred))
    # Every pair is a prediction tie, so CI is exactly 0.5 -- no ranking signal.
    assert metrics.concordance_index(y_true, y_pred) == pytest.approx(0.5)


def test_concordance_index_counts_ties_as_half():
    y_true = np.array([1.0, 2.0, 3.0])
    # Pairs: (1,2) correct, (1,3) correct, (2,3) tied in prediction.
    y_pred = np.array([1.0, 2.0, 2.0])
    assert metrics.concordance_index(y_true, y_pred) == pytest.approx(
        (2 + 0.5) / 3
    )


def test_concordance_index_of_reversed_ranking_is_zero():
    y_true = np.array([1.0, 2.0, 3.0, 4.0])
    assert metrics.concordance_index(y_true, -y_true) == pytest.approx(0.0)


def test_concordance_index_excludes_incomparable_pairs():
    """Pairs with equal true values carry no ranking information."""
    y_true = np.array([1.0, 1.0, 2.0])
    y_pred = np.array([9.0, 0.0, 5.0])
    # The (0,1) pair is tied in truth and must be excluded. Of the remaining
    # two, (0,2) is ranked wrongly and (1,2) correctly.
    assert metrics.concordance_index(y_true, y_pred) == pytest.approx(0.5)


def test_concordance_index_sampling_approximates_exact():
    rng = np.random.default_rng(0)
    n = 6000  # above CI_EXACT_MAX_N, so the sampling path is taken
    y_true = rng.normal(size=n)
    y_pred = y_true + rng.normal(scale=0.5, size=n)

    sampled = metrics.concordance_index(y_true, y_pred)
    exact = metrics.concordance_index(y_true, y_pred, exact_max_n=n)
    assert sampled == pytest.approx(exact, abs=0.01)


def test_concordance_index_sampling_is_deterministic():
    rng = np.random.default_rng(1)
    y_true = rng.normal(size=6000)
    y_pred = y_true + rng.normal(scale=0.5, size=6000)
    first = metrics.concordance_index(y_true, y_pred, sample_pairs=50_000)
    second = metrics.concordance_index(y_true, y_pred, sample_pairs=50_000)
    assert first == second


def test_pure_offset_is_removed_by_mean_shift():
    """The R4 case this decomposition exists for: ranking intact, scale shifted."""
    y_true = np.array([5.0, 6.0, 7.0, 8.0])
    y_pred = y_true - 2.0  # a perfect ranking, displaced by 2 log units

    variants = metrics.rmse_variants(y_true, y_pred)
    assert variants.raw == pytest.approx(2.0)
    assert variants.mean_shifted == pytest.approx(0.0)
    assert variants.linear_calibrated == pytest.approx(0.0)
    # The recovered shift is itself a reported result.
    assert variants.shift == pytest.approx(2.0)
    assert variants.slope == pytest.approx(1.0)


def test_scale_error_survives_mean_shift_but_not_calibration():
    y_true = np.array([2.0, 4.0, 6.0, 8.0])
    y_pred = y_true / 2.0  # right order, compressed range

    variants = metrics.rmse_variants(y_true, y_pred)
    assert variants.mean_shifted > 0.1
    assert variants.linear_calibrated == pytest.approx(0.0)
    assert variants.slope == pytest.approx(2.0)


def test_calibration_cannot_rescue_a_destroyed_ranking():
    """Calibration gives away scale, not ranking -- it must not inflate a
    model that has no ranking signal."""
    rng = np.random.default_rng(3)
    y_true = rng.normal(size=500)
    y_pred = rng.normal(size=500)  # independent of y_true

    variants = metrics.rmse_variants(y_true, y_pred)
    baseline = float(np.std(y_true))  # predicting the test mean
    assert variants.linear_calibrated == pytest.approx(baseline, rel=0.05)


def test_constant_prediction_calibrates_to_the_test_mean():
    y_true = np.array([4.0, 5.0, 6.0, 7.0])
    y_pred = np.full(4, 100.0)  # badly wrong, and carries no slope at all

    variants = metrics.rmse_variants(y_true, y_pred)
    assert variants.raw > 90
    assert variants.mean_shifted == pytest.approx(float(np.std(y_true)))
    assert variants.linear_calibrated == pytest.approx(float(np.std(y_true)))
    assert np.isnan(variants.slope)


def test_evaluate_includes_variants_only_when_asked():
    y_true = np.array([1.0, 2.0, 3.0])
    y_pred = np.array([1.5, 2.5, 3.5])

    plain = metrics.evaluate(y_true, y_pred)
    assert "rmse_raw" not in plain
    assert plain["n"] == 3

    with_variants = metrics.evaluate(y_true, y_pred, include_rmse_variants=True)
    assert with_variants["rmse_mean_shifted"] == pytest.approx(0.0)


def test_aggregate_seeds_reports_mean_and_spread():
    runs = [
        {"rmse": 1.0, "pearson": 0.5},
        {"rmse": 2.0, "pearson": 0.7},
        {"rmse": 3.0, "pearson": 0.6},
    ]
    out = metrics.aggregate_seeds(runs)
    assert out["rmse_mean"] == pytest.approx(2.0)
    assert out["rmse_std"] == pytest.approx(np.std([1.0, 2.0, 3.0]))
    assert out["rmse_n_runs"] == 3


def test_aggregate_seeds_tolerates_nan_runs():
    """A seed whose Pearson is undefined must not nan out the whole row."""
    runs = [{"pearson": float("nan")}, {"pearson": 0.4}, {"pearson": 0.6}]
    out = metrics.aggregate_seeds(runs)
    assert out["pearson_mean"] == pytest.approx(0.5)
    assert out["pearson_n_runs"] == 3
