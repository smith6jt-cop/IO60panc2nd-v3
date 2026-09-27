"""Tests for v3/confound_tests.py, all on synthetic in-memory data."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import confound_tests as ct


# --------------------------------------------------------------------------
# 1. Permutation test: null behaviour and planted-effect power
# --------------------------------------------------------------------------


def _null_donor_df(seed: int, n_per_group: int = 12) -> pd.DataFrame:
    """Metric independent of batch, batch balanced within each disease stratum."""
    rng = np.random.default_rng(seed)
    rows = []
    for disease in ("ND", "T1D"):
        for i in range(n_per_group):
            batch = "S0" if i % 2 == 0 else "S1"
            rows.append(
                {
                    "donor_id": f"{disease}_{i}",
                    "disease_status": disease,
                    "stain_group": batch,
                    "islet_density_per_mm2": rng.normal(10.0, 2.0),
                }
            )
    return pd.DataFrame(rows)


def test_permutation_pvalues_roughly_uniform_under_null():
    # Under the null (no real batch effect), p-values across independent
    # datasets should not be systematically tiny.
    pvals = []
    for seed in range(30):
        df = _null_donor_df(seed)
        result = ct.permutation_test_association(
            df, "islet_density_per_mm2", "stain_group", covariate_type="categorical",
            n_permutations=200, seed=seed,
        )
        pvals.append(result["p_value"])
    pvals = np.array(pvals)
    frac_significant = float(np.mean(pvals < 0.05))
    # Loose bound: with a true null, we expect roughly 5% false positives,
    # not e.g. 50%+.
    assert frac_significant < 0.30
    assert np.mean(pvals) > 0.25


def _planted_effect_donor_df(seed: int, n_per_group: int = 12, effect: float = 8.0) -> pd.DataFrame:
    """Metric strongly shifted by batch, within each disease stratum."""
    rng = np.random.default_rng(seed)
    rows = []
    for disease in ("ND", "T1D"):
        for i in range(n_per_group):
            batch = "S0" if i % 2 == 0 else "S1"
            base = rng.normal(10.0, 1.0)
            shift = effect if batch == "S1" else 0.0
            rows.append(
                {
                    "donor_id": f"{disease}_{i}",
                    "disease_status": disease,
                    "stain_group": batch,
                    "islet_density_per_mm2": base + shift,
                }
            )
    return pd.DataFrame(rows)


def test_permutation_pvalue_small_under_planted_effect():
    df = _planted_effect_donor_df(seed=0)
    result = ct.permutation_test_association(
        df, "islet_density_per_mm2", "stain_group", covariate_type="categorical",
        n_permutations=2000, seed=0,
    )
    assert result["p_value"] < 0.01


def test_permutation_continuous_covariate_spearman():
    rng = np.random.default_rng(0)
    n = 40
    covariate = rng.normal(size=n)
    metric = 3 * covariate + rng.normal(scale=0.2, size=n)
    df = pd.DataFrame(
        {
            "donor_id": [f"D{i}" for i in range(n)],
            "disease_status": ["ND"] * n,
            "background_anchor": covariate,
            "islet_density_per_mm2": metric,
        }
    )
    result = ct.permutation_test_association(
        df, "islet_density_per_mm2", "background_anchor", covariate_type="continuous",
        n_permutations=1000, seed=0,
    )
    assert result["p_value"] < 0.01


def test_permutation_too_few_donors_returns_note():
    df = pd.DataFrame({"donor_id": ["D1"], "disease_status": ["ND"], "stain_group": ["S0"], "islet_density_per_mm2": [1.0]})
    result = ct.permutation_test_association(df, "islet_density_per_mm2", "stain_group")
    assert np.isnan(result["p_value"])
    assert "note" in result


# --------------------------------------------------------------------------
# 2. Per-stratum sensitivity
# --------------------------------------------------------------------------


def test_per_stratum_sensitivity_consistent_direction():
    rows = []
    for batch in ("S0", "S1"):
        for i in range(5):
            rows.append({"donor_id": f"{batch}_ND_{i}", "disease_status": "ND", "stain_group": batch, "metric": 10.0})
            rows.append({"donor_id": f"{batch}_T1D_{i}", "disease_status": "T1D", "stain_group": batch, "metric": 5.0})
    df = pd.DataFrame(rows)
    result = ct.per_stratum_sensitivity(df, "metric", "disease_status", "stain_group", "ND", "T1D")
    full_row = result[result["batch"] == "full_cohort"].iloc[0]
    assert full_row["contrast"] == pytest.approx(5.0)
    for _, row in result[result["batch"] != "full_cohort"].iterrows():
        assert row["same_sign_as_full"]
        assert row["contrast"] == pytest.approx(5.0)


def test_per_stratum_sensitivity_flags_sign_flip_and_small_n():
    rows = [
        {"donor_id": "a1", "disease_status": "ND", "stain_group": "S0", "metric": 10.0},
        {"donor_id": "a2", "disease_status": "ND", "stain_group": "S0", "metric": 10.0},
        {"donor_id": "a3", "disease_status": "T1D", "stain_group": "S0", "metric": 5.0},
        {"donor_id": "a4", "disease_status": "T1D", "stain_group": "S0", "metric": 5.0},
        # S1 has only one donor per disease group (below min_donors_per_batch)
        # and the contrast direction is reversed vs the full cohort.
        {"donor_id": "b1", "disease_status": "ND", "stain_group": "S1", "metric": 1.0},
        {"donor_id": "b2", "disease_status": "T1D", "stain_group": "S1", "metric": 9.0},
    ]
    df = pd.DataFrame(rows)
    result = ct.per_stratum_sensitivity(df, "metric", "disease_status", "stain_group", "ND", "T1D", min_donors_per_batch=2)

    s1_row = result[result["batch"] == "S1"].iloc[0]
    assert bool(s1_row["same_sign_as_full"]) is False
    assert bool(s1_row["flagged_small_n"]) is True

    s0_row = result[result["batch"] == "S0"].iloc[0]
    assert bool(s0_row["same_sign_as_full"]) is True
    assert bool(s0_row["flagged_small_n"]) is False


# --------------------------------------------------------------------------
# 3. Mixed-effects model
# --------------------------------------------------------------------------


def _mixedlm_synthetic_df(seed: int = 0, n_donors: int = 20, objects_per_donor: int = 15, beta_age: float = 0.8):
    rng = np.random.default_rng(seed)
    rows = []
    for d in range(n_donors):
        donor_id = f"D{d}"
        age = rng.uniform(20, 60)
        disease = rng.choice(["ND", "T1D"])
        stain_group = rng.choice(["S0", "S1", "S2"])
        sex = rng.choice(["F", "M"])
        aab_count = rng.integers(0, 3)
        donor_intercept = rng.normal(0, 1.0)
        for _o in range(objects_per_donor):
            noise = rng.normal(0, 1.0)
            outcome = beta_age * age + donor_intercept + noise
            rows.append(
                {
                    "donor_id": donor_id,
                    "age": age,
                    "disease_status": disease,
                    "stain_group": stain_group,
                    "sex": sex,
                    "aab_count": aab_count,
                    "outcome": outcome,
                }
            )
    return pd.DataFrame(rows)


def test_mixed_model_recovers_known_fixed_effect():
    df = _mixedlm_synthetic_df(seed=1, beta_age=0.8)
    result = ct.mixed_effects_model(df, outcome_col="outcome", donor_col="donor_id")
    assert "error" not in result, result.get("error")
    summary = result["model_summary"]
    age_row = summary[summary["term"] == "age"].iloc[0]
    assert age_row["coef"] == pytest.approx(0.8, abs=0.25)
    assert age_row["p_value"] < 0.01


def test_mixed_model_too_small_design_gives_clear_message():
    df = pd.DataFrame(
        {
            "donor_id": ["D1", "D1", "D2"],
            "age": [30, 30, 40],
            "disease_status": ["ND", "ND", "T1D"],
            "stain_group": ["S0", "S0", "S1"],
            "sex": ["F", "F", "M"],
            "aab_count": [0, 0, 1],
            "outcome": [1.0, 1.1, 2.0],
        }
    )
    result = ct.mixed_effects_model(df, outcome_col="outcome", donor_col="donor_id")
    assert "error" in result
    assert "too small" in result["error"].lower() or "fewer than" in result["error"].lower()


def test_mixed_model_flags_batch_nested_in_disease():
    # Every stain_group maps to exactly one disease_status -> disease and
    # batch effects cannot be separated.
    rows = []
    rng = np.random.default_rng(0)
    for d in range(12):
        disease = "ND" if d < 6 else "T1D"
        batch = "S0" if d < 6 else "S1"  # batch fully determined by disease
        for _o in range(5):
            rows.append(
                {
                    "donor_id": f"D{d}",
                    "age": rng.uniform(20, 60),
                    "disease_status": disease,
                    "stain_group": batch,
                    "sex": rng.choice(["F", "M"]),
                    "aab_count": rng.integers(0, 3),
                    "outcome": rng.normal(),
                }
            )
    df = pd.DataFrame(rows)
    result = ct.mixed_effects_model(df, outcome_col="outcome", donor_col="donor_id")
    assert "error" in result
    assert "confound" in result["error"].lower() or "nested" in result["error"].lower()


# --------------------------------------------------------------------------
# Donor metrics computation (integration with islet_benchmark_score)
# --------------------------------------------------------------------------


def test_compute_donor_metrics_from_match_df():
    key_df = pd.DataFrame(
        {"code": ["t1", "t2"], "image": ["imgA", "imgB"], "x": [0, 0], "y": [0, 0], "w": [1000, 1000], "h": [1000, 1000], "tile_kind": ["random", "random"]}
    )
    donors_df = pd.DataFrame(
        {
            "donor_id": ["D1", "D2"],
            "image": ["imgA", "imgB"],
            "disease_status": ["ND", "T1D"],
            "stain_group": ["S0", "S1"],
            "acq_run": ["r1", "r2"],
            "age": [30, 40],
            "aab_count": [0, 2],
            "sex": ["F", "M"],
        }
    )
    match_df = pd.DataFrame(
        [
            {"image": "imgA", "outcome": "tp", "size_class": "Islet"},
            {"image": "imgA", "outcome": "fp", "size_class": "EndocrineSingle"},
            {"image": "imgB", "outcome": "fn", "size_class": "Islet"},
        ]
    )
    metrics = ct.compute_donor_metrics(match_df, key_df, donors_df, pixel_size_um=1.0)
    d1 = metrics[metrics["donor_id"] == "D1"].iloc[0]
    d2 = metrics[metrics["donor_id"] == "D2"].iloc[0]
    assert d1["n_predicted_islets"] == 1
    assert d2["n_predicted_islets"] == 0  # the imgB row is a false negative, not a prediction
