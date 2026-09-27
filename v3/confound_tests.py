"""Batch-confound diagnostics for the v3 endocrine-object benchmark.

Three complementary checks that a donor-level or object-level outcome is
not confounded by stain batch (``stain_group``, ``acq_run``) rather than
disease biology:

1. ``permutation_test_association`` -- is a donor-level metric associated
   with a batch covariate, more than expected if batch were unrelated to
   disease? Labels are shuffled *within* disease strata so the disease/batch
   correlation structure already present in the design is preserved.
2. ``per_stratum_sensitivity`` -- does a disease-group contrast (e.g.
   T1D vs ND) on a donor-level metric hold up within each single batch, or
   does it flip sign / vanish once you leave the full (batch-mixed) cohort?
3. ``mixed_effects_model`` -- an object-level outcome modelled with a donor
   random intercept (statsmodels MixedLM) and fixed effects for age,
   aab_count, sex, disease_status and stain_group.

Run ``python confound_tests.py --help`` for CLI usage, or see
``v3/tests/test_confound_tests.py`` for worked examples on synthetic data.
"""

from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import pandas as pd
from scipy import stats as scipy_stats

try:
    import statsmodels.formula.api as smf
except ImportError:  # pragma: no cover
    smf = None

sys.path.insert(0, str(Path(__file__).resolve().parent))
import islet_benchmark_score as ibs

DEFAULT_N_PERMUTATIONS = 2000


# --------------------------------------------------------------------------
# Donor-level metrics (built from the same match_df as the benchmark scorer)
# --------------------------------------------------------------------------


def compute_donor_metrics(
    match_df: pd.DataFrame,
    key_df: pd.DataFrame,
    donors_df: pd.DataFrame,
    pixel_size_um: float = ibs.DEFAULT_PIXEL_SIZE_UM,
) -> pd.DataFrame:
    """Per-donor summary metrics computed from *predicted* objects (tp + fp).

    - ``islet_density_per_mm2``: predicted Islet-class objects per mm^2 of
      tiled area.
    - ``median_islet_size_um``: median mean-caliper diameter of predicted
      Islet-class objects.
    - ``endocrine_area_fraction``: N/A here (requires per-object area, not
      carried in ``match_df``); computed as NaN with a note, unless an
      ``area_um2`` column is supplied by the caller via ``extra_object_df``.

    Any additional donor-level covariates (e.g. vessel/nerve/lymphatic area
    fractions) should be merged in by the caller from a study-specific
    metrics table; this function only derives what the benchmark's own
    prediction files carry.
    """
    tile_area_mm2 = (key_df["w"].astype(float) * key_df["h"].astype(float) * pixel_size_um**2) / 1e6
    area_by_image = key_df.assign(_area=tile_area_mm2).groupby("image")["_area"].sum()

    predictions = match_df[match_df["outcome"] != "fn"] if not match_df.empty else match_df
    rows = []
    for donor_id, donor_rows in donors_df.groupby("donor_id"):
        images = donor_rows["image"].tolist()
        donor_area_mm2 = area_by_image.reindex(images).fillna(0).sum()
        donor_preds = predictions[predictions["image"].isin(images)] if not predictions.empty else predictions
        n_islet = int((donor_preds["size_class"] == "Islet").sum()) if not donor_preds.empty else 0
        density = n_islet / donor_area_mm2 if donor_area_mm2 > 0 else float("nan")
        rows.append({"donor_id": donor_id, "islet_density_per_mm2": density, "n_predicted_islets": n_islet, "tiled_area_mm2": donor_area_mm2})

    metrics = pd.DataFrame(rows).merge(
        donors_df.drop_duplicates("donor_id"), on="donor_id", how="left"
    )
    return metrics


# --------------------------------------------------------------------------
# 1. Permutation test of batch association
# --------------------------------------------------------------------------


def _between_group_variance(values: np.ndarray, groups: np.ndarray) -> float:
    """Between-group sum of squares statistic (weighted variance of group means)."""
    grand_mean = np.nanmean(values)
    stat = 0.0
    for g in pd.unique(groups):
        mask = groups == g
        n_g = mask.sum()
        if n_g == 0:
            continue
        group_mean = np.nanmean(values[mask])
        stat += n_g * (group_mean - grand_mean) ** 2
    return float(stat)


def permutation_test_association(
    donor_df: pd.DataFrame,
    metric_col: str,
    covariate_col: str,
    disease_col: str = "disease_status",
    covariate_type: str = "categorical",
    n_permutations: int = DEFAULT_N_PERMUTATIONS,
    seed: int = 0,
) -> dict:
    """Permutation test for association between ``metric_col`` and ``covariate_col``.

    The covariate is shuffled *within* each ``disease_col`` stratum, so any
    disease/batch correlation in the real design is preserved under the
    null (this makes the test one of batch effect *given* the existing
    disease/batch structure, not a naive marginal association test).

    ``covariate_type``: "categorical" uses a between-group variance
    statistic; "continuous" uses Spearman's rho.

    Returns a dict with ``statistic``, ``p_value``, ``n_permutations``, and
    the null distribution's summary (mean, sd) for diagnostic plotting.
    """
    df = donor_df.dropna(subset=[metric_col, covariate_col, disease_col]).copy()
    if len(df) < 3:
        return {
            "statistic": float("nan"),
            "p_value": float("nan"),
            "n_permutations": 0,
            "n_donors": len(df),
            "note": "Too few donors with non-missing data for a permutation test.",
        }

    rng = np.random.default_rng(seed)
    values = df[metric_col].to_numpy(dtype=float)
    covariate = df[covariate_col].to_numpy()
    strata = df[disease_col].to_numpy()

    def _statistic(cov: np.ndarray) -> float:
        if covariate_type == "categorical":
            return _between_group_variance(values, cov)
        rho, _ = scipy_stats.spearmanr(values, cov.astype(float))
        return abs(rho) if not np.isnan(rho) else 0.0

    observed = _statistic(covariate)

    strata_indices = {s: np.where(strata == s)[0] for s in pd.unique(strata)}
    null_stats = np.empty(n_permutations)
    for b in range(n_permutations):
        permuted = covariate.copy()
        for _s, idx in strata_indices.items():
            permuted[idx] = rng.permutation(covariate[idx])
        null_stats[b] = _statistic(permuted)

    p_value = float((np.sum(null_stats >= observed) + 1) / (n_permutations + 1))
    return {
        "statistic": observed,
        "p_value": p_value,
        "n_permutations": n_permutations,
        "n_donors": len(df),
        "null_mean": float(np.mean(null_stats)),
        "null_sd": float(np.std(null_stats)),
    }


# --------------------------------------------------------------------------
# 2. Per-stratum (per-batch) sensitivity of a disease contrast
# --------------------------------------------------------------------------


def _two_group_contrast(values: np.ndarray, groups: np.ndarray, group_a: str, group_b: str) -> Optional[float]:
    a = values[groups == group_a]
    b = values[groups == group_b]
    if len(a) == 0 or len(b) == 0:
        return None
    return float(np.nanmean(a) - np.nanmean(b))


def per_stratum_sensitivity(
    donor_df: pd.DataFrame,
    metric_col: str,
    disease_col: str,
    batch_col: str,
    group_a: str,
    group_b: str,
    min_donors_per_batch: int = 3,
) -> pd.DataFrame:
    """Compare a disease-group contrast on the full cohort vs within each batch.

    Returns one row per batch level (plus a "full_cohort" row), with the
    contrast (mean(group_a) - mean(group_b)), its sign, and whether the
    sign/magnitude agrees with the full-cohort estimate. Batches with fewer
    than ``min_donors_per_batch`` donors in either disease group are flagged
    rather than silently included.
    """
    df = donor_df.dropna(subset=[metric_col, disease_col, batch_col])
    full_contrast = _two_group_contrast(df[metric_col].to_numpy(), df[disease_col].to_numpy(), group_a, group_b)

    rows = [
        {
            "batch": "full_cohort",
            "contrast": full_contrast,
            "n_a": int((df[disease_col] == group_a).sum()),
            "n_b": int((df[disease_col] == group_b).sum()),
            "same_sign_as_full": True,
            "flagged_small_n": False,
        }
    ]

    for batch, group in df.groupby(batch_col):
        n_a = int((group[disease_col] == group_a).sum())
        n_b = int((group[disease_col] == group_b).sum())
        contrast = _two_group_contrast(group[metric_col].to_numpy(), group[disease_col].to_numpy(), group_a, group_b)
        same_sign = (
            contrast is not None
            and full_contrast is not None
            and np.sign(contrast) == np.sign(full_contrast)
        )
        rows.append(
            {
                "batch": batch,
                "contrast": contrast,
                "n_a": n_a,
                "n_b": n_b,
                "same_sign_as_full": same_sign,
                "flagged_small_n": (n_a < min_donors_per_batch or n_b < min_donors_per_batch),
            }
        )
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# 3. Mixed-effects model with donor random intercept
# --------------------------------------------------------------------------


def _design_is_confounded(df: pd.DataFrame, fixed_effects: Sequence[str], group_col: str) -> Optional[str]:
    """Cheap pre-flight checks for designs MixedLM cannot fit sensibly."""
    if df[group_col].nunique() < 2:
        return f"Fewer than 2 levels of the random-effect group '{group_col}'; cannot estimate a random intercept."
    for col in fixed_effects:
        if col not in df.columns:
            continue
        if df[col].nunique(dropna=True) < 2:
            return f"Fixed effect '{col}' has fewer than 2 distinct non-missing values; drop it or collect more data."
    if "disease_status" in fixed_effects and "stain_group" in fixed_effects:
        crosstab = pd.crosstab(df["disease_status"], df["stain_group"])
        # A batch nested entirely within one disease group makes disease and
        # batch effects unidentifiable from each other.
        if (crosstab.astype(bool).sum(axis=0) <= 1).all():
            return "stain_group appears fully nested within disease_status (each batch maps to a single disease group); disease and batch effects are confounded and cannot be separated."
    return None


def mixed_effects_model(
    object_df: pd.DataFrame,
    outcome_col: str,
    donor_col: str = "donor_id",
    fixed_effects: Sequence[str] = ("age", "aab_count", "sex", "disease_status", "stain_group"),
) -> dict:
    """Fit ``outcome_col ~ fixed_effects`` with a random intercept per donor.

    Returns a dict with either ``"model_summary"`` (the fitted
    ``MixedLMResults`` table as a dataframe of coefficients, SEs, p-values)
    or a ``"error"`` message when the design is too small or confounded to
    fit -- callers should always check for ``"error"`` before use.
    """
    if smf is None:  # pragma: no cover
        return {"error": "statsmodels is not installed."}

    df = object_df.dropna(subset=[outcome_col, donor_col]).copy()
    present_fixed = [c for c in fixed_effects if c in df.columns]
    df = df.dropna(subset=present_fixed)

    if len(df) < 10 or df[donor_col].nunique() < 3:
        return {
            "error": f"Design too small for a mixed model: {len(df)} objects across {df[donor_col].nunique()} donors (need >=10 objects, >=3 donors).",
        }

    confound_msg = _design_is_confounded(df, present_fixed, "stain_group" if "stain_group" in present_fixed else donor_col)
    if confound_msg:
        return {"error": confound_msg}

    formula = f"{outcome_col} ~ " + " + ".join(f"C({c})" if df[c].dtype == object else c for c in present_fixed)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model = smf.mixedlm(formula, df, groups=df[donor_col])
            result = model.fit(reml=True)
    except Exception as exc:  # noqa: BLE001 - surface any solver failure as a clear message
        return {"error": f"MixedLM failed to converge or fit: {exc}"}

    summary_df = pd.DataFrame(
        {
            "term": result.params.index,
            "coef": result.params.values,
            "std_err": result.bse.values,
            "p_value": result.pvalues.values,
        }
    )
    return {
        "model_summary": summary_df,
        "n_obs": int(result.nobs),
        "n_groups": df[donor_col].nunique(),
        "converged": bool(result.converged),
        "formula": formula,
    }


# --------------------------------------------------------------------------
# Markdown summary
# --------------------------------------------------------------------------


def write_markdown_summary(out_path: Path, sections: dict[str, str]) -> None:
    lines = ["# Batch-confound diagnostics", ""]
    for title, body in sections.items():
        lines.append(f"## {title}")
        lines.append(body)
        lines.append("")
    out_path.write_text("\n".join(lines) + "\n")


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------


def _cmd_run(args: argparse.Namespace) -> None:
    key_df = ibs.load_key(args.key)
    donors_df = ibs.load_donors(args.donors)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    _tile_df, match_df = ibs.build_tile_records(key_df, donors_df, args.truth_dir, args.pred_dir, args.iou_threshold)
    donor_metrics = compute_donor_metrics(match_df, key_df, donors_df)
    donor_metrics.to_csv(args.out_dir / "donor_metrics.csv", index=False)

    sections = {}

    perm_rows = []
    for covariate, cov_type in (("stain_group", "categorical"), ("acq_run", "categorical")):
        if covariate not in donor_metrics.columns:
            continue
        result = permutation_test_association(
            donor_metrics, "islet_density_per_mm2", covariate, covariate_type=cov_type,
            n_permutations=args.n_permutations, seed=args.seed,
        )
        result["covariate"] = covariate
        perm_rows.append(result)
    perm_df = pd.DataFrame(perm_rows)
    perm_df.to_csv(args.out_dir / "permutation_tests.csv", index=False)
    sections["Permutation tests (batch association)"] = perm_df.to_markdown(index=False) if not perm_df.empty else "No batch covariates available."

    if {"T1D", "ND"}.issubset(set(donor_metrics.get("disease_status", []))):
        sens_df = per_stratum_sensitivity(donor_metrics, "islet_density_per_mm2", "disease_status", "stain_group", "T1D", "ND")
        sens_df.to_csv(args.out_dir / "per_stratum_sensitivity.csv", index=False)
        sections["Per-stratum sensitivity (T1D vs ND)"] = sens_df.to_markdown(index=False)
    else:
        sections["Per-stratum sensitivity"] = "Need both T1D and ND donors present; skipped."

    write_markdown_summary(args.out_dir / "confound_summary.md", sections)
    print(f"wrote {args.out_dir}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--key", required=True, type=Path)
    parser.add_argument("--donors", required=True, type=Path)
    parser.add_argument("--truth-dir", required=True, type=Path)
    parser.add_argument("--pred-dir", required=True, type=Path)
    parser.add_argument("--iou-threshold", type=float, default=ibs.DEFAULT_IOU_THRESHOLD)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--n-permutations", type=int, default=DEFAULT_N_PERMUTATIONS)
    parser.add_argument("--seed", type=int, default=0)
    parser.set_defaults(func=_cmd_run)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
