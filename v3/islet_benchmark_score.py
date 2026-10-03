"""Benchmark scorer for the v3 endocrine-object redefinition (Drotar et al. 2025).

Compares detector output (``pred/<code>.geojson``) against annotator ground
truth (``truth/<code>.geojson``) for a set of benchmark tiles described in
``key.tsv``, using object size classes recomputed from geometry:

    EndocrineSingle  <= 10 um  (mean caliper diameter)
    EndocrineCluster  10-35 um
    Islet            >= 35 um

See ``v3/tests/test_islet_benchmark_score.py`` for worked examples of every
metric, and run ``python islet_benchmark_score.py --help`` for CLI usage.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from shapely.geometry import shape as shapely_shape
from shapely.geometry.base import BaseGeometry

DEFAULT_PIXEL_SIZE_UM = 0.5078
DEFAULT_IOU_THRESHOLD = 0.5
SIZE_CLASSES = ("EndocrineSingle", "EndocrineCluster", "Islet")

STRATUM_DIMENSIONS = ("stain_group", "disease_status", "size_class", "tile_kind")


# --------------------------------------------------------------------------
# Geometry
# --------------------------------------------------------------------------


def qupath_diameters(geom: BaseGeometry) -> tuple[float, float]:
    """``(min, max)`` diameters of ``geom`` exactly as QuPath 0.7 measures "Min diameter" and
    "Max diameter", in the geometry's own units (pixels here).

    Drotar et al. 2025 averaged QuPath's exported Min and Max diameter, and the v3 detector
    classifies with the same QuPath measurements, so the scorer must use QuPath's definitions
    (read from the 0.7 bytecode), not the textbook Feret ones:

    * Max diameter = the diameter of the MINIMUM BOUNDING CIRCLE (2 x its radius). For
      rectangles and right triangles this equals the longest chord; for other shapes it is
      larger (an equilateral triangle: 2/sqrt(3) x its side).
    * Min diameter = the minimum caliper width (JTS MinimumDiameter): over every convex-hull
      edge, the largest perpendicular distance of a hull vertex from that edge's line, minimised.
    """
    import shapely

    hull = geom.convex_hull
    if hull.is_empty or hull.geom_type == "Point":
        return (0.0, 0.0)
    max_d = 2.0 * float(shapely.minimum_bounding_radius(geom))
    if hull.geom_type == "LineString":
        return (0.0, max_d)
    pts = np.asarray(hull.exterior.coords, dtype=float)[:-1]
    min_w = np.inf
    n = len(pts)
    for i in range(n):
        edge = pts[(i + 1) % n] - pts[i]
        length = float(np.hypot(edge[0], edge[1]))
        if length == 0.0:
            continue
        rel = pts - pts[i]
        width = float(np.abs(edge[0] * rel[:, 1] - edge[1] * rel[:, 0]).max() / length)
        min_w = min(min_w, width)
    return (float(min_w) if np.isfinite(min_w) else 0.0, max_d)


def mean_caliper_diameter_um(geom: BaseGeometry, pixel_size_um: float) -> float:
    """Drotar 2025 mean diameter: the mean of QuPath's Min and Max diameter, in micrometres."""
    min_d, max_d = qupath_diameters(geom)
    return ((min_d + max_d) / 2.0) * pixel_size_um


def area_equivalent_diameter_um(geom: BaseGeometry, pixel_size_um: float) -> float:
    """Diameter of a disc with the same area as ``geom``, in micrometres."""
    area_um2 = geom.area * (pixel_size_um ** 2)
    if area_um2 <= 0:
        return 0.0
    return 2.0 * np.sqrt(area_um2 / np.pi)


def classify_size(diameter_um: float) -> str:
    """Classify a mean-caliper (or equivalent) diameter into a size class."""
    if diameter_um <= 10.0:
        return "EndocrineSingle"
    if diameter_um < 35.0:
        return "EndocrineCluster"
    return "Islet"


def compute_iou(a: BaseGeometry, b: BaseGeometry) -> float:
    """Intersection-over-union of two (possibly invalid) geometries."""
    if not a.is_valid:
        a = a.buffer(0)
    if not b.is_valid:
        b = b.buffer(0)
    if a.is_empty or b.is_empty:
        return 0.0
    inter = a.intersection(b).area
    if inter <= 0.0:
        return 0.0
    union = a.area + b.area - inter
    if union <= 0.0:
        return 0.0
    return inter / union


# --------------------------------------------------------------------------
# Data loading
# --------------------------------------------------------------------------


@dataclass
class TileObject:
    """A single annotated/predicted object within one benchmark tile."""

    geometry: BaseGeometry
    mean_caliper_um: float
    equiv_diameter_um: float
    size_class_caliper: str
    size_class_equiv: str
    route: Optional[int] = None
    bridge_flag: Optional[float] = None
    raw_classification: Optional[str] = None


def _measurement_dict(props: dict) -> dict:
    """Normalize QuPath ``measurements`` (dict, or list of {name, value})."""
    m = props.get("measurements")
    if m is None:
        return {}
    if isinstance(m, dict):
        return m
    if isinstance(m, list):
        out = {}
        for item in m:
            if isinstance(item, dict) and "name" in item:
                out[item["name"]] = item.get("value")
        return out
    return {}


def load_objects(geojson_path: Path, pixel_size_um: float) -> list[TileObject]:
    """Load polygons from a QuPath-style GeoJSON, recomputing size class."""
    if not geojson_path.exists():
        return []
    with open(geojson_path) as fh:
        data = json.load(fh)
    features = data.get("features", [])
    objects: list[TileObject] = []
    for feat in features:
        geom_json = feat.get("geometry")
        if geom_json is None:
            continue
        geom = shapely_shape(geom_json)
        if geom.is_empty:
            continue
        props = feat.get("properties", {}) or {}
        measurements = _measurement_dict(props)
        route = measurements.get("Endocrine route")
        bridge_flag = measurements.get("Bridge flag")
        classification = props.get("classification")
        raw_class = classification.get("name") if isinstance(classification, dict) else None
        caliper_um = mean_caliper_diameter_um(geom, pixel_size_um)
        equiv_um = area_equivalent_diameter_um(geom, pixel_size_um)
        objects.append(
            TileObject(
                geometry=geom,
                mean_caliper_um=caliper_um,
                equiv_diameter_um=equiv_um,
                size_class_caliper=classify_size(caliper_um),
                size_class_equiv=classify_size(equiv_um),
                route=int(route) if route is not None else None,
                bridge_flag=float(bridge_flag) if bridge_flag is not None else None,
                raw_classification=raw_class,
            )
        )
    return objects


def load_key(key_path: Path) -> pd.DataFrame:
    df = pd.read_csv(key_path, sep="\t", dtype={"code": str, "image": str})
    if "pixel_size_um" not in df.columns:
        df["pixel_size_um"] = DEFAULT_PIXEL_SIZE_UM
    else:
        df["pixel_size_um"] = df["pixel_size_um"].fillna(DEFAULT_PIXEL_SIZE_UM)
    return df


def load_donors(donors_path: Path) -> pd.DataFrame:
    df = pd.read_csv(donors_path, sep="\t", dtype={"image": str, "donor_id": str})
    df["disease_status"] = df["disease_status"].replace({"AAB": "Aab+"})
    return df


# --------------------------------------------------------------------------
# Matching and per-tile scoring
# --------------------------------------------------------------------------


@dataclass
class TileScore:
    code: str
    tp: int
    fp: int
    fn: int
    matches: list[tuple[int, int, float]] = field(default_factory=list)  # (truth_idx, pred_idx, iou)
    unmatched_pred: list[int] = field(default_factory=list)
    unmatched_truth: list[int] = field(default_factory=list)


def match_tile(
    truth_objects: Sequence[TileObject],
    pred_objects: Sequence[TileObject],
    iou_threshold: float = DEFAULT_IOU_THRESHOLD,
) -> TileScore:
    """Match predictions to truth within one tile by optimal (max total) IoU."""
    n_truth, n_pred = len(truth_objects), len(pred_objects)
    if n_truth == 0 or n_pred == 0:
        return TileScore(
            code="",
            tp=0,
            fp=n_pred,
            fn=n_truth,
            unmatched_pred=list(range(n_pred)),
            unmatched_truth=list(range(n_truth)),
        )

    iou_matrix = np.zeros((n_truth, n_pred))
    for i, t in enumerate(truth_objects):
        for j, p in enumerate(pred_objects):
            iou_matrix[i, j] = compute_iou(t.geometry, p.geometry)

    cost_matrix = 1.0 - iou_matrix
    row_idx, col_idx = linear_sum_assignment(cost_matrix)

    matches: list[tuple[int, int, float]] = []
    matched_truth, matched_pred = set(), set()
    for i, j in zip(row_idx, col_idx):
        iou = iou_matrix[i, j]
        if iou >= iou_threshold:
            matches.append((int(i), int(j), float(iou)))
            matched_truth.add(int(i))
            matched_pred.add(int(j))

    unmatched_truth = [i for i in range(n_truth) if i not in matched_truth]
    unmatched_pred = [j for j in range(n_pred) if j not in matched_pred]

    return TileScore(
        code="",
        tp=len(matches),
        fp=len(unmatched_pred),
        fn=len(unmatched_truth),
        matches=matches,
        unmatched_pred=unmatched_pred,
        unmatched_truth=unmatched_truth,
    )


# --------------------------------------------------------------------------
# Building the tidy record tables
# --------------------------------------------------------------------------


def build_tile_records(
    key_df: pd.DataFrame,
    donors_df: pd.DataFrame,
    truth_dir: Path,
    pred_dir: Path,
    iou_threshold: float = DEFAULT_IOU_THRESHOLD,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Score every tile in ``key_df``.

    Returns ``(tile_df, match_df)``:
      - ``tile_df``: one row per tile with tp/fp/fn counts and stratum columns.
      - ``match_df``: one row per matched (truth, pred) pair, with both size
        classes (caliper- and equivalent-diameter-based) and the predicted route,
        for confusion-matrix and per-route analyses.
    """
    donor_by_image = donors_df.set_index("image")
    tile_rows = []
    match_rows = []

    for _, row in key_df.iterrows():
        code = row["code"]
        image = row["image"]
        pixel_size_um = float(row["pixel_size_um"])
        truth_objs = load_objects(truth_dir / f"{code}.geojson", pixel_size_um)
        pred_objs = load_objects(pred_dir / f"{code}.geojson", pixel_size_um)

        score = match_tile(truth_objs, pred_objs, iou_threshold)
        score.code = code

        donor_info = donor_by_image.loc[image] if image in donor_by_image.index else None
        donor_id = donor_info["donor_id"] if donor_info is not None else None
        disease_status = donor_info["disease_status"] if donor_info is not None else None
        stain_group = donor_info["stain_group"] if donor_info is not None else None

        tile_rows.append(
            {
                "code": code,
                "image": image,
                "donor_id": donor_id,
                "stain_group": stain_group,
                "disease_status": disease_status,
                "tile_kind": row.get("tile_kind"),
                "tp": score.tp,
                "fp": score.fp,
                "fn": score.fn,
                "n_truth": len(truth_objs),
                "n_pred": len(pred_objs),
            }
        )

        base = {
            "code": code,
            "image": image,
            "donor_id": donor_id,
            "stain_group": stain_group,
            "disease_status": disease_status,
            "tile_kind": row.get("tile_kind"),
        }

        for t_idx, p_idx, iou in score.matches:
            t_obj, p_obj = truth_objs[t_idx], pred_objs[p_idx]
            match_rows.append(
                {
                    **base,
                    "outcome": "tp",
                    "iou": iou,
                    "truth_size_class_caliper": t_obj.size_class_caliper,
                    "pred_size_class_caliper": p_obj.size_class_caliper,
                    "truth_size_class_equiv": t_obj.size_class_equiv,
                    "pred_size_class_equiv": p_obj.size_class_equiv,
                    "pred_route": p_obj.route,
                    "size_class": p_obj.size_class_caliper,
                }
            )

        # Unmatched predictions carry a size class / route for per-route precision.
        for p_idx in score.unmatched_pred:
            p_obj = pred_objs[p_idx]
            match_rows.append(
                {
                    **base,
                    "outcome": "fp",
                    "iou": 0.0,
                    "truth_size_class_caliper": None,
                    "pred_size_class_caliper": p_obj.size_class_caliper,
                    "truth_size_class_equiv": None,
                    "pred_size_class_equiv": p_obj.size_class_equiv,
                    "pred_route": p_obj.route,
                    "size_class": p_obj.size_class_caliper,
                }
            )

        # Unmatched truth carries a size class for FN attribution / recall-by-size.
        for t_idx in score.unmatched_truth:
            t_obj = truth_objs[t_idx]
            match_rows.append(
                {
                    **base,
                    "outcome": "fn",
                    "iou": 0.0,
                    "truth_size_class_caliper": t_obj.size_class_caliper,
                    "pred_size_class_caliper": None,
                    "truth_size_class_equiv": t_obj.size_class_equiv,
                    "pred_size_class_equiv": None,
                    "pred_route": None,
                    "size_class": t_obj.size_class_caliper,
                }
            )

    tile_df = pd.DataFrame(tile_rows)
    match_df = pd.DataFrame(match_rows)
    if "outcome" not in match_df.columns:
        match_df["outcome"] = pd.Series(dtype=str)
    if "size_class" not in tile_df.columns:
        tile_df["size_class"] = None
    return tile_df, match_df


# --------------------------------------------------------------------------
# Aggregate metrics
# --------------------------------------------------------------------------


def precision_recall_f1(tp: int, fp: int, fn: int) -> tuple[float, float, float]:
    """Precision/recall/F1, with the zero-division convention used for ranking.

    A stratum with no truth and no predictions (tp=fp=fn=0) is inapplicable
    and returns all-NaN so it is skipped by macro-averaging. Otherwise a
    missing denominator (no predictions, or no truth objects) counts as 0,
    following scikit-learn's ``zero_division=0`` convention -- e.g. a
    detector that predicts nothing when truth objects exist gets F1=0
    rather than an undefined score that would be silently dropped.
    """
    if tp + fp + fn == 0:
        return float("nan"), float("nan"), float("nan")
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 0.0 if (precision + recall) == 0 else 2 * precision * recall / (precision + recall)
    return precision, recall, f1


def aggregate_by_stratum(
    match_df: pd.DataFrame,
    dimensions: Sequence[str] = STRATUM_DIMENSIONS,
) -> pd.DataFrame:
    """Precision/recall/F1 per stratum, one stratum = one distinct combo of ``dimensions``.

    Operates on the object-level ``match_df`` produced by ``build_tile_records``
    (every TP/FP/FN object is one row, already carrying stain_group,
    disease_status, tile_kind and size_class), so any subset of
    ``STRATUM_DIMENSIONS`` -- including ``size_class`` -- can be used. Pass
    ``dimensions=()`` for a single overall row.
    """
    all_dims = list(STRATUM_DIMENSIONS)
    dims = [d for d in dimensions if d in all_dims]
    if match_df.empty:
        cols = list(dims) + ["tp", "fp", "fn", "precision", "recall", "f1"]
        return pd.DataFrame(columns=cols)

    def _count(outcome: str) -> pd.Series:
        subset = match_df[match_df["outcome"] == outcome]
        if subset.empty or not dims:
            return pd.Series(dtype=int)
        return subset.groupby(dims, dropna=False).size()

    if not dims:
        tp = int((match_df["outcome"] == "tp").sum())
        fp = int((match_df["outcome"] == "fp").sum())
        fn = int((match_df["outcome"] == "fn").sum())
        precision, recall, f1 = precision_recall_f1(tp, fp, fn)
        return pd.DataFrame([{"tp": tp, "fp": fp, "fn": fn, "precision": precision, "recall": recall, "f1": f1}])

    tp_counts, fp_counts, fn_counts = _count("tp"), _count("fp"), _count("fn")
    all_keys = set(tp_counts.index) | set(fp_counts.index) | set(fn_counts.index)
    rows = []
    for key in all_keys:
        keys = key if isinstance(key, tuple) else (key,)
        tp = int(tp_counts.get(key, 0))
        fp = int(fp_counts.get(key, 0))
        fn = int(fn_counts.get(key, 0))
        precision, recall, f1 = precision_recall_f1(tp, fp, fn)
        row = dict(zip(dims, keys))
        row.update({"tp": tp, "fp": fp, "fn": fn, "precision": precision, "recall": recall, "f1": f1})
        rows.append(row)
    return pd.DataFrame(rows)


def macro_f1(stratum_df: pd.DataFrame) -> float:
    """Macro-average F1 across strata (NaN rows, i.e. empty strata, are skipped)."""
    if stratum_df.empty or "f1" not in stratum_df.columns:
        return float("nan")
    valid = stratum_df["f1"].dropna()
    return float(valid.mean()) if len(valid) else float("nan")


def per_route_precision(match_df: pd.DataFrame) -> pd.DataFrame:
    """Precision of predictions broken down by 'Endocrine route' (1=CD99, 2=hormone, 3=both)."""
    predictions = match_df[match_df["outcome"] != "fn"] if not match_df.empty else match_df
    if predictions.empty:
        return pd.DataFrame(columns=["pred_route", "tp", "fp", "precision"])
    rows = []
    for route, group in predictions.groupby("pred_route", dropna=False):
        tp = int((group["outcome"] == "tp").sum())
        fp = int((group["outcome"] == "fp").sum())
        precision = tp / (tp + fp) if (tp + fp) > 0 else float("nan")
        rows.append({"pred_route": route, "tp": tp, "fp": fp, "precision": precision})
    return pd.DataFrame(rows)


def size_class_confusion(match_df: pd.DataFrame, kind: str = "caliper") -> pd.DataFrame:
    """Truth-vs-prediction size-class confusion matrix over matched (TP) pairs.

    ``kind`` is ``"caliper"`` or ``"equiv"`` (area-equivalent diameter).
    """
    col_truth = f"truth_size_class_{kind}"
    col_pred = f"pred_size_class_{kind}"
    matched = match_df[match_df["outcome"] == "tp"] if not match_df.empty else match_df
    if matched.empty:
        return pd.DataFrame(0, index=list(SIZE_CLASSES), columns=list(SIZE_CLASSES))
    confusion = pd.crosstab(matched[col_truth], matched[col_pred])
    confusion = confusion.reindex(index=SIZE_CLASSES, columns=SIZE_CLASSES, fill_value=0)
    return confusion


# --------------------------------------------------------------------------
# Donor-block bootstrap
# --------------------------------------------------------------------------


def donor_block_bootstrap(
    match_df: pd.DataFrame,
    dimensions: Sequence[str] = STRATUM_DIMENSIONS,
    n_boot: int = 500,
    seed: int = 0,
) -> pd.DataFrame:
    """Bootstrap CIs for macro-F1 and per-stratum F1 by resampling donors with replacement.

    ``match_df`` is the object-level table from ``build_tile_records``.
    Returns a tidy dataframe with one row per stratum (plus one "__macro__"
    row) giving the point estimate and (2.5, 97.5) percentile CI.
    """
    rng = np.random.default_rng(seed)
    donors = match_df["donor_id"].dropna().unique()
    if len(donors) == 0:
        return pd.DataFrame(columns=["stratum", "estimate", "ci_low", "ci_high"])

    point_strata = aggregate_by_stratum(match_df, dimensions)
    point_macro = macro_f1(point_strata)

    boot_macro = []
    boot_strata: dict[tuple, list[float]] = {}
    dims = [d for d in dimensions if d in STRATUM_DIMENSIONS]

    for _ in range(n_boot):
        sample_donors = rng.choice(donors, size=len(donors), replace=True)
        parts = [match_df[match_df["donor_id"] == d] for d in sample_donors]
        boot_match = pd.concat(parts, ignore_index=True) if parts else match_df.iloc[0:0]
        strata = aggregate_by_stratum(boot_match, dimensions)
        boot_macro.append(macro_f1(strata))
        for _, row in strata.iterrows():
            key = tuple(row[d] for d in dims)
            boot_strata.setdefault(key, []).append(row["f1"])

    rows = [
        {
            "stratum": "__macro__",
            "estimate": point_macro,
            "ci_low": float(np.nanpercentile(boot_macro, 2.5)) if boot_macro else float("nan"),
            "ci_high": float(np.nanpercentile(boot_macro, 97.5)) if boot_macro else float("nan"),
        }
    ]
    dims = [d for d in dimensions if d in point_strata.columns]
    for _, row in point_strata.iterrows():
        key = tuple(row[d] for d in dims)
        vals = [v for v in boot_strata.get(key, []) if not (isinstance(v, float) and np.isnan(v))]
        rows.append(
            {
                "stratum": "|".join(str(k) for k in key) if key else "all",
                "estimate": row["f1"],
                "ci_low": float(np.nanpercentile(vals, 2.5)) if vals else float("nan"),
                "ci_high": float(np.nanpercentile(vals, 97.5)) if vals else float("nan"),
            }
        )
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# Parameter sweep
# --------------------------------------------------------------------------


def sweep_parameter_sets(
    key_df: pd.DataFrame,
    donors_df: pd.DataFrame,
    truth_dir: Path,
    pred_dirs: dict[str, Path],
    iou_threshold: float = DEFAULT_IOU_THRESHOLD,
    dimensions: Sequence[str] = STRATUM_DIMENSIONS,
) -> pd.DataFrame:
    """Score several detector parameter sets and rank them by macro-F1.

    ``pred_dirs`` maps a parameter-set label (e.g. ``"tauH3_tauC3"``) to its
    predictions directory.
    """
    rows = []
    for label, pred_dir in pred_dirs.items():
        _tile_df, match_df = build_tile_records(key_df, donors_df, truth_dir, pred_dir, iou_threshold)
        strata = aggregate_by_stratum(match_df, dimensions)
        m_f1 = macro_f1(strata)
        if not strata.empty and strata["f1"].notna().any():
            worst = strata.loc[strata["f1"].idxmin()]
            worst_desc = ", ".join(f"{d}={worst[d]}" for d in dimensions if d in worst.index)
            worst_f1 = worst["f1"]
        else:
            worst_desc, worst_f1 = "n/a", float("nan")
        rows.append({"parameter_set": label, "macro_f1": m_f1, "worst_stratum": worst_desc, "worst_stratum_f1": worst_f1})
    result = pd.DataFrame(rows).sort_values("macro_f1", ascending=False, na_position="last").reset_index(drop=True)
    return result


# --------------------------------------------------------------------------
# Markdown summary
# --------------------------------------------------------------------------


def write_markdown_summary(
    out_path: Path,
    tile_df: pd.DataFrame,
    stratum_df: pd.DataFrame,
    macro_f1_value: float,
    route_df: pd.DataFrame,
    confusion_caliper: pd.DataFrame,
    confusion_equiv: pd.DataFrame,
) -> None:
    lines = ["# Islet benchmark summary", ""]
    tp, fp, fn = int(tile_df["tp"].sum()), int(tile_df["fp"].sum()), int(tile_df["fn"].sum())
    precision, recall, f1 = precision_recall_f1(tp, fp, fn)
    lines.append(f"- Tiles scored: {len(tile_df)}")
    lines.append(f"- Overall: precision={precision:.3f}, recall={recall:.3f}, F1={f1:.3f}")
    lines.append(f"- **Macro-F1 across strata (selection objective): {macro_f1_value:.3f}**")
    lines.append("")
    lines.append("## Per-stratum metrics")
    lines.append(stratum_df.to_markdown(index=False))
    lines.append("")
    lines.append("## Per-route precision")
    lines.append(route_df.to_markdown(index=False))
    lines.append("")
    lines.append("## Size-class confusion (caliper diameter)")
    lines.append(confusion_caliper.to_markdown())
    lines.append("")
    lines.append("## Size-class confusion (area-equivalent diameter)")
    lines.append(confusion_equiv.to_markdown())
    out_path.write_text("\n".join(lines) + "\n")


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------


def _add_common_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--key", required=True, type=Path, help="key.tsv")
    parser.add_argument("--donors", required=True, type=Path, help="donors.tsv")
    parser.add_argument("--truth-dir", required=True, type=Path, help="Directory of truth/<code>.geojson")
    parser.add_argument("--iou-threshold", type=float, default=DEFAULT_IOU_THRESHOLD)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument(
        "--stratum-dims",
        nargs="*",
        default=list(STRATUM_DIMENSIONS),
        help=f"Subset of {STRATUM_DIMENSIONS} to stratify by",
    )
    parser.add_argument("--bootstrap", type=int, default=0, help="Number of donor-block bootstrap resamples (0=skip)")
    parser.add_argument("--seed", type=int, default=0)


def _cmd_score(args: argparse.Namespace) -> None:
    key_df = load_key(args.key)
    donors_df = load_donors(args.donors)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    tile_df, match_df = build_tile_records(key_df, donors_df, args.truth_dir, args.pred_dir, args.iou_threshold)
    stratum_df = aggregate_by_stratum(match_df, args.stratum_dims)
    m_f1 = macro_f1(stratum_df)
    route_df = per_route_precision(match_df)
    confusion_caliper = size_class_confusion(match_df, "caliper")
    confusion_equiv = size_class_confusion(match_df, "equiv")

    tile_df.to_csv(args.out_dir / "tile_scores.csv", index=False)
    match_df.to_csv(args.out_dir / "matches.csv", index=False)
    stratum_df.to_csv(args.out_dir / "stratum_metrics.csv", index=False)
    route_df.to_csv(args.out_dir / "route_precision.csv", index=False)
    confusion_caliper.to_csv(args.out_dir / "size_confusion_caliper.csv")
    confusion_equiv.to_csv(args.out_dir / "size_confusion_equiv.csv")

    if args.bootstrap > 0:
        boot_df = donor_block_bootstrap(match_df, args.stratum_dims, args.bootstrap, args.seed)
        boot_df.to_csv(args.out_dir / "bootstrap_ci.csv", index=False)

    write_markdown_summary(
        args.out_dir / "summary.md", tile_df, stratum_df, m_f1, route_df, confusion_caliper, confusion_equiv
    )
    print(f"macro_f1={m_f1:.4f}  tiles={len(tile_df)}  wrote {args.out_dir}")


def _cmd_sweep(args: argparse.Namespace) -> None:
    key_df = load_key(args.key)
    donors_df = load_donors(args.donors)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    pred_dirs = {}
    for spec in args.pred_dir:
        if "=" in spec:
            label, path = spec.split("=", 1)
        else:
            path = spec
            label = Path(spec).name
        pred_dirs[label] = Path(path)

    ranking = sweep_parameter_sets(key_df, donors_df, args.truth_dir, pred_dirs, args.iou_threshold, args.stratum_dims)
    ranking.to_csv(args.out_dir / "sweep_ranking.csv", index=False)
    print(ranking.to_string(index=False))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    p_score = sub.add_parser("score", help="Score one prediction directory against truth")
    _add_common_args(p_score)
    p_score.add_argument("--pred-dir", required=True, type=Path)
    p_score.set_defaults(func=_cmd_score)

    p_sweep = sub.add_parser("sweep", help="Rank several prediction directories (parameter sets) by macro-F1")
    _add_common_args(p_sweep)
    p_sweep.add_argument(
        "--pred-dir",
        required=True,
        nargs="+",
        help="One or more [label=]path entries, e.g. tauH3_tauC3=pred_tauH3_tauC3/",
    )
    p_sweep.set_defaults(func=_cmd_sweep)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
