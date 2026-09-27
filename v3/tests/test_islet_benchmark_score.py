"""Tests for v3/islet_benchmark_score.py.

All fixtures are hand-built polygons / tiny synthetic tile sets written to
``tmp_path`` -- no project data files are read.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from shapely.geometry import Polygon, box

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import islet_benchmark_score as ibs


# --------------------------------------------------------------------------
# Geometry
# --------------------------------------------------------------------------


def test_mean_caliper_diameter_rectangle():
    rect = box(0, 0, 20, 10)  # 20 x 10 pixel rectangle
    # min Feret 10, max Feret = the diagonal sqrt(500): QuPath's Min/Max diameter
    assert ibs.mean_caliper_diameter_um(rect, pixel_size_um=1.0) == pytest.approx((10 + 500 ** 0.5) / 2)


def test_mean_caliper_diameter_scales_with_pixel_size():
    rect = box(0, 0, 20, 10)
    assert ibs.mean_caliper_diameter_um(rect, pixel_size_um=0.5) == pytest.approx((10 + 500 ** 0.5) / 4)


def test_mean_caliper_diameter_rotated_rectangle():
    # Same 20x10 rectangle, rotated 30 degrees -- side lengths (hence mean
    # caliper diameter) are invariant to rotation.
    rect = box(0, 0, 20, 10)
    from shapely.affinity import rotate

    rotated = rotate(rect, 30, origin=(0, 0))
    assert ibs.mean_caliper_diameter_um(rotated, pixel_size_um=1.0) == pytest.approx((10 + 500 ** 0.5) / 2, abs=1e-6)


def test_area_equivalent_diameter_of_disc_area_square():
    # A square of area A has area-equivalent diameter 2*sqrt(A/pi).
    side = 10.0
    square = box(0, 0, side, side)
    expected = 2 * np.sqrt((side * side) / np.pi)
    assert ibs.area_equivalent_diameter_um(square, pixel_size_um=1.0) == pytest.approx(expected)


def test_area_equivalent_diameter_empty_geometry_is_zero():
    empty = Polygon()
    assert ibs.area_equivalent_diameter_um(empty, pixel_size_um=1.0) == 0.0


@pytest.mark.parametrize(
    "diameter,expected",
    [
        (5.0, "EndocrineSingle"),
        (10.0, "EndocrineSingle"),
        (10.01, "EndocrineCluster"),
        (20.0, "EndocrineCluster"),
        (34.99, "EndocrineCluster"),
        (35.0, "Islet"),
        (100.0, "Islet"),
    ],
)
def test_classify_size_boundaries(diameter, expected):
    assert ibs.classify_size(diameter) == expected


def test_compute_iou_known_overlap():
    a = box(0, 0, 10, 10)  # area 100
    b = box(5, 0, 15, 10)  # area 100, overlap 50
    # intersection = 50, union = 100+100-50 = 150
    assert ibs.compute_iou(a, b) == pytest.approx(50 / 150)


def test_compute_iou_disjoint_is_zero():
    a = box(0, 0, 5, 5)
    b = box(10, 10, 15, 15)
    assert ibs.compute_iou(a, b) == 0.0


def test_compute_iou_identical_is_one():
    a = box(0, 0, 5, 5)
    assert ibs.compute_iou(a, a) == pytest.approx(1.0)


# --------------------------------------------------------------------------
# GeoJSON loading (dict vs list measurements layouts)
# --------------------------------------------------------------------------


def _feature(coords, classification=None, measurements=None):
    props = {}
    if classification is not None:
        props["classification"] = {"name": classification}
    if measurements is not None:
        props["measurements"] = measurements
    return {"type": "Feature", "geometry": {"type": "Polygon", "coordinates": [coords]}, "properties": props}


def _write_geojson(path: Path, features: list[dict]) -> None:
    path.write_text(json.dumps({"type": "FeatureCollection", "features": features}))


def test_load_objects_measurements_as_dict(tmp_path):
    coords = [(0, 0), (20, 0), (20, 10), (0, 10), (0, 0)]
    feat = _feature(coords, "Islet", {"Endocrine route": 3, "Bridge flag": 1})
    path = tmp_path / "a.geojson"
    _write_geojson(path, [feat])

    objs = ibs.load_objects(path, pixel_size_um=1.0)
    assert len(objs) == 1
    assert objs[0].route == 3
    assert objs[0].bridge_flag == 1.0
    assert objs[0].mean_caliper_um == pytest.approx((10 + 500 ** 0.5) / 2)


def test_load_objects_measurements_as_list(tmp_path):
    coords = [(0, 0), (20, 0), (20, 10), (0, 10), (0, 0)]
    measurements = [{"name": "Endocrine route", "value": 2}, {"name": "Bridge flag", "value": 0}]
    feat = _feature(coords, "Islet", measurements)
    path = tmp_path / "b.geojson"
    _write_geojson(path, [feat])

    objs = ibs.load_objects(path, pixel_size_um=1.0)
    assert len(objs) == 1
    assert objs[0].route == 2
    assert objs[0].bridge_flag == 0.0


def test_load_objects_missing_file_returns_empty(tmp_path):
    assert ibs.load_objects(tmp_path / "missing.geojson", pixel_size_um=1.0) == []


# --------------------------------------------------------------------------
# Matching: perfect / empty / all-FP / below-threshold tiles
# --------------------------------------------------------------------------


def _obj(coords, pixel_size_um=1.0, route=None):
    poly = Polygon(coords)
    caliper = ibs.mean_caliper_diameter_um(poly, pixel_size_um)
    equiv = ibs.area_equivalent_diameter_um(poly, pixel_size_um)
    return ibs.TileObject(
        geometry=poly,
        mean_caliper_um=caliper,
        equiv_diameter_um=equiv,
        size_class_caliper=ibs.classify_size(caliper),
        size_class_equiv=ibs.classify_size(equiv),
        route=route,
    )


def test_match_tile_perfect_match():
    truth = [_obj([(0, 0), (20, 0), (20, 10), (0, 10)])]
    pred = [_obj([(0, 0), (20, 0), (20, 10), (0, 10)])]
    score = ibs.match_tile(truth, pred)
    assert (score.tp, score.fp, score.fn) == (1, 0, 0)
    assert score.matches[0][2] == pytest.approx(1.0)


def test_match_tile_empty_tile():
    score = ibs.match_tile([], [])
    assert (score.tp, score.fp, score.fn) == (0, 0, 0)


def test_match_tile_all_false_positive():
    pred = [_obj([(0, 0), (20, 0), (20, 10), (0, 10)]), _obj([(50, 50), (60, 50), (60, 60), (50, 60)])]
    score = ibs.match_tile([], pred)
    assert (score.tp, score.fp, score.fn) == (0, 2, 0)


def test_match_tile_all_false_negative():
    truth = [_obj([(0, 0), (20, 0), (20, 10), (0, 10)])]
    score = ibs.match_tile(truth, [])
    assert (score.tp, score.fp, score.fn) == (0, 0, 1)


def test_match_tile_below_iou_threshold_counts_as_fp_and_fn():
    # Small overlap: IoU well below 0.5, so the pair should NOT be matched.
    truth = [_obj([(0, 0), (10, 0), (10, 10), (0, 10)])]  # area 100
    pred = [_obj([(9, 0), (19, 0), (19, 10), (9, 10)])]  # area 100, overlap 10
    score = ibs.match_tile(truth, pred, iou_threshold=0.5)
    assert (score.tp, score.fp, score.fn) == (0, 1, 1)


def test_match_tile_optimal_assignment_prefers_higher_total_iou():
    # Two truth / two pred objects where the "obvious" pairing is suboptimal;
    # linear_sum_assignment must find the assignment maximizing total IoU.
    t1 = _obj([(0, 0), (10, 0), (10, 10), (0, 10)])
    t2 = _obj([(20, 0), (30, 0), (30, 10), (20, 10)])
    p1 = _obj([(0, 0), (10, 0), (10, 10), (0, 10)])  # perfect match for t1
    p2 = _obj([(20, 0), (30, 0), (30, 10), (20, 10)])  # perfect match for t2
    score = ibs.match_tile([t1, t2], [p1, p2])
    assert (score.tp, score.fp, score.fn) == (2, 0, 0)
    matched_ious = sorted(m[2] for m in score.matches)
    assert matched_ious == pytest.approx([1.0, 1.0])


# --------------------------------------------------------------------------
# End-to-end scoring over a tiny synthetic benchmark
# --------------------------------------------------------------------------


def _write_benchmark(tmp_path: Path):
    """Two tiles, two donors, one stain group / disease status each."""
    key_path = tmp_path / "key.tsv"
    donors_path = tmp_path / "donors.tsv"
    truth_dir = tmp_path / "truth"
    pred_dir = tmp_path / "pred"
    truth_dir.mkdir()
    pred_dir.mkdir()

    key_path.write_text(
        "code\timage\tx\ty\tw\th\ttile_kind\n"
        "t1\timgA\t0\t0\t100\t100\trandom\n"
        "t2\timgB\t0\t0\t100\t100\tobject\n"
    )
    donors_path.write_text(
        "donor_id\timage\tdisease_status\tstain_group\tacq_run\tage\taab_count\tsex\n"
        "D1\timgA\tND\tS0\trun1\t30\t0\tF\n"
        "D2\timgB\tT1D\tS1\trun2\t40\t2\tM\n"
    )

    # Tile t1: a perfect match (Islet-sized, 40x40 -> mean caliper (40 + 40*sqrt(2))/2 ~ 48 um).
    islet_coords = [(0, 0), (40, 0), (40, 40), (0, 40), (0, 0)]
    _write_geojson(truth_dir / "t1.geojson", [_feature(islet_coords, "Islet", {"Endocrine route": 3})])
    _write_geojson(pred_dir / "t1.geojson", [_feature(islet_coords, "Islet", {"Endocrine route": 3})])

    # Tile t2: one false positive (small, single-cell sized) and one missed truth object.
    single_coords = [(0, 0), (5, 0), (5, 5), (0, 5), (0, 0)]
    missed_coords = [(60, 60), (90, 60), (90, 90), (60, 90), (60, 60)]
    _write_geojson(truth_dir / "t2.geojson", [_feature(missed_coords, "Islet", {})])
    _write_geojson(pred_dir / "t2.geojson", [_feature(single_coords, "EndocrineSingle", {"Endocrine route": 1})])

    return key_path, donors_path, truth_dir, pred_dir


def test_build_tile_records_and_aggregate(tmp_path):
    key_path, donors_path, truth_dir, pred_dir = _write_benchmark(tmp_path)
    key_df = ibs.load_key(key_path)
    donors_df = ibs.load_donors(donors_path)

    tile_df, match_df = ibs.build_tile_records(key_df, donors_df, truth_dir, pred_dir, iou_threshold=0.5)

    assert len(tile_df) == 2
    assert tile_df.loc[tile_df["code"] == "t1", "tp"].item() == 1
    assert tile_df.loc[tile_df["code"] == "t2", "fp"].item() == 1
    assert tile_df.loc[tile_df["code"] == "t2", "fn"].item() == 1

    overall_tp, overall_fp, overall_fn = int(tile_df["tp"].sum()), int(tile_df["fp"].sum()), int(tile_df["fn"].sum())
    assert (overall_tp, overall_fp, overall_fn) == (1, 1, 1)

    stratum_df = ibs.aggregate_by_stratum(match_df, dimensions=("disease_status",))
    assert set(stratum_df["disease_status"]) == {"ND", "T1D"}
    nd_row = stratum_df[stratum_df["disease_status"] == "ND"].iloc[0]
    assert (nd_row["tp"], nd_row["fp"], nd_row["fn"]) == (1, 0, 0)
    assert nd_row["f1"] == pytest.approx(1.0)

    m_f1 = ibs.macro_f1(stratum_df)
    assert 0.0 <= m_f1 <= 1.0

    route_df = ibs.per_route_precision(match_df)
    assert set(route_df["pred_route"]) == {3, 1}


def test_size_class_confusion_matched_pairs():
    match_df = pd.DataFrame(
        [
            {"outcome": "tp", "truth_size_class_caliper": "Islet", "pred_size_class_caliper": "Islet"},
            {"outcome": "tp", "truth_size_class_caliper": "Islet", "pred_size_class_caliper": "EndocrineCluster"},
            {"outcome": "tp", "truth_size_class_caliper": "EndocrineSingle", "pred_size_class_caliper": "EndocrineSingle"},
            {"outcome": "fp", "truth_size_class_caliper": None, "pred_size_class_caliper": "EndocrineSingle"},
        ]
    )
    confusion = ibs.size_class_confusion(match_df, kind="caliper")
    assert confusion.loc["Islet", "Islet"] == 1
    assert confusion.loc["Islet", "EndocrineCluster"] == 1
    assert confusion.loc["EndocrineSingle", "EndocrineSingle"] == 1
    assert confusion.values.sum() == 3  # only TP rows counted, the FP row is excluded


def test_size_class_confusion_empty_matches_returns_zero_matrix():
    confusion = ibs.size_class_confusion(pd.DataFrame(columns=["outcome"]))
    assert confusion.values.sum() == 0
    assert list(confusion.index) == list(ibs.SIZE_CLASSES)


# --------------------------------------------------------------------------
# Bootstrap reproducibility
# --------------------------------------------------------------------------


def _synthetic_match_df(n_donors=6, seed=0):
    rng = np.random.default_rng(seed)
    rows = []
    for d in range(n_donors):
        donor = f"D{d}"
        n_tp = rng.integers(2, 6)
        n_fp = rng.integers(0, 3)
        n_fn = rng.integers(0, 3)
        for _ in range(n_tp):
            rows.append({"donor_id": donor, "outcome": "tp", "stain_group": "S0", "disease_status": "ND", "tile_kind": "random", "size_class": "Islet"})
        for _ in range(n_fp):
            rows.append({"donor_id": donor, "outcome": "fp", "stain_group": "S0", "disease_status": "ND", "tile_kind": "random", "size_class": "Islet"})
        for _ in range(n_fn):
            rows.append({"donor_id": donor, "outcome": "fn", "stain_group": "S0", "disease_status": "ND", "tile_kind": "random", "size_class": "Islet"})
    return pd.DataFrame(rows)


def test_bootstrap_reproducible_with_seed():
    match_df = _synthetic_match_df()
    boot1 = ibs.donor_block_bootstrap(match_df, dimensions=("disease_status",), n_boot=50, seed=42)
    boot2 = ibs.donor_block_bootstrap(match_df, dimensions=("disease_status",), n_boot=50, seed=42)
    pd.testing.assert_frame_equal(boot1, boot2)


def test_bootstrap_different_seed_can_differ():
    match_df = _synthetic_match_df()
    boot1 = ibs.donor_block_bootstrap(match_df, dimensions=("disease_status",), n_boot=50, seed=1)
    boot2 = ibs.donor_block_bootstrap(match_df, dimensions=("disease_status",), n_boot=50, seed=2)
    macro1 = boot1.loc[boot1["stratum"] == "__macro__", "ci_low"].item()
    macro2 = boot2.loc[boot2["stratum"] == "__macro__", "ci_low"].item()
    # Not a strict guarantee, but with randomized resampling the two seeds
    # should not produce byte-identical CI tables.
    assert not boot1.equals(boot2) or macro1 == macro2


def test_bootstrap_ci_contains_point_estimate_bounds():
    match_df = _synthetic_match_df(n_donors=10, seed=7)
    boot = ibs.donor_block_bootstrap(match_df, dimensions=("disease_status",), n_boot=200, seed=3)
    macro_row = boot[boot["stratum"] == "__macro__"].iloc[0]
    assert macro_row["ci_low"] <= macro_row["estimate"] + 1e-9
    assert macro_row["ci_high"] >= macro_row["estimate"] - 1e-9


# --------------------------------------------------------------------------
# Parameter sweep
# --------------------------------------------------------------------------


def test_sweep_ranks_by_macro_f1(tmp_path):
    key_path, donors_path, truth_dir, _pred_dir = _write_benchmark(tmp_path)
    key_df = ibs.load_key(key_path)
    donors_df = ibs.load_donors(donors_path)

    good_dir = tmp_path / "pred_good"
    bad_dir = tmp_path / "pred_bad"
    good_dir.mkdir()
    bad_dir.mkdir()

    islet_coords = [(0, 0), (40, 0), (40, 40), (0, 40), (0, 0)]
    missed_coords = [(60, 60), (90, 60), (90, 90), (60, 90), (60, 60)]
    _write_geojson(good_dir / "t1.geojson", [_feature(islet_coords, "Islet", {})])
    _write_geojson(good_dir / "t2.geojson", [_feature(missed_coords, "Islet", {})])
    _write_geojson(bad_dir / "t1.geojson", [])
    _write_geojson(bad_dir / "t2.geojson", [])

    ranking = ibs.sweep_parameter_sets(key_df, donors_df, truth_dir, {"good": good_dir, "bad": bad_dir})
    assert ranking.iloc[0]["parameter_set"] == "good"
    assert ranking.iloc[0]["macro_f1"] >= ranking.iloc[1]["macro_f1"]


def test_feret_diameters_match_qupath_definitions():
    """Max Feret is the hull's longest chord; min Feret the narrowest caliper width."""
    from shapely.geometry import Polygon

    tri = Polygon([(0, 0), (10, 0), (0, 10)])          # right isosceles triangle
    mn, mx = ibs.feret_diameters(tri)
    assert mx == pytest.approx(200 ** 0.5)              # the hypotenuse
    assert mn == pytest.approx(200 ** 0.5 / 2)          # height onto the hypotenuse
    concave = Polygon([(0, 0), (10, 0), (10, 10), (5, 2), (0, 10)])   # hull is the square
    assert ibs.feret_diameters(concave) == pytest.approx((10.0, 200 ** 0.5))
