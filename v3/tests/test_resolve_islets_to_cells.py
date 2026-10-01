"""Tests for scripts/v3/resolve_islets_to_cells.py (synthetic geometry only; never touches the project).

Run:  <Senior>/.venv/bin/python -m pytest -q scripts/v3/tests/test_resolve_islets_to_cells.py
"""
from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

import pytest
from shapely.geometry import Polygon, box, mapping
from shapely.geometry.polygon import orient

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import resolve_islets_to_cells as r  # noqa: E402

SCRIPT = Path(__file__).resolve().parents[1] / "resolve_islets_to_cells.py"


# ---------------------------------------------------------------- diameters

@pytest.mark.parametrize("sign", [1.0, -1.0])
def test_min_caliper_width_does_not_depend_on_ring_direction(sign):
    # GEOS returns the hull clockwise; a normal pointing out of the hull once made every width 0
    rect = orient(box(0, 0, 10, 4), sign=sign)
    assert r.min_caliper_width(rect) == pytest.approx(4.0)


def test_min_caliper_width_of_a_rotated_rectangle_and_a_triangle():
    c, s = math.cos(0.5), math.sin(0.5)
    rot = Polygon([(x * c - y * s, x * s + y * c) for x, y in [(0, 0), (12, 0), (12, 3), (0, 3)]])
    assert r.min_caliper_width(rot) == pytest.approx(3.0)
    tri = Polygon([(0, 0), (6, 0), (0, 8)])          # smallest altitude: onto the hypotenuse (10)
    assert r.min_caliper_width(tri) == pytest.approx(4.8)


def test_min_caliper_width_uses_the_hull_of_every_part():
    two = box(0, 0, 2, 2).union(box(8, 0, 10, 2))
    assert r.min_caliper_width(two) == pytest.approx(2.0)


def test_size_class_boundaries_are_drotars():
    assert r.size_class(10.0) == "EndocrineSingle"
    assert r.size_class(10.01) == "EndocrineCluster"
    assert r.size_class(34.99) == "EndocrineCluster"
    assert r.size_class(35.0) == "Islet"


# ---------------------------------------------------------------- joining cells

def test_cells_one_pixel_apart_become_one_region_with_the_seam_closed():
    a, b = box(0, 0, 10, 10), box(10.98, 0, 20, 10)
    plain = r.join_cells([a, b], 0.0)
    assert plain.geom_type == "MultiPolygon"
    joined = r.join_cells([a, b], 1.0)
    assert joined.geom_type == "Polygon"
    assert joined.area == pytest.approx(200.0)                  # the two cells + the 0.98 x 10 seam
    assert joined.bounds == pytest.approx((0, 0, 20, 10))       # nothing added at the outer boundary


def test_cells_two_pixels_apart_stay_apart_and_keep_their_outlines():
    a, b = box(0, 0, 10, 10), box(12, 0, 20, 10)
    joined = r.join_cells([a, b], 1.0)
    assert joined.geom_type == "MultiPolygon"
    assert joined.area == pytest.approx(a.area + b.area)


def test_a_single_cell_is_returned_unchanged():
    a = Polygon([(0, 0), (5, 0), (5, 1), (6, 1), (6, 5), (0, 5)])
    assert r.join_cells([a], 1.0).equals(a)


def test_outer_boundary_closes_enclosed_holes_only():
    ring = box(0, 0, 10, 10).difference(box(4, 4, 6, 6))
    assert r.outer_boundary(ring).area == pytest.approx(100.0)
    horseshoe = box(0, 0, 10, 10).difference(box(4, 4, 6, 10))
    assert r.outer_boundary(horseshoe).area == pytest.approx(horseshoe.area)


# ---------------------------------------------------------------- the rule, end to end

def _cell(i, geom, nucleus=True):
    f = {"type": "Feature", "id": f"cell-{i}", "geometry": mapping(geom), "properties": {"objectType": "cell"}}
    if nucleus:
        f["nucleusGeometry"] = mapping(geom.buffer(-1))
    return f


def test_membership_is_majority_of_area_with_a_nucleus(tmp_path):
    obj = box(0, 0, 30, 14)                                    # one pixel object
    cells = [
        _cell(0, box(0, 0, 10, 10)),                           # wholly inside            -> member
        _cell(1, box(10.98, 0, 20, 10)),                       # wholly inside            -> member
        _cell(2, box(24, 0, 34, 10)),                          # 60 % inside              -> member
        _cell(3, box(0, 10.98, 10, 20.98)),                    # 30 % inside              -> not
        _cell(4, box(11, 15, 20, 25)),                         # outside                  -> not
        _cell(5, box(21, 0, 23, 10), nucleus=False),           # inside, no nucleus       -> not a cell
    ]
    objects = {"type": "FeatureCollection", "features": [
        {"type": "Feature", "id": "obj-1", "geometry": mapping(obj),
         "properties": {"objectType": "annotation", "classification": {"name": "Islet"}}},
        {"type": "Feature", "id": "tissue", "geometry": mapping(box(-50, -50, 100, 100)),
         "properties": {"objectType": "annotation", "classification": {"name": "Tissue"}}},
        {"type": "Feature", "id": "obj-2", "geometry": mapping(box(60, 60, 62, 62)),      # covers no cell
         "properties": {"objectType": "annotation", "classification": {"name": "EndocrineSingle"}}}]}
    (tmp_path / "objects.geojson").write_text(json.dumps(objects))
    (tmp_path / "cells.geojson").write_text(json.dumps({"type": "FeatureCollection", "features": cells}))
    subprocess.run([sys.executable, str(SCRIPT), "--objects", str(tmp_path / "objects.geojson"),
                    "--cells", str(tmp_path / "cells.geojson"), "--out", str(tmp_path / "out.geojson"),
                    "--members", str(tmp_path / "members.tsv"), "--summary", str(tmp_path / "summary.json"),
                    "--pixel-um", "1.0"], check=True, capture_output=True)
    summary = json.loads((tmp_path / "summary.json").read_text())
    members = [line.split("\t")[0] for line in (tmp_path / "members.tsv").read_text().splitlines()[1:]]
    assert sorted(members) == ["cell-0", "cell-1", "cell-2"]
    assert summary["majority_cells_without_nucleus"] == 1
    assert summary["pixel_objects"] == 2 and summary["resolved_objects"] == 1
    assert summary["dropped_no_member_by_pixel_class"] == {"EndocrineSingle": 1}
    out = json.loads((tmp_path / "out.geojson").read_text())["features"]
    assert [f["id"] for f in out] == ["obj-1"]                  # the pixel object's id is kept
    m = out[0]["properties"]["measurements"]
    assert m["n_cells"] == 3
    # cells 0 and 1 join across the one-pixel seam; cell 2 is 4 px away and stays a separate part
    assert m["n_parts"] == 2
    assert m["Min diameter um"] == pytest.approx(10.0)
    assert m["Max diameter um"] == pytest.approx(math.hypot(34, 10), rel=1e-6)
    assert out[0]["properties"]["classification"]["name"] == "EndocrineCluster"   # (35.4 + 10) / 2 = 22.7
