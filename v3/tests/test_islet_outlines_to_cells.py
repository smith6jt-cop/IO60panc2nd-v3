"""Tests for the islet step's second half (maintainer, 2026-10-02): outlines erode to the cell outlines; an outline no
whole cell fits in is a Hormone Fragment. Synthetic shapes; never touches the project.

Run:  <Senior>/.venv/bin/python -m pytest -q scripts/v3/tests/test_islet_outlines_to_cells.py
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np
import pytest
import shapely
from shapely.geometry import box

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import islet_classifier_fit as fit  # noqa: E402
import islet_outlines_to_cells as oc  # noqa: E402

V3 = Path(__file__).resolve().parents[1]
PX = oc.FULL_RES_UM


def grid(nx: int, ny: int, size: float = 20.0, gap: float = 1.0, x0: float = 0.0, y0: float = 0.0) -> list:
    """Cells the way InstanSeg leaves them: squares of `size` px with a `gap` px strip of nothing between neighbours."""
    return [box(x0 + i * (size + gap), y0 + j * (size + gap), x0 + i * (size + gap) + size, y0 + j * (size + gap) + size)
            for j in range(ny) for i in range(nx)]


def run(objs, cells, nucleus=None, **kw):
    nucleus = np.ones(len(cells), dtype=bool) if nucleus is None else np.asarray(nucleus, dtype=bool)
    return oc.erode_to_cells(np.array(objs, dtype=object), np.array(cells, dtype=object), nucleus, **kw)


# ----------------------------------------------------------------------------------------------- the erosion

def test_the_outline_ends_on_the_cells_that_fit_and_the_cut_cells_are_out():
    cells = grid(5, 5)                                        # the middle 3 x 3 lies whole in the outline
    outline = box(15, 15, 89, 89)                             # cuts through the 16 cells around them
    r = run([outline], cells)
    islet, = r["islets"]
    inner = [i for i, c in enumerate(cells) if outline.covers(c)]
    assert len(inner) == 9 and sorted(islet["cells"].tolist()) == inner
    assert r["fragments"] == []
    g = islet["geometry"]
    for i, c in enumerate(cells):
        if i in inner:
            assert g.covers(c), "a cell that fits is kept whole"
        else:
            assert shapely.intersection(g, c).area < 1e-9, "a cell the outline cut is taken out of it"
    # the edge now runs along the outlines of the cut cells: the nine cells (21..83), the gaps between them and the
    # 1-px gap up to the cut cells on every side; the strips between two cut cells, out to the old outline, are gone
    # (to within a tenth of a pixel: where a removed strip met the islet, a disc-shaped opening leaves a nub that small)
    assert g.bounds == pytest.approx((20.0, 20.0, 84.0, 84.0), abs=0.15)
    assert g.area == pytest.approx(64.0 * 64.0, abs=3.0)        # and it rounds the four corners of the gap ring
    assert r["islet_of"][inner].tolist() == [0] * 9 and (np.delete(r["islet_of"], inner) == -1).all()


def test_one_whole_cell_takes_the_cells_own_outline():
    cells = grid(3, 3)
    outline = box(19, 19, 43, 43)                             # holds the centre cell (21..41) and cuts the eight around it
    r = run([outline], cells)
    islet, = r["islets"]
    assert islet["cells"].tolist() == [4]
    assert islet["geometry"].equals(cells[4]), '"if one whole cell fits in islet outline, use cell segmentation"'


def test_an_outline_no_whole_cell_fits_in_is_a_hormone_fragment():
    cells = grid(3, 3)
    blob = box(24, 24, 36, 36)                                # inside the centre cell: smaller than a cell
    straddle = box(38, 5, 46, 15)                             # across two cells, holding neither
    r = run([blob, straddle], cells)
    assert r["islets"] == [] and r["fragments"] == [0, 1]
    assert (r["home"] == -1).all() and (r["islet_of"] == -1).all()


def test_a_cell_fits_when_it_reaches_no_farther_than_one_classifier_pixel_beyond_the_outline():
    cell = box(0, 0, 20, 20)
    one_px = oc.FIT_TOLERANCE_UM / PX                          # 2 full-resolution pixels
    just = box(0, 0, 20 - (one_px - 0.1), 20)                  # the cell reaches 1.9 px beyond it
    over = box(0, 0, 20 - (one_px + 0.1), 20)                  # 2.1 px
    assert len(run([just], [cell])["islets"]) == 1
    assert run([over], [cell])["fragments"] == [0]
    assert run([just], [cell], tol_um=0.0)["fragments"] == [0], "tolerance 0 is the strict reading"
    # most of the cell must lie inside as well: a sliver of an outline beside a small cell does not hold it
    assert run([box(-1.9, 0, 0.4, 2)], [box(0, 0, 2, 2)])["fragments"] == [0]


def test_a_detection_without_a_nucleus_is_not_a_whole_cell():
    cells = grid(2, 1)
    outline = box(-1, -1, 42, 21)
    with_nuclei = run([outline], cells)
    assert len(with_nuclei["islets"]) == 1 and len(with_nuclei["islets"][0]["cells"]) == 2
    one = run([outline], cells, nucleus=[True, False])
    islet, = one["islets"]
    assert islet["cells"].tolist() == [0] and islet["geometry"].equals(cells[0])
    none = run([outline], cells, nucleus=[False, False])
    assert none["islets"] == [] and none["fragments"] == [0], "an outline holding only detections without a nucleus"


def test_an_outline_splits_where_only_cut_cells_join_its_parts():
    cells = grid(5, 1)                                         # a row: 0 1 | 2 | 3 4
    outline = shapely.union_all([box(-1, -1, 42, 21), box(62, -1, 105, 21), box(40, 8, 64, 12)])   # a thin neck across cell 2
    r = run([outline], cells)
    assert sorted(sorted(i["cells"].tolist()) for i in r["islets"]) == [[0, 1], [3, 4]]
    assert r["home"][2] == -1 and all(shapely.intersection(i["geometry"], cells[2]).area < 1e-9 for i in r["islets"])


def test_a_detection_fits_one_outline_only():
    cell = box(0, 0, 20, 20)
    a, b = box(-1, -1, 19, 21), box(0.5, -1, 22, 21)           # both hold it within the tolerance; b covers more of it
    r = run([a, b], [cell])
    assert r["home"].tolist() == [1]
    assert r["fragments"] == [0] and len(r["islets"]) == 1 and r["islets"][0]["source"] == 1


def test_no_whisker_is_left_between_two_cut_cells():
    cells = grid(4, 4)
    outline = box(15, 15, 68, 68)                             # the inner 2 x 2 fits; twelve cells are cut
    islet, = run([outline], cells)["islets"]
    g = islet["geometry"]
    assert g.geom_type == "Polygon" and len(g.interiors) == 0
    assert g.area == pytest.approx(43.0 * 43.0, abs=3.0)       # the four cells (21..62) and the 1-px gaps around them; no strips outward
    assert run([outline], cells, whisker_um=0.0)["islets"][0]["geometry"].area > g.area + 20, "without the opening the strips stay"


# ----------------------------------------------------------------------------------------------- the marker classifier

THRESHOLDS = {"CD99": 800.0, "INS": 500.0, "GCG": None, "SST": 250.0}


def test_the_marker_classifier_is_the_hormone_branch_and_nothing_else():
    j = fit.marker_json(THRESHOLDS)
    assert [c["channelName"] for c in j["op"]["colorTransforms"]] == ["INS", "SST"], "no CD99, and no hormone that carries no threshold"
    ops = j["op"]["op"]["ops"]
    assert [o["type"] for o in ops] == ["op.core.multiply", "op.filters.gaussian", "op.channels.maximum", "op.threshold.constant", "op.core.convert"]
    assert ops[0]["values"] == [1 / 500.0, 1 / 250.0] and ops[3]["thresholds"] == [1.0]
    assert ops[1]["sigmaX"] == pytest.approx(fit.SIGMA_UM / fit.CLASSIFIER_UM)
    assert j["metadata"]["classificationLabels"]["1"]["name"] == oc.CLASS_FRAGMENT
    assert j["metadata"]["inputResolution"]["pixelWidth"]["value"] == fit.CLASSIFIER_UM
    # the filters QuPath 0.7 cannot run in parallel (qupath/qupath#2202) are not in it
    assert not any(o["type"] in ("op.filters.maximum", "op.filters.minimum", "op.filters.morph-open", "op.filters.morph-close") for o in ops)
    one = fit.marker_json({"CD99": 1.0, "INS": None, "GCG": 40.0, "SST": None})
    assert [o["type"] for o in one["op"]["op"]["ops"]] == ["op.core.multiply", "op.filters.gaussian", "op.threshold.constant", "op.core.convert"]
    with pytest.raises(SystemExit, match="no hormone"):
        fit.marker_json({"CD99": 1.0, "INS": None, "GCG": None, "SST": None})


def test_the_marker_mask_is_what_the_classifier_file_says():
    smoothed = {"CD99": np.full((4, 4), 5000.0), "INS": np.zeros((4, 4)), "GCG": np.full((4, 4), 9999.0), "SST": np.zeros((4, 4))}
    smoothed["INS"][0, 0] = 501.0; smoothed["INS"][0, 1] = 500.0; smoothed["SST"][3, 3] = 251.0
    mask = fit.marker_mask(smoothed, THRESHOLDS)               # GCG carries no threshold and CD99 is not a hormone
    assert mask.sum() == 2 and mask[0, 0] and mask[3, 3] and not mask[0, 1], "strictly above, as op.threshold.constant"
    assert fit.marker_classifier_path(Path("x/Islets_6539_Scan1.json")).name == "IsletHormone_6539_Scan1.json"
    with pytest.raises(SystemExit):
        fit.marker_classifier_path(Path("x/Other.json"))


# ----------------------------------------------------------------------------------------------- the step in the workflow

def steps() -> dict:
    lines = [ln for ln in (V3 / "steps.tsv").read_text().splitlines() if ln and not ln.startswith("#")]
    return {r["step"]: r for r in csv.DictReader(lines, delimiter="\t")}


def test_the_step_runs_after_the_islets_and_the_cells_and_before_everything_that_reads_islets():
    s = steps()
    assert s["1b"]["script"] == "scripts/v3/qupath/1b_islet_outlines_to_cells.groovy" and "I" in s["1b"]["passes"].split(",")
    assert s["1b"]["after"].split(",") == ["1", "0?", "1_outlines?"] and s["1b"]["save"] == "footer-saves" and s["1b"]["postcondition"] == "v3_result"
    # script 1's outlines are exported as it drew them before 1b changes them
    assert s["1_outlines"]["args"] == "islets|{OUT}/script1_outlines/{TAG}.geojson" and s["1_outlines"]["save"] == "never"
    assert list(s).index("1") < list(s).index("1_outlines") < list(s).index("1b") < list(s).index("4")
    for later in ("3", "8", "J_annotations", "J_islets"):      # pass J cannot run on outlines that were not eroded
        assert "1b" in s[later]["after"].split(","), later
    for later in ("4", "R_export"):                            # the Islets2 comparison never runs 1b
        assert "1b?" in s[later]["after"].split(","), later


def test_the_groovy_step_carries_the_same_numbers_and_names_as_the_reference():
    g = (V3 / "qupath" / "1b_islet_outlines_to_cells.groovy").read_text()
    assert f"double FIT_TOLERANCE_UM = {oc.FIT_TOLERANCE_UM}" in g and f"double WHISKER_UM = {oc.WHISKER_UM}" in g
    assert f"String CLASS_FRAGMENT = '{oc.CLASS_FRAGMENT}'" in g and "String MARKER_CLASSIFIER_PREFIX = 'IsletHormone_'" in g
    assert "JOIN_MITRE" not in g, "mitre joins collapse some outlines in JTS 1.20; the opening uses round joins"
    assert "openings_skipped" in g and "0.9 * want" in g, "an opening that takes the cells away is refused"
    # no measurement of any cell is read: the step is geometry
    assert "getMeasurementList().get" not in g and "measurements[" not in g.replace("obj.measurements[key]", "")
    # a Hormone Fragment is never an islet: script 8 parents cells to these classes only
    parenting = (V3 / "qupath" / "8_hierarchy_parenting_v3.groovy").read_text()
    order = parenting[parenting.index("List<String> PARENT_ORDER"):].splitlines()[0]
    assert oc.CLASS_FRAGMENT not in order and "'Islet'" in order


def test_script_1_itself_is_still_the_maintainers_with_one_line_changed():
    mine = (V3 / "qupath" / "1_panc_objects_perimage.groovy").read_text().splitlines()
    theirs = (V3.parent / "1_panc_objects.groovy").read_text().splitlines()
    body = mine[len(mine) - len(theirs):]
    changed = [(a, b) for a, b in zip(theirs, body) if a != b]
    assert len(changed) == 1 and 'loadPixelClassifier("Islets2")' in changed[0][0]


# ----------------------------------------------------------------------------------------------- the check's own reading

def test_the_check_takes_each_fragment_back_to_the_outline_it_lay_in(tmp_path):
    """A fragment shares an edge with a large neighbouring outline: only the outline it lay IN is its source (the
    first version counted every cell of the neighbour as lying under a fragment)."""
    import islet_classifier_check as chk  # noqa: PLC0415
    big, small = box(0, 0, 100, 100), box(100, 40, 110, 50)
    fragment = box(100, 42, 104, 46)                            # inside `small`, touching `big` along x = 100
    x = np.array([50.0, 105.0, 300.0]); y = np.array([50.0, 45.0, 300.0])
    lab = {"confirmed": np.array([True, False, False]), "endocrine": np.array([True, False, False]),
           "other": np.array([False, True, False]), "enclosed": np.zeros(3, bool), "left_out": np.zeros(3, bool)}
    out = chk.fragments_report([fragment], [big, small], x, y, lab, np.ones(3, bool))
    assert out["hormone_fragments"] == 1
    assert out["cells_with_centre_in_a_fragments_old_outline"] == {"all": 1, "confirmed_endocrine": 0, "endocrine_by_cd99": 0,
                                                                   "other": 1, "left_out": 0}


def test_the_export_is_read_by_class_and_only_islets_count_as_islets(tmp_path):
    import islet_classifier_check as chk  # noqa: PLC0415
    from shapely.geometry import mapping  # noqa: PLC0415
    feats = [{"type": "Feature", "geometry": mapping(box(0, 0, 10, 10)),
              "properties": {"classification": {"name": "Islet"}, "measurements": {"Whole cells": 3.0}}},
             {"type": "Feature", "geometry": mapping(box(20, 0, 23, 3)), "properties": {"classification": {"name": oc.CLASS_FRAGMENT}}},
             {"type": "Feature", "geometry": mapping(box(0, 0, 500, 500)), "properties": {"classification": {"name": "Tissue"}}}]
    f = tmp_path / "x.geojson"
    f.write_text(json.dumps({"type": "FeatureCollection", "features": feats}))
    got = chk.read_by_class(f)
    assert sorted(got["by_class"]) == [oc.CLASS_FRAGMENT, "Islet"] and got["whole_cells"] == [3.0]
