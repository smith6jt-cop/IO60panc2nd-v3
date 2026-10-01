"""Tests for scripts/v3/islet_classifier_fit.py (synthetic data; never touches the project).

Run:  <Senior>/.venv/bin/python -m pytest -q scripts/v3/tests/test_islet_classifier_fit.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import islet_classifier_fit as f  # noqa: E402

SCRIPTS = Path(__file__).resolve().parents[2]
PLAIN = {"seed_px": 6, "close_px": 0, "open_px": 0, "margin_px": 0}       # no shaping: the test is about something else


# ----------------------------------------------------------------------------------------------- the classifier's operations

def test_the_disks_are_the_ones_qupaths_sum_filter_uses():
    # printed by QuPath 0.7 itself for op.filters.sum on a single bright pixel (2026-10-01)
    assert f.disk(1).tolist() == [[0, 1, 0], [1, 1, 1], [0, 1, 0]]
    assert f.disk(2).tolist() == [[0, 0, 1, 0, 0], [0, 1, 1, 1, 0], [1] * 5, [0, 1, 1, 1, 0], [0, 0, 1, 0, 0]]
    assert f.disk(3).tolist() == [[0, 0, 0, 1, 0, 0, 0], [0, 1, 1, 1, 1, 1, 0], [0, 1, 1, 1, 1, 1, 0], [1] * 7,
                                  [0, 1, 1, 1, 1, 1, 0], [0, 1, 1, 1, 1, 1, 0], [0, 0, 0, 1, 0, 0, 0]]
    assert [int(f.disk(r).sum()) for r in (1, 2, 3, 4)] == [5, 13, 29, 49]


def _blob(shape, cy, cx, r, value):
    yy, xx = np.ogrid[:shape[0], :shape[1]]
    return np.where((yy - cy) ** 2 + (xx - cx) ** 2 <= r * r, np.float32(value), np.float32(0))


def test_a_hormone_makes_an_islet_and_cd99_only_widens_one():
    shape = (80, 240)
    zero = np.zeros(shape, np.float32)
    planes = {"CD99": _blob(shape, 40, 40, 8, 1000)             # CD99 alone, far from any hormone
                      + _blob(shape, 40, 150, 8, 1000),         # CD99 beside the INS cell at x = 135
              "INS": _blob(shape, 40, 135, 6, 300) + _blob(shape, 40, 210, 6, 300), "GCG": zero, "SST": zero}
    thr = {"CD99": 500.0, "INS": 150.0, "GCG": 100.0, "SST": 100.0}
    m = f.classifier_mask(planes, thr, **PLAIN)
    assert m[40, 210]                                           # a hormone alone is an islet
    assert not m[40, 40]                                        # CD99 alone is not
    assert m[40, 135] and m[40, 146]                            # the INS cell, and the CD99 beside it within 6 px
    assert not m[40, 157]                                       # the same CD99 blob, farther than 6 px from the hormone
    assert f.classifier_mask(planes, thr, **{**PLAIN, "seed_px": 20})[40, 157]
    # each channel is judged against its OWN threshold
    assert not f.classifier_mask(planes, {**thr, "INS": 500.0}, **PLAIN).any()
    # a channel with no threshold is not used; without CD99 the islet is the hormone alone
    assert not f.classifier_mask(planes, {**thr, "CD99": None}, **PLAIN)[40, 146]
    with pytest.raises(SystemExit, match="no hormone"):
        f.classifier_mask(planes, {"CD99": 500.0, "INS": None, "GCG": None, "SST": None}, **PLAIN)


def test_closing_covers_the_nucleus_and_the_margin_is_the_same_for_dim_and_bright():
    shape = (80, 200)
    ring = _blob(shape, 40, 50, 7, 1000) - _blob(shape, 40, 50, 2, 1000)              # cytoplasm around a dark nucleus
    sm = {"INS": ring + _blob(shape, 40, 110, 6, 400) + _blob(shape, 40, 160, 6, 40000)}   # no smoothing: the geometry is the test
    thr = {"INS": 200.0}
    assert not f.classifier_mask(sm, thr, **PLAIN)[40, 50]                            # the centroid sits in the nucleus
    assert f.classifier_mask(sm, thr, **{**PLAIN, "close_px": 3})[40, 50]
    wide_nucleus = {"INS": _blob(shape, 40, 50, 9, 1000) - _blob(shape, 40, 50, 3, 1000)}
    assert not f.classifier_mask(wide_nucleus, thr, **{**PLAIN, "close_px": 3})[40, 50]   # a hole the disk itself fits in stays open
    assert f.classifier_mask(wide_nucleus, thr, **{**PLAIN, "close_px": 4})[40, 50]
    tight = f.classifier_mask(sm, thr, **{**PLAIN, "close_px": 3})
    wide = f.classifier_mask(sm, thr, **{**PLAIN, "close_px": 3, "margin_px": 1})
    dim = lambda m: int(m[:, 95:125].sum())       # noqa: E731
    bright = lambda m: int(m[:, 145:175].sum())   # noqa: E731
    assert dim(tight) == bright(tight)                                                # a hundred times brighter, the same outline
    assert dim(wide) - dim(tight) == bright(wide) - bright(tight) > 0


def test_opening_removes_a_speck_but_not_a_cell():
    shape = (60, 120)
    img = {"INS": _blob(shape, 30, 30, 1, 1000) + _blob(shape, 30, 80, 4, 1000)}
    kept = f.classifier_mask(img, {"INS": 500.0}, **{**PLAIN, "open_px": 2, "margin_px": 1})
    assert not kept[:, :60].any() and kept[30, 80]
    assert f.classifier_mask(img, {"INS": 500.0}, **{**PLAIN, "margin_px": 1})[30, 30]


def _types(ops):
    return [o["type"] for o in ops]


def test_the_classifier_file_says_what_the_mask_does():
    thr = {"CD99": 800.0, "INS": 1000.0, "GCG": None, "SST": 500.0}
    j = f.classifier_json(thr, seed_px=6, close_px=3, open_px=2, margin_px=1)
    assert j["pixel_classifier_type"] == "OpenCVPixelClassifier"
    assert j["metadata"]["inputResolution"]["pixelWidth"] == {"value": 1.0156, "unit": "µm"}
    assert {k: v["name"] for k, v in j["metadata"]["classificationLabels"].items()} == {"0": "Ignore*", "1": "Islet"}
    assert [t["channelName"] for t in j["op"]["colorTransforms"]] == ["CD99", "INS", "SST"]      # GCG is not used, not zero
    ops = j["op"]["op"]["ops"]
    assert _types(ops[:5]) == ["op.core.multiply", "op.filters.gaussian", "op.core.split-merge", "op.channels.sum", "op.threshold.constant"]
    assert ops[0]["values"] == pytest.approx([1 / 800.0, 1 / 1000.0, 1 / 500.0])
    assert ops[1]["sigmaX"] == pytest.approx(1.0 / 1.0156) == ops[1]["sigmaY"]
    hormone, anything = (branch["ops"] for branch in ops[2]["ops"])
    # the hormone branch: channels 1.. (not CD99), any above its threshold, grown by the seed distance
    assert _types(hormone) == ["op.channels.extract-channels", "op.channels.maximum", "op.threshold.constant", "op.filters.sum", "op.threshold.constant"]
    assert hormone[0]["channels"] == [1, 2] and hormone[2]["thresholds"] == [1.0] and hormone[3]["radius"] == 6 and hormone[4]["thresholds"] == [0.5]
    assert _types(anything) == ["op.channels.maximum", "op.threshold.constant"] and anything[1]["thresholds"] == [1.0]
    assert ops[4]["thresholds"] == [1.5]                                                         # both branches on
    sum_then = lambda radius, level: [("op.filters.sum", radius), ("op.threshold.constant", level)]  # noqa: E731
    shaping = (sum_then(3, 0.5) + sum_then(3, 28.5)        # closing: any of the 29 on, then all 29 on
               + sum_then(2, 12.5) + sum_then(2, 0.5)      # opening: all 13 on, then any on
               + sum_then(1, 0.5))                         # margin
    assert _types(ops[5:]) == [name for name, _ in shaping] + ["op.core.convert"]
    assert [o.get("radius", o.get("thresholds", [None])[0]) for o in ops[5:-1]] == [value for _, value in shaping]

    def every(o):
        yield o
        for child in o.get("ops", []):
            yield from every(child)
    # QuPath's closing / minimum / maximum filters share a kernel between threads (see the module docstring): never written
    unsafe = ("op.filters.morph-close", "op.filters.morph-open", "op.filters.minimum", "op.filters.maximum")
    assert not any(o["type"] in unsafe for top in ops for o in every(top))

    # without CD99 there is nothing to seed: one hormone, no branches, no maximum over channels
    one = f.classifier_json({"CD99": None, "INS": 800.0, "GCG": None, "SST": None}, seed_px=6, close_px=0, open_px=0, margin_px=0)
    assert [t["channelName"] for t in one["op"]["colorTransforms"]] == ["INS"]
    assert _types(one["op"]["op"]["ops"]) == ["op.core.multiply", "op.filters.gaussian", "op.threshold.constant", "op.core.convert"]
    with pytest.raises(SystemExit, match="no hormone"):
        f.classifier_json({"CD99": 800.0, "INS": None, "GCG": None, "SST": None})


# ----------------------------------------------------------------------------------------------- the fit

def _no_seed(direct):
    seed = direct.copy()
    seed[:, 0] = np.nan
    return seed


def test_fit_finds_each_hormones_own_threshold_and_leaves_useless_channels_out():
    rng = np.random.default_rng(1)
    n_other, n_a, n_b = 4000, 300, 300
    n = n_other + n_a + n_b
    cd99 = rng.uniform(1, 100, n)                                                                              # says nothing here
    ins = np.concatenate([rng.uniform(1, 10, n_other), rng.uniform(100, 200, n_a), rng.uniform(1, 10, n_b)])    # marks group a
    gcg = np.concatenate([rng.uniform(5, 50, n_other), rng.uniform(5, 50, n_a), rng.uniform(900, 1500, n_b)])   # marks group b
    sst = rng.uniform(1, 100, n)                                                                               # says nothing
    direct = np.column_stack([cd99, ins, gcg, sst])
    endocrine = np.arange(n) >= n_other
    t, loss = f.fit_thresholds(direct, _no_seed(direct), endocrine, ~endocrine, miss_weight=2.0)
    assert loss == 0
    assert 10 < t[1] < 100 and 50 < t[2] < 900                 # in the gap, not at its edge
    assert t[1] == pytest.approx(np.sqrt(10 * 100), rel=0.35)  # the middle of the run of equal costs
    assert np.isinf(t[0]) and np.isinf(t[3])


def test_cd99_is_fitted_only_where_a_hormone_seeds_it():
    rng = np.random.default_rng(4)
    n_other, n_far, n_ins, n_beside = 3000, 400, 300, 300
    n = n_other + n_far + n_ins + n_beside
    grp = np.repeat(["other", "cd99_far", "ins", "beside"], [n_other, n_far, n_ins, n_beside])
    cd99 = np.where(np.isin(grp, ["cd99_far", "beside"]), rng.uniform(800, 1200, n), rng.uniform(1, 50, n))
    ins = np.where(grp == "ins", rng.uniform(500, 900, n), rng.uniform(1, 10, n))
    ins_near = np.where(np.isin(grp, ["ins", "beside"]), rng.uniform(500, 900, n), rng.uniform(1, 10, n))      # an INS pixel within reach
    low = rng.uniform(1, 10, n)
    direct = np.column_stack([cd99, ins, low, low])
    seed = np.column_stack([np.full(n, np.nan), ins_near, low, low])
    endocrine = np.isin(grp, ["ins", "beside"])                # the CD99 cells beside an INS cell are endocrine ...
    other = ~endocrine                                         # ... the equally bright CD99 cells far from any hormone are not
    t, loss = f.fit_thresholds(direct, seed, endocrine, other, miss_weight=2.0)
    assert loss == 0 and 50 < t[0] < 800 and 10 < t[1] < 500
    c = f.covered(direct, seed, t)
    assert c[grp == "beside"].all() and c[grp == "ins"].all() and not c[grp == "cd99_far"].any()


def test_the_weight_decides_which_way_an_overlap_goes():
    rng = np.random.default_rng(2)
    n = 5500
    ins = np.concatenate([rng.normal(100, 20, 5000), rng.normal(160, 20, 500)])
    direct = np.column_stack([np.ones(n), ins, np.ones(n), np.ones(n)])
    e = np.arange(n) >= 5000
    low, _ = f.fit_thresholds(direct, _no_seed(direct), e, ~e, miss_weight=8.0)
    high, _ = f.fit_thresholds(direct, _no_seed(direct), e, ~e, miss_weight=1.0)
    assert low[1] < high[1]
    assert (ins[e] <= low[1]).sum() < (ins[e] <= high[1]).sum()


def test_a_threshold_at_the_floor_is_refused():
    rng = np.random.default_rng(3)
    n = 6000
    direct = np.column_stack([np.ones(n), rng.uniform(1, 100, n), np.ones(n), np.ones(n)])     # endocrine and other read the same
    e = np.arange(n) >= 3000
    with pytest.raises(SystemExit, match="floor"):
        f.fit_thresholds(direct, _no_seed(direct), e, ~e, miss_weight=50.0)


def test_cd99_alone_cannot_set_a_classifier():
    rng = np.random.default_rng(5)
    n = 4000
    e = np.arange(n) >= 3500
    cd99 = np.where(e, rng.uniform(800, 1200, n), rng.uniform(1, 50, n))
    low = rng.uniform(1, 10, n)
    direct = np.column_stack([cd99, low, low, low])            # no hormone tells the two groups apart
    with pytest.raises(SystemExit, match="no hormone ends up with a threshold"):
        f.fit_thresholds(direct, _no_seed(direct), e, ~e, miss_weight=0.01)


# ----------------------------------------------------------------------------------------------- script 1, emulated

def test_script_1_removes_small_pieces_and_cuts_a_bridge():
    m = np.zeros((120, 200), dtype=bool)
    m[10:16, 10:16] = True                                     # 6 x 6 px = 37 um^2: under the 50 um^2 minimum
    m[40:60, 20:40] = True                                     # two 20 x 20 px bodies joined by a 1 px bridge
    m[40:60, 60:80] = True
    m[50, 40:60] = True
    m[90:110, 100:102] = True                                  # a strip 2 px wide: the 1 um erosion leaves nothing
    lab = f.script1_objects(m)
    assert lab[12, 12] == 0
    assert lab[50, 30] > 0 and lab[50, 70] > 0 and lab[50, 30] != lab[50, 70] and lab[50, 50] == 0
    assert not lab[90:110, 100:102].any()
    assert lab.max() == 2
    # erode 1 um, dilate 0.85 um: the body keeps its place and loses only its rim
    body = lab[40:60, 20:40] > 0
    assert body[5:15, 5:15].all() and 0.75 < body.mean() <= 1.0


def test_script_1_erodes_by_the_objects_area():
    m = np.zeros((200, 200), dtype=bool)
    m[20:80, 20:80] = True                                     # 60 x 60 px = 3,713 um^2: the 2 um class
    m[50:53, 80:95] = True                                     # a spur 3 px wide survives 1 um but not 2 um
    assert not f.script1_objects(m)[51, 90]
    small = np.zeros((200, 200), dtype=bool)
    small[20:50, 20:50] = True                                 # 30 x 30 px = 928 um^2: the 1 um class
    small[30:33, 50:65] = True
    assert f.script1_objects(small)[31, 60] > 0


# ----------------------------------------------------------------------------------------------- the cells

def _cells():
    """An islet of confirmed cells at (0..56, 0..56) um, ordinary tissue far away, and the cells the rules are about."""
    rows = []

    def add(name, x, y, **m):
        rows.append({"id": name, "x": float(x), "y": float(y), "INS": 5.0, "GCG": 5.0, "SST": 5.0, "CD99": 20.0, "Ker": 300.0,
                     "CD31": 5.0, **m})
    k = 0
    for i in range(8):
        for j in range(8):
            add(f"beta{k}", 8 * i, 8 * j, INS=3000.0, CD99=1500.0 + k); k += 1
    for i in range(300):                                       # ordinary tissue: CD99 1..300, far from the islet
        add(f"t{i}", 400 + 12 * (i % 20), 400 + 12 * (i // 20), CD99=float(i + 1))
    add("capillary", 28, 28, CD31=900.0, CD99=30.0)            # enclosed by the islet's own cells
    add("ring", 70, 30, CD99=40.0)                             # the acinar neighbour
    add("hormone_no_cd99", 200, 200, SST=900.0, CD99=10.0)     # a hormone highest, no CD99
    add("cd99_beside_islet", 64, 64, CD99=2500.0, Ker=800.0)   # CD99 at the confirmed cells' level, beside the islet
    add("cd99_far_epithelial", 300, 100, CD99=2500.0, Ker=800.0)   # the same cell far from any hormone: 6482's membranes
    add("cd99_far_endothelium", 250, 300, CD99=2500.0, CD31=4000.0)
    add("vessel_beside_islet", 64, 20, CD99=2500.0, CD31=4000.0)   # CD99 on a vessel cell at the islet's edge
    add("between_levels", 66, 40, CD99=700.0)                  # above the reference level, below the confirmed cells'
    add("unmeasured", 66, 10, CD99=np.nan)                     # beside the islet, CD99 not measured: cannot be judged
    return pd.DataFrame(rows).set_index("id")


def _label(d, eligible=None, in_structure=None):
    eligible = np.ones(len(d), dtype=bool) if eligible is None else eligible
    in_structure = np.zeros(len(d), dtype=bool) if in_structure is None else in_structure
    lab = f.label_cells(d[["INS", "GCG", "SST", "CD99", "Ker", "CD31"]], d["x"].to_numpy(), d["y"].to_numpy(), eligible, in_structure,
                        other_lineage_markers={"CD31"}, quantile=99.0)
    return lab, {name: i for i, name in enumerate(d.index)}


def test_cells_are_endocrine_other_or_left_out():
    d = _cells()
    lab, at = _label(d)
    ordinary = np.array([n.startswith("t") for n in d.index])
    # the level is the ordinary tissue's alone: nothing within 30 um of a hormone-highest cell is in it
    assert lab["reference"][ordinary].all() and not lab["reference"][at["ring"]] and not lab["reference"][at["beta0"]]
    assert lab["cd99_reference_level"] == pytest.approx(np.percentile(d["CD99"][lab["reference"]].dropna(), 99.0))
    assert 250 < lab["cd99_reference_level"] < 700 < lab["cd99_confirmed_level"]
    assert lab["confirmed"][[at[f"beta{k}"] for k in range(64)]].all() and int(lab["confirmed"].sum()) == 64
    # CD99 at the confirmed cells' level makes a cell endocrine only beside an islet
    assert lab["by_cd99"][at["cd99_beside_islet"]] and lab["endocrine"][at["cd99_beside_islet"]]
    for name in ("cd99_far_epithelial", "cd99_far_endothelium"):       # ordinary tissue, whatever its CD99
        assert lab["other"][at[name]] and lab["reference"][at[name]] and not lab["endocrine"][at[name]], name
    for name in ("hormone_no_cd99", "vessel_beside_islet", "between_levels", "unmeasured"):
        assert lab["left_out"][at[name]] and not lab["endocrine"][at[name]] and not lab["other"][at[name]], name
    assert lab["other"][at["ring"]] and not lab["enclosed"][at["ring"]]
    assert lab["enclosed"][at["capillary"]] and not lab["other"][at["capillary"]]      # the islet's own: not counted against
    assert lab["other"][ordinary].all()
    # every cell is in exactly one of the groups
    groups = lab["endocrine"].astype(int) + lab["other"] + lab["enclosed"] + lab["left_out"]
    assert (groups == 1).all()


def test_a_cell_without_a_nucleus_or_outside_tissue_is_in_no_group():
    d = _cells()
    eligible = np.ones(len(d), dtype=bool)
    _, at = _label(d)
    eligible[at["beta0"]] = False
    eligible[at["ring"]] = False
    lab, _ = _label(d, eligible)
    for name in ("beta0", "ring"):
        assert not any(lab[k][at[name]] for k in ("endocrine", "other", "enclosed", "left_out", "reference", "confirmed"))


def test_structure_cells_leave_the_reference_population():
    d = _cells()
    lab0, at = _label(d)
    in_structure = np.zeros(len(d), dtype=bool)
    in_structure[[at[f"t{i}"] for i in range(250, 300)]] = True          # the CD99-brightest ordinary cells sit in a vessel object
    in_structure[[at["cd99_far_epithelial"], at["cd99_far_endothelium"]]] = True
    lab, _ = _label(d, in_structure=in_structure)
    assert lab["cd99_reference_level"] < lab0["cd99_reference_level"]


def test_too_few_confirmed_cells_is_a_refusal_not_a_classifier():
    d = _cells()
    d.loc[[f"beta{k}" for k in range(30, 64)], "INS"] = 5.0             # 30 confirmed cells left
    with pytest.raises(SystemExit, match="cannot set a classifier"):
        _label(d)


# ----------------------------------------------------------------------------------------------- the one changed line

def test_the_per_image_script_is_the_maintainers_script_with_one_line_changed():
    theirs = (SCRIPTS / "1_panc_objects.groovy").read_text().splitlines()
    ours = (SCRIPTS / "v3" / "qupath" / "1_panc_objects_perimage.groovy").read_text().splitlines()
    start = len(ours) - len(theirs)
    assert start > 0 and ours[0] == "/*" and ours[start - 1] == " */"     # only a header comment before their first line
    assert all(line.startswith(("/*", " *")) for line in ours[:start])
    body = ours[start:]
    changed = [(a, b) for a, b in zip(theirs, body) if a != b]
    assert len(changed) == 1
    assert changed[0][0].strip() == 'def isletClassifier = loadPixelClassifier("Islets2")'
    assert changed[0][1].strip().startswith('def isletClassifier = loadPixelClassifier("Islets_" + getProjectEntry().getImageName()')
