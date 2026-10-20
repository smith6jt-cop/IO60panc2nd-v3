"""Tests for scripts/v3/islet_cell_test.py (synthetic data; never touches the project).

Run:  <Senior>/.venv/bin/python -m pytest -q scripts/v3/tests/test_islet_cell_test.py
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from shapely.geometry import box, mapping

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import islet_cell_test as t  # noqa: E402

SCRIPT = Path(__file__).resolve().parents[1] / "islet_cell_test.py"


def _objects(path: Path, geoms: dict) -> None:
    feats = [{"type": "Feature", "id": f"o{i}", "geometry": mapping(g),
              "properties": {"objectType": "annotation", "classification": {"name": cls}}}
             for i, (cls, g) in enumerate(geoms)]
    path.write_text(json.dumps({"type": "FeatureCollection", "features": feats}))


def _cells(path: Path, n_ordinary: int = 200) -> pd.DataFrame:
    """Ordinary tissue cells with CD99 1..n at (1000+, 1000+) px, then the cells the tests are about."""
    rows = [{"Object ID": f"t{i}", "Parent": "Annotation (Tissue)", "x": 1000.0 + i, "y": 1000.0, "nuc": 30.0,
             "CD99": float(i + 1), "INS": 5.0, "GCG": 5.0, "SST": 5.0, "E-cadherin": 50.0, "Ker8/18": 400.0}
            for i in range(n_ordinary)]
    special = [
        # id,            parent,                 x,     y,  nucleus, CD99,   INS,   GCG, SST, E-cad, Ker8/18
        ("beta",         "Annotation (Tissue)",  5.0,  5.0, 30.0,  2000.0, 3000.0, 5.0, 5.0, 60.0, 80.0),    # in the islet object
        ("halo",         "Annotation (Tissue)",  8.0,  8.0, 30.0,   120.0,  190.0, 5.0, 5.0, 60.0, 1500.0),  # lit by its neighbour
        ("delta_no99",   "Annotation (Tissue)",  6.0,  6.0, 30.0,    20.0,    5.0, 5.0, 900.0, 60.0, 80.0),  # hormone, no CD99
        ("ecad_top",     "Annotation (Tissue)",  7.0,  7.0, 30.0,  1500.0,    5.0, 5.0, 800.0, 6000.0, 80.0),  # brighter channel on a real cell
        ("lymphocyte",   "Annotation (pLN)",   500.0, 500.0, 30.0,  9000.0,   5.0, 5.0, 5.0, 60.0, 80.0),    # lymph node: not reference
        ("no_nucleus",   "Annotation (Tissue)", 600.0, 600.0, None,  9000.0,   5.0, 5.0, 5.0, 60.0, 80.0),   # not a cell
        ("endothelium",  "Annotation (Tissue)", 305.0, 305.0, 30.0,  9000.0,   5.0, 5.0, 5.0, 60.0, 80.0),   # in a vessel object
    ]
    for cid, parent, x, y, nuc, cd99, ins, gcg, sst, ecad, ker in special:
        rows.append({"Object ID": cid, "Parent": parent, "x": x, "y": y, "nuc": nuc, "CD99": cd99, "INS": ins, "GCG": gcg,
                     "SST": sst, "E-cadherin": ecad, "Ker8/18": ker})
    d = pd.DataFrame(rows)
    out = pd.DataFrame({"Object ID": d["Object ID"], "Parent": d["Parent"], "Centroid X µm": d["x"], "Centroid Y µm": d["y"],
                        "Nucleus: Area µm^2": d["nuc"]})
    for m in ("CD99", "INS", "GCG", "SST", "E-cadherin", "Ker8/18"):
        out[f"Cell: {m}: Mean"] = d[m]
    pq.write_table(pa.Table.from_pandas(out, preserve_index=False), path)
    return d


def _run(tmp_path: Path, rule: str, *extra: str, check: bool = True):
    cmd = [sys.executable, str(SCRIPT), "--rule", rule, "--cells", str(tmp_path / "cells.parquet"),
           "--objects", str(tmp_path / "objects.geojson"), "--pixel-um", "1.0",
           "--out", str(tmp_path / f"{rule}.parquet"), "--summary", str(tmp_path / f"{rule}.json"), *extra]
    done = subprocess.run(cmd, capture_output=True, text=True)
    if check:
        assert done.returncode == 0, done.stderr
        table = pq.read_table(tmp_path / f"{rule}.parquet")
        frame = table.to_pandas().set_index("object_id")
        return frame, json.loads((tmp_path / f"{rule}.json").read_text()), table.schema.metadata
    return done


@pytest.fixture()
def image(tmp_path):
    _cells(tmp_path / "cells.parquet")
    _objects(tmp_path / "objects.geojson", [("Islet", box(0, 0, 10, 10)), ("Tissue", box(0, 0, 5000, 5000)), ("pLN", box(400, 400, 520, 520))])
    for name in t.STRUCTURES:
        (tmp_path / "structures" / name).mkdir(parents=True)
        _objects(tmp_path / "structures" / name / f"{name}__img.geojson",
                 [(name, box(300, 300, 310, 310))] if name == "Blood_vessel" else [])
    return tmp_path


def test_the_level_is_the_p99_of_ordinary_tissue_only(image):
    # without the structure objects the vessel's cell is ordinary tissue and pulls the level up
    _, loose, _ = _run(image, "cd99-level")
    assert loose["reference_population"] == 201
    _, s, meta = _run(image, "cd99-level", "--structures", str(image / "structures"), "--tag", "img")
    # 200 ordinary cells with CD99 1..200: not the lymph node cell, not the cell without a nucleus,
    # not the cell in the vessel object, not the four cells whose centroid is in the islet object
    assert s["reference_population"] == 200
    assert s["level"] == pytest.approx(np.percentile(np.arange(1, 201), 99))
    assert loose["level"] > s["level"]
    assert meta[b"rule"] == b"cd99-level" and float(meta[b"level"]) == pytest.approx(s["level"])
    assert s["reference_quantile"] == 99.0


def test_each_rule_on_the_cells_that_tell_them_apart(image):
    args = ("--structures", str(image / "structures"), "--tag", "img")
    top, s, _ = _run(image, "top-marker", *args)
    cd99, _, _ = _run(image, "cd99-level", *args)
    both, _, _ = _run(image, "cd99-or-hormone-top", *args)
    want = {            # cell:        top-marker, cd99-level, cd99-or-hormone-top
        "beta":        (True,  True,  True),
        "halo":        (False, False, False),    # the neighbour is out under every rule
        "delta_no99":  (True,  False, True),     # a hormone cell without CD99 needs the hormone
        "ecad_top":    (False, True,  True),     # a brighter channel on a real cell needs CD99
    }
    for cell, expected in want.items():
        assert (bool(top.loc[cell, "passes"]), bool(cd99.loc[cell, "passes"]), bool(both.loc[cell, "passes"])) == expected, cell
    assert top.loc["halo", "why"] == "highest marker Ker8_18"                 # names are the registry's, not the channel's
    assert top.loc["ecad_top", "top"] == "E_cadherin"
    assert both.loc["halo", "why"] == "CD99 below the level; highest marker Ker8_18"
    assert both.loc["delta_no99", "why"] == "highest marker SST" and both.loc["beta", "why"] == "CD99 at or above the level"
    assert s["markers_not_in_this_panel"] and "Ker8_18" not in s["markers_not_in_this_panel"]


def test_a_cell_without_the_measurement_cannot_pass_on_it(image):
    cells = pq.read_table(image / "cells.parquet").to_pandas()
    cells.loc[cells["Object ID"] == "beta", "Cell: CD99: Mean"] = np.nan
    cells.loc[cells["Object ID"] == "ecad_top", "Cell: CD99: Mean"] = np.nan
    pq.write_table(pa.Table.from_pandas(cells, preserve_index=False), image / "cells.parquet")
    both, _, _ = _run(image, "cd99-or-hormone-top")
    assert bool(both.loc["beta", "passes"]) and both.loc["beta", "why"] == "highest marker INS"   # the hormone still speaks
    assert not bool(both.loc["ecad_top", "passes"])
    assert both.loc["ecad_top", "why"] == "CD99 not measured; highest marker E_cadherin"


def test_rules_that_need_cd99_refuse_a_panel_without_it(image):
    cells = pq.read_table(image / "cells.parquet").to_pandas().drop(columns=["Cell: CD99: Mean"])
    pq.write_table(pa.Table.from_pandas(cells, preserve_index=False), image / "cells.parquet")
    done = _run(image, "cd99-or-hormone-top", check=False)
    assert done.returncode != 0 and "cannot be asked" in done.stderr
    top, _, _ = _run(image, "top-marker")                                      # the argmax needs no level
    assert bool(top.loc["beta", "passes"])


def test_a_level_at_the_floor_is_refused():
    with pytest.raises(SystemExit, match="at the floor"):
        t.reference_level(np.zeros(500), 99.0)
    with pytest.raises(SystemExit, match="no measurement"):
        t.reference_level(np.full(10, np.nan), 99.0)
    assert t.reference_level(np.arange(1.0, 201.0), 99.0) == pytest.approx(198.01)


def test_redsea_values_replace_the_raw_ones_and_fit_their_own_level(image):
    d = pq.read_table(image / "cells.parquet").to_pandas()
    red = pd.DataFrame({"object_id": d["Object ID"], "Cell__CD99": d["Cell: CD99: Mean"] * 0.5, "Cell__INS": d["Cell: INS: Mean"],
                        "Cell__GCG": d["Cell: GCG: Mean"], "Cell__SST": d["Cell: SST: Mean"], "Cell__Ker8_18": d["Cell: Ker8/18: Mean"],
                        "Cell__E_cadherin": d["Cell: E-cadherin: Mean"], "Cell__CD3e": np.nan})
    red = red[red["object_id"] != "t0"]                                        # a cell the source workflow did not keep
    pq.write_table(pa.Table.from_pandas(red, preserve_index=False), image / "redsea.parquet")
    both, s, meta = _run(image, "cd99-or-hormone-top", "--values", "redsea", "--redsea", str(image / "redsea.parquet"),
                         "--structures", str(image / "structures"), "--tag", "img")
    assert meta[b"values"] == b"redsea"
    assert s["level"] == pytest.approx(np.percentile(np.arange(2, 201) * 0.5, 99))     # the level is fitted on the same values
    assert s["cells_without_values"] == 1 and not bool(both.loc["t0", "passes"])
    assert "CD3e" in s["markers_not_in_this_panel"]                            # a column with no value is not a marker here
    assert bool(both.loc["beta", "passes"]) and not bool(both.loc["halo", "passes"])
    done = _run(image, "cd99-level", "--values", "redsea", check=False)
    assert done.returncode != 0 and "go together" in done.stderr
