"""Tests for scripts/v3/islet_export_adapter.py (synthetic data; never touches the project).

The last two tests read the adapter's output with the Senior typing lane's OWN readers, from the v3
worktree, so a format the lane cannot read fails here and not in a cohort run.

Run:  <Senior>/.venv/bin/python -m pytest -q scripts/v3/tests/test_islet_export_adapter.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from shapely.geometry import box, mapping, shape

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import islet_export_adapter as a  # noqa: E402

SENIOR = Path("/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior-v3")
D = a.DISTANCE


def _feature(qclass: str, name: str, geom, uid: str = "u", **meas) -> dict:
    return {"type": "Feature", "id": uid, "geometry": mapping(geom),
            "properties": {"objectType": "annotation", "name": name, "classification": {"name": qclass},
                           "measurements": meas}}


# ---------------------------------------------------------------- numbering and selection

def test_typing_numbers_keep_the_qupath_object_readable():
    assert a.typing_number("Islet", "Islet_12") == 12
    assert a.typing_number("EndocrineCluster", "Cluster_7") == 100_007
    assert a.typing_number("EndocrineSingle", "Single_3") == 200_003
    for qclass, bad in (("Islet", "Cluster_7"), ("EndocrineCluster", None), ("Islet", "Islet_0"),
                        ("EndocrineSingle", "Single_100000")):
        with pytest.raises(ValueError):
            a.typing_number(qclass, bad)


def test_which_objects_count_is_the_ruling_not_a_default():
    feats = [_feature("Islet", "Islet_1", box(0, 0, 10, 10)), _feature("EndocrineCluster", "Cluster_1", box(50, 0, 55, 5)),
             _feature("EndocrineSingle", "Single_1", box(90, 0, 92, 2)), _feature("Pixel object (before cells)", "pixel Islet", box(0, 0, 1, 1))]
    assert sorted(a.select_islets(feats, "all")) == [1, 100_001, 200_001]
    assert sorted(a.select_islets(feats, "islet")) == [1]
    with pytest.raises(ValueError, match="two objects"):
        a.select_islets(feats + [_feature("Islet", "Islet_1", box(0, 20, 10, 30))], "all")
    args = ["--export", "x", "--senior", "y"]
    for missing in (["--islets", "all"], ["--expansion", "none"]):       # both rulings must be stated
        with pytest.raises(SystemExit):
            sys.argv = ["islet_export_adapter.py", *args, *missing]
            a.main()


# ---------------------------------------------------------------- expansions

def test_bands_exclude_every_islet_and_split_at_equal_distance():
    cores = {1: box(0, 0, 100, 100), 2: box(120, 0, 220, 100), 3: box(1000, 1000, 1010, 1010)}
    bands = a.nearest_islet_bands(cores, width_px=40.0)
    # an islet that meets no other keeps the whole 40 px zone around it
    assert bands[3].area == pytest.approx(box(1000, 1000, 1010, 1010).buffer(40).area - 100, rel=1e-6)
    # the two neighbours are 20 px apart: their zones meet and are cut at x = 110
    for n in (1, 2):
        assert bands[n].intersection(cores[1]).area == pytest.approx(0, abs=1e-6)
        assert bands[n].intersection(cores[2]).area == pytest.approx(0, abs=1e-6)
    assert bands[1].intersection(bands[2]).area == pytest.approx(0, abs=1e-6)
    assert bands[1].bounds[2] == pytest.approx(110, abs=1.5) and bands[2].bounds[0] == pytest.approx(110, abs=1.5)
    assert bands[1].area == pytest.approx(bands[2].area, rel=0.02)          # mirror images
    assert bands[1].contains(box(101, 40, 109, 60)) and bands[2].contains(box(111, 40, 119, 60))


def test_a_cell_sized_islet_beside_a_large_one_gets_the_side_nearer_to_it():
    cores = {1: box(0, 0, 400, 400), 200_001: box(420, 190, 440, 210)}
    bands = a.nearest_islet_bands(cores, width_px=40.0)
    assert bands[200_001].contains(box(442, 195, 470, 205))                  # beyond the single: its own
    assert bands[1].contains(box(402, 0, 408, 50))                           # far along the big one: the big one's
    assert bands[1].intersection(bands[200_001]).area == pytest.approx(0, abs=1e-6)


# ---------------------------------------------------------------- the table

def _tables(tmp_path: Path, drop_context_row: bool = False):
    cells = pd.DataFrame({
        "Object ID": ["c1", "c2", "c3", "c4", "c5"], "Image": ["img - resolution #1"] * 5,
        "Classification": [None] * 5, "Centroid X µm": [10.0, 20.0, 30.0, 40.0, 50.0],
        "Centroid Y µm": [1.0, 2.0, 3.0, 4.0, 5.0], "Nucleus: Area µm^2": [30.0, None, 28.0, 31.0, 29.0],
        "Cell: Area µm^2": [60.0, 55.0, 58.0, 61.0, 59.0], "Cell: Solidity": [0.9] * 5,
        "Nucleus/Cell area ratio": [0.5, None, 0.48, 0.5, 0.49], "Cell: DAPI: Mean": [1.0] * 5})
    context = pd.DataFrame({
        "Object ID": ["c1", "c2", "c3", "c4", "c5"],
        "Parent": ["Islet_4", "Cluster_2", "Single_9", "Annotation (Tissue)", "Annotation (pLN)"],
        "Centroid X µm": [10.0, 20.0, 30.0, 40.0, 50.0], "Centroid Y µm": [1.0, 2.0, 3.0, 4.0, 5.0],
        D.format("Islet"): [-3.0, 80.0, 200.0, 15.0, float("nan")],
        D.format("EndocrineCluster"): [40.0, -2.0, 90.0, 60.0, float("nan")],
        D.format("EndocrineSingle"): [70.0, 30.0, -1.0, 12.0, float("nan")]})
    if drop_context_row:
        context = context.iloc[:4]
    pq.write_table(pa.Table.from_pandas(cells), tmp_path / "cells.parquet")
    pq.write_table(pa.Table.from_pandas(context), tmp_path / "context.parquet")
    return tmp_path / "cells.parquet", tmp_path / "context.parquet"


def _chosen(rule: str) -> dict:
    feats = [_feature("Islet", "Islet_4", box(0, 0, 10, 10)), _feature("EndocrineCluster", "Cluster_2", box(50, 0, 55, 5)),
             _feature("EndocrineSingle", "Single_9", box(90, 0, 92, 2))]
    return a.select_islets(feats, rule)


def test_table_with_every_endocrine_object_as_an_islet(tmp_path):
    cells, context = _tables(tmp_path)
    stats = a.write_tsv(cells, context, tmp_path / "d.tsv", _chosen("all"), "all", scale=0.5)
    t = pd.read_csv(tmp_path / "d.tsv", sep="\t").set_index("Object ID")
    assert list(t.columns) == a.TSV_COLUMNS[1:]
    assert t["Parent"].tolist() == ["Islet_4", "Islet_100002", "Islet_200009", "Annotation (Tissue)", "Annotation (pLN)"]
    # the distance to the nearest object that counts, rescaled; a cell with no distance stays without one
    assert t[D.format("Islet")].tolist()[:4] == [-1.5, -1.0, -0.5, 6.0]
    assert pd.isna(t.loc["c5", D.format("Islet")])
    assert t["Centroid X µm"].tolist() == [5.0, 10.0, 15.0, 20.0, 25.0]
    assert pd.isna(t.loc["c2", "Nucleus: Area µm^2"])                        # no nucleus stays unavailable
    assert stats == {"rows": 5, "parent_is_islet": 3, "parent_is_pln": 1, "distance_inside_islet": 3, "distance_missing": 1}


def test_table_with_only_large_islets_puts_the_other_cells_in_tissue(tmp_path):
    cells, context = _tables(tmp_path)
    a.write_tsv(cells, context, tmp_path / "d.tsv", _chosen("islet"), "islet", scale=1.0)
    t = pd.read_csv(tmp_path / "d.tsv", sep="\t").set_index("Object ID")
    assert t["Parent"].tolist() == ["Islet_4", "Annotation (Tissue)", "Annotation (Tissue)", "Annotation (Tissue)", "Annotation (pLN)"]
    assert t[D.format("Islet")].tolist()[:4] == [-3.0, 80.0, 200.0, 15.0]    # the Islet column alone


def test_exports_of_different_cells_are_refused(tmp_path):
    cells, context = _tables(tmp_path, drop_context_row=True)
    with pytest.raises(ValueError, match="not of the same cells"):
        a.write_tsv(cells, context, tmp_path / "d.tsv", _chosen("all"), "all", scale=1.0)
    assert not (tmp_path / "d.tsv").exists() and not (tmp_path / "d.tsv.part").exists()


# ---------------------------------------------------------------- read back by the Senior lane itself

@pytest.fixture()
def senior_typing():
    typing = SENIOR / "core/phenocycler/scripts/typing"
    if not typing.is_dir():
        pytest.skip("the v3 Senior worktree is not on this machine")
    sys.path.insert(0, str(typing))
    return typing


def test_the_lanes_own_readers_take_the_islet_file(tmp_path, monkeypatch, senior_typing):
    import build_islet_geometry as geometry
    import islet_membership as membership
    chosen = a.select_islets([_feature("Islet", "Islet_4", box(0, 0, 100, 100), n_cells=30),
                              _feature("EndocrineSingle", "Single_9", box(130, 0, 140, 10), n_cells=1)], "all")
    bands = a.nearest_islet_bands({n: o["geom"] for n, o in chosen.items()}, 40.0)
    folder = tmp_path / "typing_islets"
    folder.mkdir()
    (folder / "Islets__6539_Scan1.er.qptiff_-_resolution_1.geojson").write_text(
        json.dumps({"type": "FeatureCollection", "features": a.islet_features(chosen, bands)}))
    monkeypatch.setattr(geometry.S, "ISLETS", folder)
    core, expansion = membership.outlines("6539_Scan1.er")
    assert sorted(core) == ["200009", "4"] and sorted(expansion) == ["200009", "4"]
    assert core["4"].area == pytest.approx(10_000)
    assert sorted(geometry.islet_outlines("6539_Scan1.er")) == ["Islet_200009", "Islet_4"]   # cores only


def test_the_lanes_own_reader_takes_the_table(tmp_path, senior_typing):
    import rebuild_donor_from_revised_export as rebuild
    cells, context = _tables(tmp_path)
    a.write_tsv(cells, context, tmp_path / "6539.tsv", _chosen("all"), "all", scale=1.0)
    frame = rebuild.read_export(tmp_path / "6539.tsv")
    assert rebuild.read_first_image(tmp_path / "6539.tsv") == "img - resolution #1"
    assert rebuild.islet_number(frame["parent_raw"]).tolist()[:3] == [4, 100_002, 200_009]
    is_pln, census = rebuild.classify_parents(frame["parent_raw"])
    assert is_pln.tolist() == [False, False, False, False, True] and census["islet_parent_rows"] == 3
    assert rebuild.region(frame["signed_dist_islet_um"]).tolist()[:4] == ["core", "core", "core", "expansion"]
