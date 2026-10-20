"""Tests for scripts/v3/senior/pixel_levels.py and cell_mean_check.py (synthetic data; never touches the project).

Run:  <Senior>/.venv/bin/python -m pytest -q scripts/v3/tests/test_preprocessing_audit.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "senior"))
import cell_mean_check as c  # noqa: E402
import pixel_levels as p  # noqa: E402


# ----------------------------------------------------------------------------------------------- stored pixel levels

def test_the_histogram_gives_the_exact_levels_of_the_pixels_it_counts():
    rng = np.random.default_rng(0)
    values = np.concatenate([np.zeros(4000, dtype=np.uint16), rng.integers(1, 3000, 6000).astype(np.uint16), np.array([65535], dtype=np.uint16)])
    got = p.from_histogram(np.bincount(values, minlength=65536))
    assert got["n_pixels"] == len(values)
    assert got["frac_zero"] == pytest.approx(4000 / len(values))
    assert got["max"] == 65535 and got["mean"] == pytest.approx(values.mean())
    for q in p.PERCENTILES:                                # the smallest stored value with that share of pixels at or below it
        assert got[f"p{q:g}"] == np.percentile(values, q, method="inverted_cdf")


def test_an_empty_region_reports_no_levels_rather_than_zeros():
    assert p.from_histogram(np.zeros(65536, dtype=np.int64)) == {"n_pixels": 0}


def test_tiles_are_spread_through_the_region_and_the_tag_is_the_exports_own():
    tiles = [(i, 0) for i in range(100)]
    picked = p.evenly(tiles, 5)
    assert picked == [(0, 0), (25, 0), (50, 0), (74, 0), (99, 0)] and p.evenly(tiles[:3], 5) == tiles[:3]
    assert p.image_tag("0080502_6580_pLNpanc_Scan1.er.qptiff - resolution #1") == "0080502_6580_pLNpanc_Scan1.er.qptiff_-_resolution_1"


# ----------------------------------------------------------------------------------------------- cell means against pixels

def _feature(i: int, name: str = "cell") -> dict:
    return {"type": "Feature", "id": f"id-{i}", "geometry": {"type": "Polygon", "coordinates": [[[i, 0], [i + 4, 0], [i + 4, 4], [i, 4], [i, 0]]]},
            "properties": {"objectType": "cell", "name": name}}


def test_whole_features_are_read_from_the_middle_of_a_large_export(tmp_path):
    feats = [_feature(0), _feature(1, name='a "quoted" } brace'), _feature(2, name="{ open"), _feature(3), _feature(4)]
    text = '{"type":"FeatureCollection","features":[' + ",".join(json.dumps(f, separators=(",", ":")) for f in feats) + "]}"
    path = tmp_path / "cells.geojson"
    path.write_text(text)
    assert [f["id"] for f in c.features_at(path, 0, 10)] == [f"id-{i}" for i in range(5)]
    # braces and escaped quotes inside a string do not end a feature early
    assert c.features_at(path, 0, 10)[1]["properties"]["name"] == 'a "quoted" } brace'
    # starting inside feature 1: the first WHOLE feature after the offset is feature 2
    inside = text.index('"id-1"')
    assert [f["id"] for f in c.features_at(path, inside, 2)] == ["id-2", "id-3"]


def test_a_feature_cut_off_by_the_end_of_the_chunk_is_not_returned(tmp_path, monkeypatch):
    feats = [_feature(i) for i in range(3)]
    text = '{"type":"FeatureCollection","features":[' + ",".join(json.dumps(f, separators=(",", ":")) for f in feats) + "]}"
    path = tmp_path / "cells.geojson"
    path.write_text(text)
    monkeypatch.setattr(c, "CHUNK", text.index('"id-2"') + 20)             # the chunk ends inside the third feature
    assert [f["id"] for f in c.features_at(path, 0, 10)] == ["id-0", "id-1"]


def test_the_check_reads_each_marker_from_the_page_of_its_own_name(tmp_path):
    """Two pages in an order the measurement table does not share: a mean taken from the wrong page would show."""
    tifffile = pytest.importorskip("tifffile")
    pytest.importorskip("zarr")
    rng = np.random.default_rng(1)
    planes = {"DAPI": rng.integers(0, 500, (64, 64)).astype(np.uint16), "CD8": rng.integers(0, 500, (64, 64)).astype(np.uint16)}
    qptiff = tmp_path / "image.qptiff"
    with tifffile.TiffWriter(qptiff) as tw:
        for name, plane in planes.items():
            tw.write(plane, tile=(32, 32), description=f'<?xml version="1.0"?><x><Biomarker>{name}</Biomarker></x>', metadata=None)
    cells = {"a": (2, 3, 12, 14), "b": (20, 22, 31, 40), "c": (40, 5, 55, 21)}                    # x0, y0, x1, y1
    feats = [{"type": "Feature", "id": k, "geometry": {"type": "Polygon", "coordinates": [[[x0, y0], [x1, y0], [x1, y1], [x0, y1], [x0, y0]]]},
              "properties": {"objectType": "cell"}} for k, (x0, y0, x1, y1) in cells.items()]
    geojson = tmp_path / "cells.geojson"
    geojson.write_text('{"type":"FeatureCollection","features":[' + ",".join(json.dumps(f, separators=(",", ":")) for f in feats) + "]}")
    table = pd.DataFrame({"Object ID": list(cells),
                          "Cell: DAPI: Mean": [planes["DAPI"][y0:y1, x0:x1].mean() for (x0, y0, x1, y1) in cells.values()],
                          "Cell: CD8: Mean": [planes["CD8"][y0:y1, x0:x1].mean() for (x0, y0, x1, y1) in cells.values()]})
    parquet = tmp_path / "cells.parquet"
    table.to_parquet(parquet, index=False)

    rows = {r["marker"]: r for r in c.one_image(("d", str(qptiff), str(geojson), str(parquet), 24, False))}
    assert set(rows) == {"DAPI", "CD8"} and rows["CD8"]["page"] == 1
    for r in rows.values():
        assert r["n_cells"] == 3 and r["median_ratio_qupath_over_pixels"] == pytest.approx(1.0) and r["share_within_2pct"] == 1.0
        assert r["median_abs_difference"] == pytest.approx(0.0, abs=1e-9)

    # the table's CD8 column filled from the DAPI page: the check must say so
    table["Cell: CD8: Mean"] = table["Cell: DAPI: Mean"]
    table.to_parquet(parquet, index=False)
    wrong = {r["marker"]: r for r in c.one_image(("d", str(qptiff), str(geojson), str(parquet), 24, False))}
    assert wrong["CD8"]["share_within_2pct"] < 1.0 and wrong["DAPI"]["share_within_2pct"] == 1.0
