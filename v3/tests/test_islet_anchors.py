"""Tests for scripts/v3/islet_anchors.py (synthetic data only; never touches the project).

Run:  /scratch/io60/venv-v3/bin/python -m pytest -q scripts/v3/tests/test_islet_anchors.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import islet_anchors as ia  # noqa: E402

RNG = np.random.default_rng(1234)


# ---------------------------------------------------------------- names, channels, levels

def test_image_tag_matches_export_convention():
    assert ia.image_tag("6374_Scan1.er.qptiff - resolution #1") == "6374_Scan1.er.qptiff_-_resolution_1"
    assert ia.image_tag("a/b c") == "a_b_c"


def test_resolve_channel_exact_alias_and_forbidden():
    names = ["DAPI", "Insulin (C02)", "CD99", "Glucagon.tif", "SST", "IAPP", "HLA-A"]
    assert ia.resolve_channel(names, "CD99") == 2
    assert ia.resolve_channel(names, "INS") == 1      # alias + "(Cxx)" suffix stripped
    assert ia.resolve_channel(names, "GCG") == 3      # alias + ".tif" stripped
    with pytest.raises(ValueError):
        ia.resolve_channel(names, "IAPP")             # failed QC
    with pytest.raises(ValueError):
        ia.resolve_channel(names, "HLA-A")            # disease-associated
    with pytest.raises(KeyError):
        ia.resolve_channel(["DAPI"], "CD99")


def test_choose_level_picks_closest_to_2um():
    full = (51120, 35520)
    shapes = [full, (25560, 17760), (12780, 8880), (6390, 4440), (3195, 2220), (1597, 1110)]
    li, um, dx, dy = ia.choose_level(shapes, full, 0.5078, 2.0)
    assert li == 2 and dx == pytest.approx(4.0) and um == pytest.approx(2.0312, abs=1e-4)


# ---------------------------------------------------------------- rasterisation

def _square(x0, y0, x1, y1):
    return np.array([[x0, y0], [x1, y0], [x1, y1], [x0, y1], [x0, y0]], float)


def test_rasterize_square_pixel_centres():
    m = ia.rasterize_even_odd([_square(10, 10, 30, 20)], (40, 40))
    assert m.sum() == 20 * 10
    assert m[10, 10] and m[19, 29] and not m[20, 10] and not m[10, 30]
    # at downsample 2 the square covers level pixels 5..14 x 5..9
    m2 = ia.rasterize_even_odd([_square(10, 10, 30, 20)], (20, 20), 2.0, 2.0)
    assert m2.sum() == 10 * 5


def test_rasterize_hole_island_and_multipolygon_match_shapely():
    shapely = pytest.importorskip("shapely")
    from shapely.geometry import MultiPolygon, Polygon
    outer = _square(5, 5, 95, 95)
    hole = _square(30, 30, 70, 70)
    island = _square(40, 40, 60, 60)          # a separate polygon inside the hole
    other = np.array([[100, 10], [140, 20], [120, 60], [100, 10]], float)
    geom = MultiPolygon([Polygon(outer, [hole]), Polygon(island), Polygon(other)])
    m = ia.rasterize_even_odd([outer, hole, island, other], (100, 150))
    yy, xx = np.mgrid[0:100, 0:150]
    ref = shapely.contains_xy(geom, xx + 0.5, yy + 0.5)
    assert np.array_equal(m, ref)


def test_rasterize_random_star_matches_shapely():
    shapely = pytest.importorskip("shapely")
    from shapely.geometry import Polygon
    ang = np.sort(RNG.uniform(0, 2 * np.pi, 40))
    rad = RNG.uniform(20, 60, 40)
    ring = np.c_[80 + rad * np.cos(ang) + 0.123, 70 + rad * np.sin(ang) + 0.377]
    ring = np.vstack([ring, ring[:1]])
    m = ia.rasterize_even_odd([ring], (150, 160))
    yy, xx = np.mgrid[0:150, 0:160]
    ref = shapely.contains_xy(Polygon(ring).buffer(0), xx + 0.5, yy + 0.5)
    assert (m != ref).sum() <= 2   # only pixel centres lying on an edge may differ


# ---------------------------------------------------------------- QPTIFF metadata

QPI_XML = """<?xml version="1.0" encoding="utf-16"?>
<PerkinElmer-QPI-ImageDescription>
  <DescriptionVersion>6</DescriptionVersion>
  <ImageType>FullResolution</ImageType>
  <ExposureTime>{exp}</ExposureTime>
  <SignalUnits>69</SignalUnits>
  <Name>ATTO 550</Name>
  <Responsivity>
    <Filter>
      <Name>ATTO 550</Name>
      <Response>215.656290166283</Response>
      <Date>2026-06-11T17:53:27.2379148Z</Date>
    </Filter>
  </Responsivity>
  <Biomarker>{bio}</Biomarker>
  <ScanProfile>{{"channels":[{{"name":"X","defaultExposure":150}}]}}</ScanProfile>
</PerkinElmer-QPI-ImageDescription>"""


def test_parse_qpi_description_units_and_fields():
    d = ia.parse_qpi_description(QPI_XML.format(exp=75000, bio="INS"))
    assert d["exposure_ms"] == pytest.approx(75.0)          # microseconds -> ms
    assert d["xml_biomarker"] == "INS"
    assert d["responsivity"] == pytest.approx(215.656290166283)
    assert d["filter_date"] == "2026-06-11"
    assert d["signal_units"] == 69
    # bytes in UTF-16 decode too; non-QPI descriptions give nothing
    assert ia.parse_qpi_description(QPI_XML.format(exp=1000, bio="DAPI").encode("utf-16"))["exposure_ms"] == 1.0
    assert ia.parse_qpi_description('{"shape": [1, 2]}') == {}
    assert ia.parse_qpi_description(None) == {}


# ---------------------------------------------------------------- SMO port

def test_smo_null_selects_five_percent_of_pure_noise():
    thr = ia.smo_null_rv().ppf(0.05)
    assert 0.02 < thr < 0.1
    noise = np.random.default_rng(7).uniform(size=(512, 512))
    s, _ = ia.smo_operator(noise)
    assert 0.035 < np.mean(s <= thr) < 0.065


def test_smo_mask_propagation_matches_package_semantics():
    # one masked pixel masks its whole 9x9 Gaussian window, then the 7x7 averaging window
    img = RNG.uniform(size=(40, 40))
    mask = np.zeros_like(img, bool)
    mask[20, 20] = True
    _, m = ia.smo_operator(img, mask)
    assert m[20, 20] and m[20 + 7, 20] and not m[20 + 8, 20]


def _blobby(shape=(256, 256), bg=20.0, amp=400.0, n=25, seed=3):
    rng = np.random.default_rng(seed)
    yy, xx = np.mgrid[0:shape[0], 0:shape[1]]
    sig = np.zeros(shape)
    for _ in range(n):
        cy, cx = rng.uniform(20, shape[0] - 20), rng.uniform(20, shape[1] - 20)
        sig += amp * np.exp(-((yy - cy) ** 2 + (xx - cx) ** 2) / (2 * 6.0 ** 2))
    return rng.poisson(bg + sig).astype(np.uint16)


def test_smo_background_recovers_flat_level():
    img = _blobby()
    tissue = np.ones(img.shape, bool)
    r = ia.smo_background(img, None, tissue & (img > 0))
    assert r["smo_bg_fraction"] > 0.01
    assert r["bg_smo_mean"] == pytest.approx(20.0, rel=0.15)


def test_smo_background_sample_excludes_zeros():
    rng = np.random.default_rng(5)
    img = np.clip(np.round(rng.normal(0, 6, (256, 256))), 0, None).astype(np.uint16)  # ~50% zeros
    tissue = np.ones(img.shape, bool)
    r = ia.smo_background(img, None, tissue & (img > 0))
    assert np.mean(img == 0) > 0.4
    assert r["bg_smo_mean"] > 1.0        # a zero-including mean would be far lower
    assert r["bg_smo_median"] >= 1.0


# ---------------------------------------------------------------- log-space background

def test_log_density_dequantises_low_integer_counts():
    # clamped half-normal noise, sigma 5: density of log x peaks at x = sigma
    x = np.round(np.abs(np.random.default_rng(11).normal(0, 5.0, 400_000))).astype(np.int64)
    lb = ia.log_background(x)
    assert lb["bg_logmode"] == pytest.approx(5.0, rel=0.3)
    assert lb["n_peaks"] == 1


def test_upper_substantial_peak_is_the_tissue_mode():
    rng = np.random.default_rng(2)
    empty = rng.lognormal(np.log(8), 0.35, 450_000)
    tissue = rng.lognormal(np.log(150), 0.30, 550_000)
    x = np.round(np.r_[empty, tissue]).astype(np.int64)
    lb = ia.log_background(x)
    assert lb["n_peaks"] == 2
    assert lb["bg_logmode"] == pytest.approx(150, rel=0.1)
    assert lb["logmode_lower"] == pytest.approx(8, rel=0.15)
    assert lb["lower_mass_fraction"] == pytest.approx(0.45, abs=0.05)
    assert lb["bg_halfmax"] > lb["bg_logmode"]


def test_halfmax_of_lognormal():
    mu, s = np.log(200.0), 0.4
    x = np.round(np.random.default_rng(3).lognormal(mu, s, 600_000)).astype(np.int64)
    lb = ia.log_background(x)
    s_eff = np.sqrt(s ** 2 + ia.LOG_SMOOTH ** 2)      # density is smoothed by LOG_SMOOTH
    assert lb["bg_halfmax"] == pytest.approx(np.exp(mu + s_eff * np.sqrt(2 * np.log(2))), rel=0.05)


def test_anchor_is_scale_equivariant():
    base = np.random.default_rng(4).lognormal(np.log(60), 0.5, 500_000)
    a1 = ia.log_background(np.round(base).astype(np.int64))["bg_halfmax"]
    a3 = ia.log_background(np.round(3 * base).astype(np.int64))["bg_halfmax"]
    assert a3 / a1 == pytest.approx(3.0, rel=0.03)


def test_anchor_ignores_a_few_percent_of_bright_endocrine_pixels():
    rng = np.random.default_rng(6)
    bg = rng.lognormal(np.log(40), 0.45, 600_000)
    islets = rng.lognormal(np.log(3000), 0.5, 18_000)          # 3% of tissue
    a0 = ia.log_background(np.round(bg).astype(np.int64))["bg_halfmax"]
    a1 = ia.log_background(np.round(np.r_[bg, islets]).astype(np.int64))["bg_halfmax"]
    assert a1 / a0 == pytest.approx(1.0, abs=0.03)


# ---------------------------------------------------------------- CSV upsert

def test_upsert_csv_is_idempotent(tmp_path):
    out = tmp_path / "anchors.csv"
    r = lambda img, ch, a: dict(image=img, tag=ia.image_tag(img), channel=ch, anchor=a)
    df = ia.upsert_csv([r("A", "CD99", 1.0), r("A", "INS", 2.0)], out)
    assert len(df) == 2 and list(df.columns[:13]) == ia.COLUMNS[:13]
    df = ia.upsert_csv([r("A", "INS", 5.0)], out)
    assert len(df) == 2 and df.set_index("channel").loc["INS", "anchor"] == 5.0
    df = ia.upsert_csv([r("B", "CD99", 3.0)], out)
    assert len(df) == 3
    df = ia.upsert_csv([r("B", "CD99", 3.0)], out)
    assert len(df) == 3
    assert ia.existing_images(out, method_version="nope") == set()


# ---------------------------------------------------------------- end to end on a synthetic QPTIFF

def _write_qptiff(path: Path, chans: list[np.ndarray], names: list[str], exposures_us: list[int]):
    tifffile = pytest.importorskip("tifffile")
    levels = 3
    with tifffile.TiffWriter(path) as tw:
        def page(a, bio, exp, kind):
            desc = QPI_XML.format(exp=exp, bio=bio).replace("FullResolution", kind)
            tw.write(a, tile=(16, 16), software="PerkinElmer-QPI", description=desc,
                     metadata=None, contiguous=False)
        for a, n, e in zip(chans, names, exposures_us):
            page(a, n, e, "FullResolution")
        tw.write(np.zeros((8, 8, 3), np.uint8), software="PerkinElmer-QPI", metadata=None,
                 description=QPI_XML.format(exp=0, bio="thumb").replace("FullResolution", "Thumbnail"))
        for lvl in range(1, levels):
            f = 2 ** lvl
            for a, n, e in zip(chans, names, exposures_us):
                red = a.reshape(a.shape[0] // f, f, a.shape[1] // f, f).mean(axis=(1, 3)).astype(np.uint16)
                page(red, n, e, "ReducedResolution")


def test_anchors_for_image_end_to_end(tmp_path):
    rng = np.random.default_rng(9)
    H = W = 256
    names = ["DAPI", "CD99", "INS", "GCG", "SST"]
    tissue_bg = {"DAPI": 300, "CD99": 120, "INS": 30, "GCG": 90, "SST": 60}
    chans = []
    for n in names:
        a = rng.poisson(tissue_bg[n], (H, W)).astype(np.float64)
        a[:, :40] = rng.poisson(3, (H, 40))            # outside the tissue polygon
        if n == "INS":
            a[100:130, 100:130] += 4000                # an "islet"
        chans.append(np.clip(a, 0, 65535).astype(np.uint16))
    tif = tmp_path / "synthetic.qptiff"
    _write_qptiff(tif, chans, names, [1000, 40000, 75000, 50000, 150000])
    gj = tmp_path / "manual"
    gj.mkdir()
    name = "synthetic.qptiff - resolution #1"
    ring = [[48, 0], [256, 0], [256, 256], [48, 256], [48, 0]]
    feats = [{"type": "Feature", "geometry": {"type": "Polygon", "coordinates": [ring]},
              "properties": {"classification": {"name": "Tissue"}}}]
    (gj / f"{ia.image_tag(name)}.geojson").write_text(json.dumps({"type": "FeatureCollection", "features": feats}))
    im = ia.ProjectImage(image=name, tag=ia.image_tag(name), entry_id=1, path=tif, channels=names,
                         pixel_um=0.5, width=W, height=H)
    rows = ia.anchors_for_image(im, manual_dir=gj, target_um=2.0, with_l0=False)
    by = {r["channel"]: r for r in rows}
    assert set(by) == {"CD99", "INS", "GCG", "SST"}
    assert by["CD99"]["level_um"] == pytest.approx(2.0)
    assert by["INS"]["exposure_ms"] == pytest.approx(75.0)
    assert by["SST"]["biomarker_match"] is True
    for ch in by:
        # level-2 block means of Poisson(bg): the anchor sits just above the tissue level,
        # and the bright "islet" does not move it
        assert tissue_bg[ch] < by[ch]["anchor"] < 1.5 * tissue_bg[ch], (ch, by[ch]["anchor"])
        assert by[ch]["n_pixels"] == 52 * 64                  # columns 12..63 at 4x downsample
        assert np.isfinite(by[ch]["bg_smo_mean"])
    out = ia.upsert_csv(rows, tmp_path / "a.csv")
    assert len(out) == 4 and out["anchor"].notna().all()
    # level-0 zero-fraction sampler on the same file (tiles must lie inside tissue)
    import tifffile
    with tifffile.TiffFile(tif) as tf:
        s = tf.series[0]
        tissue = ia.rasterize_even_odd(ia.read_annotation_rings(gj / f"{im.tag}.geojson"), (64, 64), 4, 4)
        zf = ia._sample_level0_zero_fraction(s, 2, tissue, 4.0, 4.0, n_tiles=4, tile=32)
        assert 0.0 <= zf < 0.01
