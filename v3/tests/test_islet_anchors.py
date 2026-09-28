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
    rows = ia.anchors_for_image(im, manual_dir=gj, resolution_um=2.0, n_l0_tiles=16, l0_tile=32)
    by = {r["channel"]: r for r in rows}
    assert set(by) == {"CD99", "INS", "GCG", "SST"}
    assert by["CD99"]["level_um"] == pytest.approx(2.0)
    assert by["INS"]["exposure_ms"] == pytest.approx(75.0)
    assert by["SST"]["biomarker_match"] is True
    # rows 2..61 x columns 12..61 at 4x downsample: Tissue from x = 48, and 5 um (10 px) inside the image edge
    assert all(r["n_pixels"] == 60 * 50 for r in rows)
    assert len({r["annotation_hash"] for r in rows}) == 1 and len({r["params_hash"] for r in rows}) == 1
    assert by["INS"]["excluded_fraction"] > 0.05          # the islet (+25 um) is excluded in pass 2
    for ch in by:
        # blurred level-2 block means of Poisson(bg): the anchor sits just above the tissue level,
        # and the bright "islet" does not move it
        assert tissue_bg[ch] < by[ch]["anchor"] < 1.5 * tissue_bg[ch], (ch, by[ch]["anchor"])
        assert np.isfinite(by[ch]["bg_smo_mean"])
        assert by[ch]["clamp_fraction"] == 0.0 and by[ch]["zero_fraction_l0"] < 0.01   # Poisson(>=30): no zeros
    out = ia.upsert_csv(rows, tmp_path / "a.csv")
    assert len(out) == 4 and out["anchor"].notna().all()
    # --skip-existing: same export -> skipped; an added Artifact changes the hash -> recomputed
    cur = {name: ia.current_annotation_hash(im, gj, None)}
    phash = rows[0]["params_hash"]
    assert ia.existing_images(tmp_path / "a.csv", ia.METHOD_VERSION, ia.CHANNELS, cur, phash) == {name}
    feats.append({"type": "Feature", "geometry": {"type": "Polygon", "coordinates": [
        [[150, 150], [200, 150], [200, 200], [150, 200], [150, 150]]]},
        "properties": {"classification": {"name": "Artifact"}}})
    (gj / f"{ia.image_tag(name)}.geojson").write_text(json.dumps({"type": "FeatureCollection", "features": feats}))
    cur2 = {name: ia.current_annotation_hash(im, gj, None)}
    assert cur2[name] != cur[name]
    assert ia.existing_images(tmp_path / "a.csv", ia.METHOD_VERSION, ia.CHANNELS, cur2, phash) == set()
    rows2 = ia.anchors_for_image(im, manual_dir=gj, resolution_um=2.0, with_l0=False, with_smo=False)
    assert rows2[0]["annotation_hash"] == cur2[name]
    assert rows2[0]["n_pixels"] == 60 * 50 - 13 * 13      # Artifact covers level pixels 37..49 x 37..49
    assert rows2[0]["n_excluded_polygons"] == 1


# ---------------------------------------------------------------- censored background model

def _clamped(mu, sigma, n, seed=0):
    y = np.random.default_rng(seed).normal(mu, sigma, n)
    return np.maximum(0, np.rint(y)).astype(np.int64)


@pytest.mark.parametrize("mu,sigma", [(-27.0, 27.0), (-5.0, 40.0), (20.0, 35.0), (90.0, 50.0)])
def test_clamped_normal_recovers_latent_background(mu, sigma):
    x = _clamped(mu, sigma, 2_000_000, seed=int(abs(mu)))
    r = ia.fit_clamped_normal(x)
    assert r["clamp_converged"]
    assert r["clamp_mu"] == pytest.approx(mu, abs=0.03 * sigma + 0.3)
    assert r["clamp_sigma"] == pytest.approx(sigma, rel=0.03)
    from scipy.special import ndtr
    assert r["clamp_fraction"] == pytest.approx(ndtr((0.5 - mu) / sigma), abs=0.002)
    assert r["clamp_fit_ks"] < 0.005


def test_clamped_normal_ignores_bright_signal_tail():
    rng = np.random.default_rng(3)
    bg = _clamped(-20.0, 30.0, 1_500_000, seed=4)
    signal = np.rint(rng.lognormal(np.log(600), 0.6, 80_000)).astype(np.int64)   # ~5% bright pixels
    r0, r1 = ia.fit_clamped_normal(bg), ia.fit_clamped_normal(np.r_[bg, signal])
    assert r1["clamp_mu"] == pytest.approx(r0["clamp_mu"], abs=1.5)
    assert r1["clamp_sigma"] == pytest.approx(r0["clamp_sigma"], rel=0.05)


def test_clamped_normal_is_scale_equivariant():
    y = np.random.default_rng(8).normal(-10.0, 25.0, 1_000_000)
    r1 = ia.fit_clamped_normal(np.maximum(0, np.rint(y)))
    r3 = ia.fit_clamped_normal(np.maximum(0, np.rint(3 * y)))
    assert r3["clamp_mu"] == pytest.approx(3 * r1["clamp_mu"], abs=1.0)
    assert r3["clamp_sigma"] == pytest.approx(3 * r1["clamp_sigma"], rel=0.02)
    assert r3["clamp_fraction"] == pytest.approx(r1["clamp_fraction"], abs=0.005)


def test_clamped_normal_unidentifiable_or_small_gives_nan():
    r = ia.fit_clamped_normal(_clamped(-60.0, 20.0, 500_000))      # ~99.9% zeros
    assert r["clamp_fraction"] > ia.CLAMP_MAX_FRACTION and np.isnan(r["clamp_mu"]) and np.isnan(r["clamp_sigma"])
    r = ia.fit_clamped_normal(_clamped(0.0, 10.0, 100))
    assert np.isnan(r["clamp_mu"]) and np.isnan(r["clamp_fraction"])


def test_positivity_threshold_is_max_of_fold_and_floor():
    assert ia.positivity_threshold(10.0, 5.0) == 50.0
    assert ia.positivity_threshold(10.0, 5.0, -20.0, 30.0, 3.0) == pytest.approx(70.0)   # floor binds
    assert ia.positivity_threshold(10.0, 5.0, -20.0, 30.0, 2.0) == pytest.approx(50.0)   # fold binds (40 < 50)
    assert ia.positivity_threshold(10.0, 5.0, -20.0, 30.0, 0.0) == 50.0                  # k = 0: floor off
    assert ia.positivity_threshold(10.0, 5.0, np.nan, np.nan, 3.0) == 50.0               # no model: floor off


# ---------------------------------------------------------------- two passes: the halo around islets

def _distance_to_centres(shape, centres):
    yy, xx = np.mgrid[0:shape[0], 0:shape[1]]
    d = np.full(shape, np.inf)
    for cy, cx in centres:
        d = np.minimum(d, np.hypot(yy - cy, xx - cx))
    return d


def test_two_pass_anchor_excludes_islet_halo():
    """Signal spill-over decays away from each islet (here 2x the background at the islet edge, 15 um
    decay length, out to 30 um): it fattens the right shoulder of the background peak, so a one-pass
    half-max rises; pass 2 drops the 25 um-dilated endocrine mask of pass 1 and recovers the exocrine
    background."""
    rng = np.random.default_rng(12)
    shape = (600, 600)                                   # 1 um pixels
    centres = [(100, 100), (100, 300), (100, 500), (300, 200), (300, 400), (500, 100), (500, 300), (500, 500)]
    d = _distance_to_centres(shape, centres)
    core, halo = d <= 20, (d > 20) & (d <= 50)           # islets r = 20 um
    bg = lambda s: rng.lognormal(np.log(s), 0.25, shape)
    ins = bg(40.0)
    ins[halo] *= 1 + np.exp(-(d[halo] - 20) / 15.0)
    ins[core] = 3000.0
    blurred = {"CD99": bg(300.0), "INS": ins, "GCG": bg(120.0), "SST": bg(80.0)}
    ones = np.ones(shape, bool)
    truth = ia.log_background(bg(40.0).ravel())["bg_halfmax"]
    res = ia.compute_anchors(blurred, ones, ones, 1.0, None, two_pass=True, excl_dilate_um=25.0)
    one = ia.compute_anchors(blurred, ones, ones, 1.0, None, two_pass=False)
    assert one["INS"]["anchor"] > 1.03 * truth                       # the halo inflates a one-pass anchor
    assert res["INS"]["anchor_pass1"] == pytest.approx(one["INS"]["anchor"])
    assert res["INS"]["anchor"] == pytest.approx(truth, rel=0.01)     # pass 2 recovers the exocrine background
    assert res["_meta"]["excluded_fraction"] == pytest.approx((d <= 45).mean(), abs=0.04)
    assert res["_meta"]["n_sampled_pass2"] < res["_meta"]["n_sampled"]
    for ch in ("CD99", "GCG", "SST"):                                  # channels without a halo do not move
        assert res[ch]["anchor"] == pytest.approx(one[ch]["anchor"], rel=0.01)


def test_dilate_um_covers_the_radius_conservatively():
    m = np.zeros((200, 200), bool)
    m[100, 100] = True
    d = ia.dilate_um(m, 25.0, 1.0, block=4)
    yy, xx = np.mgrid[0:200, 0:200]
    r = np.hypot(yy - 100, xx - 100)
    assert d[r <= 25].all()                              # every pixel within 25 um is covered
    assert not d[r > 25 + 2 * 4 * np.sqrt(2)].any()      # conservative by at most ~ one block diagonal x 2


def test_endocrine_mask_uses_the_floor():
    f = {"CD99": np.array([[0.0, 0.0]]), "INS": np.array([[60.0, 80.0]])}
    st = {"CD99": dict(anchor=100.0, clamp_mu=np.nan, clamp_sigma=np.nan),
          "INS": dict(anchor=10.0, clamp_mu=-20.0, clamp_sigma=30.0)}
    assert ia.endocrine_mask(f, st, 3.0, 5.0, 0.0).tolist() == [[True, True]]     # 60, 80 >= 50
    assert ia.endocrine_mask(f, st, 3.0, 5.0, 3.0).tolist() == [[False, True]]    # floor -20 + 90 = 70


def test_occupied_threshold_finds_the_empty_space_valley():
    rng = np.random.default_rng(21)
    empty = rng.lognormal(np.log(15), 0.3, 300_000)
    tissue = rng.lognormal(np.log(400), 0.3, 700_000)
    thr = ia.occupied_threshold(np.r_[empty, tissue])
    assert 30 < thr < 200
    assert ia.occupied_threshold(tissue) == 0.0          # unimodal: everything is occupied


# ---------------------------------------------------------------- analysis mask: union of Tissue minus pLN/Artifact

def test_analysis_mask_unions_tissue_and_subtracts_excluded_classes():
    shapely = pytest.importorskip("shapely")
    from shapely.geometry import Polygon
    from shapely.ops import unary_union
    t1, t2 = _square(10, 10, 110, 90), _square(60, 40, 180, 150)       # two overlapping Tissue polygons
    pln = _square(150, 120, 200, 200)                                  # a pLN partly outside Tissue
    art = np.array([[20, 20], [50, 25], [35, 60], [20, 20]], float)     # an Artifact inside Tissue
    tissue = [("Tissue", [t1]), ("Tissue", [t2])]
    excluded = [("pLN", [pln]), ("Artifact", [art])]
    m = ia.analysis_mask(tissue, excluded, (200, 200), 1.0, 1.0)
    yy, xx = np.mgrid[0:200, 0:200]
    ref_geom = unary_union([Polygon(t1), Polygon(t2)]).difference(unary_union([Polygon(pln), Polygon(art)]))
    inside = shapely.contains_xy(ref_geom, xx + 0.5, yy + 0.5)
    covered = shapely.intersects_xy(ref_geom, xx + 0.5, yy + 0.5)
    assert ((m == inside) | (m == covered)).all()         # only pixel centres ON an edge may go either way
    assert m[65, 80]                                      # the Tissue overlap is kept (even-odd would cancel it)
    assert not m[130, 170] and not m[35, 35]              # pLN and Artifact are cut out
    # even-odd over all rings (v1.0) cancels the overlap
    assert not ia.rasterize_even_odd([t1, t2], (200, 200))[65, 80]
    # downsampled level + 5 um edge margin (pixel centres), as the detector's inner rectangle
    m2 = ia.analysis_mask([("Tissue", [_square(0, 0, 200, 200)])], [], (100, 100), 2.0, 2.0, (200, 200), 10.0)
    assert m2[5:95, 5:95].all() and not m2[:4].any() and not m2[:, 96:].any()


def test_polygons_area_counts_overlaps_once():
    pytest.importorskip("shapely")
    a = ia.polygons_area_um2([("Tissue", [_square(0, 0, 10, 10)]), ("Tissue", [_square(5, 0, 15, 10)])], 0.5)
    assert a == pytest.approx(150 * 0.25)


# ---------------------------------------------------------------- annotation hash (detector parity)

def _fc(*feats):
    return {"type": "FeatureCollection", "features": [
        {"type": "Feature", "geometry": g, "properties": {"classification": {"name": c}}} for c, g in feats]}


def test_annotation_hash_format_is_pinned(tmp_path):
    # The Groovy port (1_panc_objects_v3.groovy, annotationHash) must produce exactly this text.
    sq = [[0, 0], [10.5, 0], [10.5, 10.25], [0, 10.25], [0, 0]]
    f = tmp_path / "a.geojson"
    f.write_text(json.dumps(_fc(("Tissue", {"type": "Polygon", "coordinates": [sq]}),
                                ("pLN", {"type": "MultiPolygon", "coordinates": [[[[20, 20], [30, 20], [30, 30], [20, 20]]]]}),
                                ("Islet", {"type": "Polygon", "coordinates": [sq]}))))
    polys = ia.read_annotation_polygons(f, ("Tissue", "pLN", "Artifact"))
    assert [c for c, _ in polys] == ["Tissue", "pLN"]
    lines = ia.annotation_signature(polys)
    assert lines == ["Tissue|0.00,0.00;10.50,0.00;10.50,10.25;0.00,10.25;0.00,0.00",
                     "pLN|20.00,20.00;30.00,20.00;30.00,30.00;20.00,20.00"]
    import hashlib
    text = "classes=Artifact,Tissue,pLN\n" + "\n".join(lines)
    assert ia.annotation_hash(polys, ("Tissue", "pLN", "Artifact")) == hashlib.sha256(text.encode()).hexdigest()[:16]


def test_annotation_hash_sensitivity_and_invariance():
    t = ("Tissue", [_square(0, 0, 100, 100), _square(40, 40, 60, 60)])
    p = ("pLN", [_square(200, 0, 250, 50)])
    cls = ("Tissue", "pLN", "Artifact")
    h = ia.annotation_hash([t, p], cls)
    assert ia.annotation_hash([p, t], cls) == h                                  # feature order
    moved = ("Tissue", [_square(0, 0, 100, 100.01), _square(40, 40, 60, 60)])
    assert ia.annotation_hash([moved, p], cls) != h                              # a vertex moved by 0.01 px
    art = ("Artifact", [_square(10, 10, 20, 20)])
    assert ia.annotation_hash([t, p, art], cls) != h                             # an Artifact was drawn
    assert ia.annotation_hash([t, p], ("Tissue", "pLN")) != h                    # the class list is part of it
    tiny = ("Tissue", [_square(0, 0, 100, 100) + 0.004, _square(40, 40, 60, 60)])
    assert ia.annotation_hash([tiny, p], cls) == h                               # below the export's 2 decimals
