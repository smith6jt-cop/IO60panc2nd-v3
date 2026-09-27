#!/usr/bin/env python3
"""islet_anchors.py -- per-image, per-channel background anchors for the v3 islet detector.

The v3 detector (qupath/1_panc_objects_v3.groovy, ISLET_METHOD='drotar_cd99') thresholds
FOLDS, i.e. raw counts divided by one scalar per image and channel (the "anchor"). This
module estimates those scalars for CD99, INS, GCG and SST and writes them to
results/v3/islet_anchors.csv, one row per (image, channel).

Governing rule (Phenocycler_Analysis docs/INTENSITY_SCALE.md): *standardize the unit, never
the distribution*. An anchor is a property of the measurement (instrument scale, background,
censoring), estimated from the image's own pixels in log space and applied as a LINEAR divisor.
Nothing here matches distributions between images, and no anchor is a quantile of the
positive (endocrine) population.

What is measured, per channel, on the pyramid level closest to 2 um/px inside the manual
Tissue annotation (pLN is a separate annotation and is never part of Tissue):

  n_pixels, zero_fraction   in-tissue pixel count and the fraction that is exactly 0.
                            The scanner subtracts offset + autofluorescence and stores UINT16,
                            so every negative excursion is clamped to 0: a zero is "<= 0", not
                            "no signal". Zeros are EXCLUDED from every background estimate.
                            (The pyramid is block-averaged, so the ~2 um level has fewer exact
                            zeros than full resolution; zero_fraction_l0 is sampled at level 0.)
  bg_smo_mean, bg_smo_std   Silver Mountain Operator (SMO) background, ported from the `smo`
                            package that KINTSUGI's signal/smo.py wraps, used as
                            SMO(sigma=1.0, size=7), threshold 0.05, then bg_rv().mean()/.std().
                            The operator runs on the 2D level image (only out-of-tissue and
                            saturated pixels masked); zeros are dropped from the background
                            SAMPLE, not from the operator input -- SMO dilates every masked
                            pixel to a ~15x15 block, so masking 30-80% zeros would discard
                            nearly every pixel next to a zero and bias the sample bright.
  bg_logmode                mode of the log-intensity density of the nonzero tissue pixels
                            (integer counts spread uniformly over [k-0.5, k+0.5] before the log,
                            so low counts are not quantised). When the density is bimodal the
                            UPPER substantial peak is taken: the lower one is empty space inside
                            the Tissue outline (lumens, fat, gaps; ~5-13 counts for CD99/GCG),
                            the upper one is the exocrine tissue the endocrine signal must exceed.
  bg_halfmax                right half-maximum point of that peak (in log space): the tissue
                            background level PLUS its width. This is the ANCHOR (see below).

Anchor choice (ANCHOR_SOURCE = 'halfmax'):
  * SMO is recorded but not used. At ~2 um the level image is smooth relative to SMO's iid-noise
    null, so it accepts only ~1-2% of tissue pixels, roughly uniformly across intensity; the
    mean then includes islet signal (6450panc INS: mean 21 vs median 3), i.e. it depends on the
    endocrine population, which the governing rule forbids. At full resolution SMO does
    discriminate, but what it isolates is the empty-space noise floor, not the exocrine level.
  * The mode alone ignores background width: on 6450panc the SST exocrine background is broad
    (mode 5 counts, exocrine p90 ~7x the mode), so the same fold means very different false-
    positive rates on different channels. The right half-maximum scales with level and width
    together, stays a pure scale estimate (anchor(k*I) = k*anchor(I)), and is insensitive to
    the few-percent endocrine tail. The evaluation behind this choice is in ISLETS.md.

Exposure: ExposureTime (microseconds) is read from the channel's level-0 page XML
(PerkinElmer-QPI-ImageDescription) and written as exposure_ms, with the filter responsivity and
its calibration date (epochs 2025-09-10, 2026-03-12, 2026-06-11 differ by up to 1.67x in
exposure). The XML <Biomarker> is checked against the project channel name.

Usage:
  islet_anchors.py --images 6450panc 6374_Scan1         # substring match on image names
  islet_anchors.py --all [--skip-existing]
  islet_anchors.py --all --out /scratch/.../anchors.csv # anywhere; default results/v3/islet_anchors.csv

Idempotent: rows are upserted by (image, channel); rerunning an image replaces its rows, and
the values are deterministic (only computed_at changes).
"""
from __future__ import annotations

import argparse
import datetime as _dt
import html
import json
import logging
import os
import re
import sys
import tempfile
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

log = logging.getLogger("islet_anchors")

PROJ = Path("/home/smith6jt/IO60panc2nd")
DEFAULT_PROJECT = PROJ / "project.qpproj"
DEFAULT_MANUAL_DIR = PROJ / "results" / "v3" / "manual_annotations"
DEFAULT_OUT = PROJ / "results" / "v3" / "islet_anchors.csv"

CHANNELS = ("CD99", "INS", "GCG", "SST")
# Never used for detection: IAPP failed QC; HLA-A is disease-associated (rises in T1D exocrine).
FORBIDDEN_CHANNELS = frozenset({"IAPP", "HLA-A"})
# Same alias map as 1_panc_objects*.groovy, so a channel resolves to the same index in both.
CHANNEL_ALIASES = {
    "SST": ("SST", "Somatostatin", "Stt"),
    "GCG": ("GCG", "Glucagon", "Gcg", "Opal 520"),
    "INS": ("INS", "Insulin", "Ins", "Opal 690"),
    "CD99": ("CD99",),
}

METHOD_VERSION = "1.0"
TARGET_UM = 2.0

SMO_SIGMA = 1.0
SMO_SIZE = 7
SMO_THRESHOLD = 0.05          # probability under the SMO null (smo.SMO.bg_rv default)
SMO_NULL_SHAPE = (1024, 1024)  # null depends only on sigma/size; 1M samples fix its 5% quantile

LOG_STEP = 0.01               # log-density grid (natural-log units)
LOG_SMOOTH = 0.10             # Gaussian smoothing of the log-density (log units)
PEAK_MIN_HEIGHT = 0.25        # a "substantial" peak is >= 25% of the tallest one ...
PEAK_MIN_PROMINENCE = 0.05    # ... and stands out by >= 5% of it
ANCHOR_SOURCE = "halfmax"     # which estimator becomes the anchor column

L0_TILES = 32                 # level-0 tiles sampled for zero_fraction_l0
L0_TILE = 512
SATURATION = 65535

COLUMNS = [
    # required
    "image", "tag", "channel", "level_um", "n_pixels", "zero_fraction",
    "bg_smo_mean", "bg_smo_std", "bg_logmode", "anchor", "exposure_ms",
    "method_version", "computed_at",
    # provenance / diagnostics
    "anchor_source", "bg_halfmax", "bg_smo_median", "smo_bg_fraction",
    "logmode_global", "logmode_lower", "lower_mass_fraction", "n_peaks",
    "median_nonzero", "p99_nonzero", "zero_fraction_l0", "saturated_fraction",
    "channel_index", "xml_biomarker", "biomarker_match", "responsivity", "filter_date",
    "signal_units", "entry_id", "level_index", "image_path",
]


# ----------------------------------------------------------------------------------------
# Project / image helpers
# ----------------------------------------------------------------------------------------

def image_tag(name: str) -> str:
    """File tag used by every v3 export: runs of characters outside [A-Za-z0-9._-] -> '_'."""
    return re.sub(r"[^A-Za-z0-9._-]+", "_", name)


def normalize_channel_name(name: str | None) -> str | None:
    """Mirror 1_panc_objects.groovy: strip a trailing '(Cxx)' and '.tif(f)', trim."""
    if name is None:
        return None
    s = re.sub(r"\s*\(C\d+\)\s*$", "", name)
    s = re.sub(r"(?i)\.(tif|tiff)\s*$", "", s)
    return s.strip()


def resolve_channel(names: list[str], canonical: str) -> int:
    """Index of `canonical` in the project channel names (exact first, then aliases)."""
    if canonical in FORBIDDEN_CHANNELS:
        raise ValueError(f"channel {canonical} must never be used for islet detection")
    norm = [normalize_channel_name(n) for n in names]
    if canonical in norm:
        return norm.index(canonical)
    for alias in CHANNEL_ALIASES.get(canonical, ()):
        for i, n in enumerate(norm):
            if n is not None and n.lower() == alias.lower():
                return i
    raise KeyError(f"channel {canonical!r} not found in {norm}")


@dataclass
class ProjectImage:
    image: str
    tag: str
    entry_id: int
    path: Path
    channels: list[str]
    pixel_um: float
    width: int
    height: int
    extra: dict = field(default_factory=dict)


def read_project(project: Path = DEFAULT_PROJECT) -> list[ProjectImage]:
    """Entries of a QuPath project with channel names and pixel size from serverBuilder.metadata."""
    with open(project, encoding="utf-8") as fh:
        proj = json.load(fh)
    out = []
    for e in proj["images"]:
        sb = e["serverBuilder"]
        md = sb["metadata"]
        uri = sb["uri"]
        path = Path(re.sub(r"^file:", "", uri))
        pc = md.get("pixelCalibration", {}).get("pixelWidth", {})
        out.append(ProjectImage(
            image=e["imageName"], tag=image_tag(e["imageName"]), entry_id=int(e["entryID"]),
            path=path, channels=[c["name"] for c in md["channels"]],
            pixel_um=float(pc.get("value", 0.5078)), width=int(md["width"]), height=int(md["height"])))
    return out


def select_images(images: list[ProjectImage], patterns: list[str] | None) -> list[ProjectImage]:
    """All images, or those whose name contains any of the patterns."""
    if not patterns:
        return images
    sel = [im for im in images if any(p in im.image for p in patterns)]
    missing = [p for p in patterns if not any(p in im.image for im in images)]
    if missing:
        raise SystemExit(f"no project image matches: {missing}")
    return sel


def choose_level(level_shapes: list[tuple[int, int]], full_shape: tuple[int, int],
                 pixel_um: float, target_um: float = TARGET_UM) -> tuple[int, float, float, float]:
    """Pyramid level whose pixel size is closest to target_um (in log ratio).

    Returns (level index, level pixel size um, downsample x, downsample y). The downsample is
    taken per axis from the actual dimensions, since odd sizes are floored by the scanner.
    """
    best = None
    for i, (h, w) in enumerate(level_shapes):
        dx = full_shape[1] / w
        dy = full_shape[0] / h
        um = pixel_um * np.sqrt(dx * dy)
        score = abs(np.log(um / target_um))
        if best is None or score < best[0] - 1e-12:
            best = (score, i, um, dx, dy)
    _, i, um, dx, dy = best
    return i, float(um), float(dx), float(dy)


# ----------------------------------------------------------------------------------------
# Tissue rasterisation (even-odd scanline; exact at pixel centres, O(perimeter + pixels))
# ----------------------------------------------------------------------------------------

def read_annotation_rings(geojson: Path, classes=("Tissue",)) -> list[np.ndarray]:
    """All rings (exteriors and holes) of the features whose classification is in `classes`.

    Coordinates are QuPath full-resolution pixels. Even-odd filling of ALL rings reproduces
    Polygon/MultiPolygon semantics for valid geometries (holes flip parity, islands flip back).
    """
    with open(geojson, encoding="utf-8") as fh:
        gj = json.load(fh)
    feats = gj["features"] if gj.get("type") == "FeatureCollection" else [gj]
    rings = []
    for ft in feats:
        cls = ((ft.get("properties") or {}).get("classification") or {}).get("name")
        if cls not in classes:
            continue
        g = ft["geometry"]
        if g["type"] == "Polygon":
            polys = [g["coordinates"]]
        elif g["type"] == "MultiPolygon":
            polys = g["coordinates"]
        else:
            continue
        for poly in polys:
            for ring in poly:
                rings.append(np.asarray(ring, dtype=np.float64))
    return rings


def rasterize_even_odd(rings: list[np.ndarray], shape: tuple[int, int],
                       dx: float = 1.0, dy: float = 1.0) -> np.ndarray:
    """Boolean mask of `shape`; pixel (r, c) is inside if its centre ((c+.5)dx, (r+.5)dy) is.

    `rings` are closed coordinate arrays in full-resolution pixels; dx, dy are the level
    downsamples. A crossing at level-x position xc switches pixels with c + 0.5 >= xc.
    """
    H, W = shape
    cnt = np.zeros((H, W + 1), dtype=np.int32)
    for ring in rings:
        if len(ring) < 3:
            continue
        x = ring[:, 0] / dx
        y = ring[:, 1] / dy
        if not (x[0] == x[-1] and y[0] == y[-1]):
            x = np.append(x, x[0])
            y = np.append(y, y[0])
        x0, y0, x1, y1 = x[:-1], y[:-1], x[1:], y[1:]
        keep = y0 != y1
        x0, y0, x1, y1 = x0[keep], y0[keep], x1[keep], y1[keep]
        ylo = np.minimum(y0, y1)
        yhi = np.maximum(y0, y1)
        # rows whose centre r + 0.5 lies in [ylo, yhi)
        r0 = np.clip(np.ceil(ylo - 0.5).astype(np.int64), 0, H)
        r1 = np.clip(np.ceil(yhi - 0.5).astype(np.int64), 0, H)
        n = r1 - r0
        tot = int(n.sum())
        if tot == 0:
            continue
        idx = np.repeat(np.arange(len(n)), n)
        rows = r0[idx] + (np.arange(tot) - np.repeat(np.cumsum(n) - n, n))
        t = (rows + 0.5 - y0[idx]) / (y1[idx] - y0[idx])
        xc = x0[idx] + t * (x1[idx] - x0[idx])
        cols = np.clip(np.ceil(xc - 0.5).astype(np.int64), 0, W)
        np.add.at(cnt, (rows, cols), 1)
    return (np.cumsum(cnt, axis=1)[:, :W] & 1).astype(bool)


# ----------------------------------------------------------------------------------------
# QPTIFF acquisition metadata
# ----------------------------------------------------------------------------------------

def _xml_field(xml: str, tag: str) -> str | None:
    m = re.search(rf"<{tag}>(.*?)</{tag}>", xml, flags=re.S)
    return html.unescape(m.group(1).strip()) if m else None


def parse_qpi_description(description: str | bytes | None) -> dict:
    """Channel acquisition fields from a PerkinElmer-QPI-ImageDescription page XML.

    ExposureTime is in microseconds in these files (SST page 150000 vs ScanProfile 150 ms;
    DAPI 1000 vs 1 ms), so exposure_ms = ExposureTime / 1000. Missing fields -> None.
    """
    if description is None:
        return {}
    if isinstance(description, bytes):
        for enc in ("utf-8", "utf-16", "latin-1"):
            try:
                description = description.decode(enc)
                break
            except UnicodeDecodeError:
                continue
    xml = re.sub(r"<ScanProfile>.*?</ScanProfile>", "", description, flags=re.S)
    if "PerkinElmer-QPI-ImageDescription" not in xml:
        return {}
    out = {}
    exp = _xml_field(xml, "ExposureTime")
    out["exposure_ms"] = float(exp) / 1000.0 if exp not in (None, "") else None
    out["xml_biomarker"] = _xml_field(xml, "Biomarker")
    resp = re.search(r"<Responsivity>.*?<Response>(.*?)</Response>.*?<Date>(.*?)</Date>", xml, flags=re.S)
    out["responsivity"] = float(resp.group(1)) if resp else None
    out["filter_date"] = resp.group(2).strip()[:10] if resp else None
    su = _xml_field(xml, "SignalUnits")
    out["signal_units"] = int(su) if su not in (None, "") and su.isdigit() else su
    return out


# ----------------------------------------------------------------------------------------
# SMO (Silver Mountain Operator) -- port of `smo` 2.0.2 (M. Silberberg, MIT licence,
# github.com/maurosilber/smo), the package KINTSUGI src/kintsugi/signal/smo.py wraps.
# Same filters, same mask propagation, same null (uniform image, seed 42), same rv_histogram.
# ----------------------------------------------------------------------------------------

def _masked_filter(filt, data: np.ndarray, mask: np.ndarray | None, **kw):
    """scipy.ndimage filter respecting a mask exactly as smo._filter does (mode='mirror').

    Masked pixels are filled with 0; the output mask is ~filter(~mask) evaluated on a BOOLEAN
    array (scipy truncates the weighted sum), which masks every pixel whose window touches a
    masked pixel.
    """
    if mask is not None and mask.any():
        out = filt(np.where(mask, 0.0, data), mode="mirror", **kw)
        m = ~filt(~mask, mode="mirror", **kw)
        return out, m
    return filt(data, mode="mirror", **kw), mask


def smo_operator(image: np.ndarray, mask: np.ndarray | None = None,
                 sigma: float = SMO_SIGMA, size: int = SMO_SIZE):
    """SMO image: norm of the window-averaged unit gradient of the Gaussian-smoothed image.

    Flat (background) regions have random gradient directions -> small values; structure has
    aligned directions -> values near 1. Returns (smo, output_mask or None).
    """
    from scipy.ndimage import gaussian_filter, uniform_filter
    data = np.asarray(image, dtype=np.float64)
    smoothed, m = _masked_filter(gaussian_filter, data, mask, sigma=sigma)
    grad = np.gradient(smoothed)
    norm = np.sqrt(sum(g * g for g in grad))
    ok = norm > 0
    if m is not None:
        ok &= ~m
    for g in grad:
        np.divide(g, norm, where=ok, out=g)
    outs = []
    out_mask = m
    for g in grad:
        o, out_mask = _masked_filter(uniform_filter, g, m, size=size)
        outs.append(o)
    return np.sqrt(sum(o * o for o in outs)), out_mask


def _rv_histogram(x: np.ndarray):
    """scipy rv_histogram of np.histogram(x, bins='auto'), as smo.smo._rv builds it."""
    from scipy.stats import rv_histogram
    return rv_histogram(np.histogram(np.ravel(x), bins="auto"), density=False)


_NULL_CACHE: dict = {}


def smo_null_rv(sigma: float = SMO_SIGMA, size: int = SMO_SIZE, shape=SMO_NULL_SHAPE, seed: int = 42):
    """Distribution of SMO on an iid uniform image (the SMO null hypothesis)."""
    key = (float(sigma), int(size), tuple(shape), int(seed))
    if key not in _NULL_CACHE:
        rng = np.random.default_rng(seed)
        s, _ = smo_operator(rng.uniform(size=shape), None, sigma, size)
        _NULL_CACHE[key] = _rv_histogram(s)
    return _NULL_CACHE[key]


def smo_background(image: np.ndarray, operator_mask: np.ndarray | None, sample_mask: np.ndarray,
                   sigma: float = SMO_SIGMA, size: int = SMO_SIZE,
                   threshold: float = SMO_THRESHOLD) -> dict:
    """SMO background statistics, KINTSUGI-style: SMO(sigma, size) -> bg_rv().mean()/.std().

    operator_mask: pixels hidden from the operator (out of tissue, saturated).
    sample_mask:   pixels eligible for the background sample (in tissue AND nonzero).
    """
    smo, m = smo_operator(image, operator_mask, sigma, size)
    thr = smo_null_rv(sigma, size).ppf(threshold)
    bg = smo <= thr
    if m is not None:
        bg &= ~m
    sel = bg & sample_mask
    n_eligible = int(sample_mask.sum())
    x = np.asarray(image)[sel].astype(np.float64)
    if x.size < 10:
        return dict(bg_smo_mean=np.nan, bg_smo_std=np.nan, bg_smo_median=np.nan,
                    smo_bg_fraction=x.size / max(n_eligible, 1), smo_n=int(x.size))
    rv = _rv_histogram(x)
    return dict(bg_smo_mean=float(rv.mean()), bg_smo_std=float(rv.std()),
                bg_smo_median=float(rv.median()), smo_bg_fraction=x.size / max(n_eligible, 1),
                smo_n=int(x.size))


# ----------------------------------------------------------------------------------------
# Log-space background: exact integer-aware log density, peaks, half-maximum
# ----------------------------------------------------------------------------------------

def log_density(values: np.ndarray, step: float = LOG_STEP, smooth: float = LOG_SMOOTH):
    """Smoothed density of log(intensity) for positive integer counts.

    Each integer count k carries its mass uniformly over [k-0.5, k+0.5] (what dithering
    would do, but exact and deterministic), so the density at low counts is not a comb.
    Returns (bin centres in log units, density normalised to sum 1).
    """
    from scipy.ndimage import gaussian_filter1d
    v = np.asarray(values)
    v = v[v > 0]
    if v.size == 0:
        return np.array([]), np.array([])
    if np.issubdtype(v.dtype, np.integer) or np.all(np.mod(v, 1) == 0):
        counts = np.bincount(v.astype(np.int64))
    else:  # non-integer data: plain histogram of log values (no dequantisation needed)
        vmax = float(v.max())
        edges = np.arange(np.log(v.min()), np.log(vmax) + step, step)
        h, _ = np.histogram(np.log(v), bins=edges)
        hs = gaussian_filter1d(h.astype(float), smooth / step)
        return 0.5 * (edges[1:] + edges[:-1]), hs / hs.sum()
    kmax = counts.size - 1
    cum = np.cumsum(counts).astype(np.float64)  # cum[n] = #values <= n
    edges = np.arange(np.log(0.5), np.log(kmax + 0.5) + step, step)
    y = np.exp(edges) - 0.5                      # CDF in linear space at e^edge
    n = np.clip(np.floor(y).astype(np.int64), 0, kmax)
    frac = np.clip(y - n, 0.0, 1.0)
    nxt = np.clip(n + 1, 0, kmax)
    F = cum[n] + np.where(n + 1 <= kmax, counts[nxt], 0) * frac
    F = np.minimum(F, cum[-1])
    h = np.diff(F)
    hs = gaussian_filter1d(h, smooth / step)
    return 0.5 * (edges[1:] + edges[:-1]), hs / max(hs.sum(), 1e-300)


def background_peaks(centres: np.ndarray, dens: np.ndarray,
                     min_height: float = PEAK_MIN_HEIGHT,
                     min_prominence: float = PEAK_MIN_PROMINENCE) -> dict:
    """Peaks of the log density and the tissue-background estimates derived from them.

    tissue peak = the highest-intensity SUBSTANTIAL peak (>= min_height of the tallest, with
    prominence >= min_prominence of it). The endocrine tail never forms one (islets are a few
    percent of tissue); a lower substantial peak is empty space within the Tissue outline.
    """
    from scipy.signal import find_peaks
    out = dict(logmode_global=np.nan, bg_logmode=np.nan, bg_halfmax=np.nan,
               logmode_lower=np.nan, lower_mass_fraction=0.0, n_peaks=0)
    if dens.size == 0:
        return out
    top = float(dens.max())
    pk, props = find_peaks(dens, height=min_height * top, prominence=min_prominence * top)
    ig = int(np.argmax(dens))
    if pk.size == 0:
        pk = np.array([ig])
    out["n_peaks"] = int(pk.size)
    out["logmode_global"] = float(np.exp(centres[ig]))
    it = int(pk.max())                      # upper substantial peak
    il = int(pk.min())
    out["bg_logmode"] = float(np.exp(centres[it]))
    # right half-maximum of the tissue peak
    half = 0.5 * dens[it]
    right = np.nonzero(dens[it:] <= half)[0]
    ir = it + int(right[0]) if right.size else dens.size - 1
    if right.size and ir > it:   # interpolate between ir-1 and ir
        d0, d1 = dens[ir - 1], dens[ir]
        t = (d0 - half) / (d0 - d1) if d0 != d1 else 0.0
        x_r = centres[ir - 1] + t * (centres[ir] - centres[ir - 1])
    else:
        x_r = centres[ir]
    out["bg_halfmax"] = float(np.exp(x_r))
    if il != it:
        valley = il + int(np.argmin(dens[il:it + 1]))
        out["logmode_lower"] = float(np.exp(centres[il]))
        out["lower_mass_fraction"] = float(dens[:valley + 1].sum() / dens.sum())
    return out


def log_background(values: np.ndarray) -> dict:
    """bg_logmode / bg_halfmax etc. for positive values (see background_peaks)."""
    c, d = log_density(values)
    return background_peaks(c, d)


# ----------------------------------------------------------------------------------------
# Per image
# ----------------------------------------------------------------------------------------

def _crop_to_mask(mask: np.ndarray, margin: int):
    rows = np.nonzero(mask.any(axis=1))[0]
    cols = np.nonzero(mask.any(axis=0))[0]
    if rows.size == 0:
        return slice(0, 0), slice(0, 0)
    r0 = max(int(rows[0]) - margin, 0)
    r1 = min(int(rows[-1]) + margin + 1, mask.shape[0])
    c0 = max(int(cols[0]) - margin, 0)
    c1 = min(int(cols[-1]) + margin + 1, mask.shape[1])
    return slice(r0, r1), slice(c0, c1)


def _sample_level0_zero_fraction(series, ci: int, tissue: np.ndarray, dx: float, dy: float,
                                 n_tiles: int = L0_TILES, tile: int = L0_TILE, seed: int = 0):
    """Exact-zero fraction at full resolution on random tiles lying fully inside tissue."""
    import zarr
    store = series.aszarr(level=0)
    try:
        z = zarr.open(store, mode="r")
        H0, W0 = z.shape[-2:]
        rng = np.random.default_rng(seed)
        cand = np.argwhere(tissue)
        if cand.size == 0:
            return np.nan
        zeros = total = got = tries = 0
        while got < n_tiles and tries < n_tiles * 50:
            tries += 1
            r, c = cand[rng.integers(len(cand))]
            y0 = int((r + 0.5) * dy) - tile // 2
            x0 = int((c + 0.5) * dx) - tile // 2
            if y0 < 0 or x0 < 0 or y0 + tile > H0 or x0 + tile > W0:
                continue
            rs = slice(int(y0 / dy), int(np.ceil((y0 + tile) / dy)))
            cs = slice(int(x0 / dx), int(np.ceil((x0 + tile) / dx)))
            if tissue[rs, cs].mean() < 0.999:
                continue
            a = np.asarray(z[ci, y0:y0 + tile, x0:x0 + tile])
            zeros += int((a == 0).sum())
            total += a.size
            got += 1
        return zeros / total if total else np.nan
    finally:
        try:
            store.close()
        except Exception:
            pass


def anchors_for_image(im: ProjectImage, channels=CHANNELS, manual_dir: Path = DEFAULT_MANUAL_DIR,
                      target_um: float = TARGET_UM, tissue_geojson: Path | None = None,
                      with_smo: bool = True, with_l0: bool = True) -> list[dict]:
    """One row per channel for one project image (see module docstring for the columns)."""
    import tifffile
    t0 = time.time()
    gj = tissue_geojson or (manual_dir / f"{im.tag}.geojson")
    if not gj.exists():
        raise FileNotFoundError(f"no Tissue annotation export for {im.image}: {gj}")
    rings = read_annotation_rings(gj, ("Tissue",))
    if not rings:
        raise ValueError(f"{gj} has no Tissue polygon")
    rows = []
    with tifffile.TiffFile(im.path) as tf:
        s = tf.series[0]
        full = s.levels[0].shape[-2:]
        if (full[0], full[1]) != (im.height, im.width):
            raise ValueError(f"{im.image}: TIFF {full} != project {im.height}x{im.width}")
        shapes = [lv.shape[-2:] for lv in s.levels]
        li, level_um, dx, dy = choose_level(shapes, full, im.pixel_um, target_um)
        lv = s.levels[li]
        tissue = rasterize_even_odd(rings, lv.shape[-2:], dx, dy)
        rs, cs = _crop_to_mask(tissue, margin=4 * SMO_SIZE)
        tissue_c = tissue[rs, cs]
        n_tissue = int(tissue_c.sum())
        log.info("%s: level %d (%.3f um/px, %s), tissue %.1f mm2", im.image, li, level_um,
                 "x".join(map(str, lv.shape[-2:])), n_tissue * level_um ** 2 / 1e6)
        n_pages = len(lv.pages)
        for ch in channels:
            ci = resolve_channel(im.channels, ch)
            if ci >= n_pages:
                raise ValueError(f"{im.image}: channel {ch} index {ci} >= {n_pages} pages")
            a = lv.pages[ci].asarray()[rs, cs]
            meta = parse_qpi_description(s.levels[0].pages[ci].description)
            v = a[tissue_c]
            nz = v[v > 0]
            zero_fraction = 1.0 - nz.size / max(v.size, 1)
            row = dict(image=im.image, tag=im.tag, channel=ch, level_um=round(level_um, 4),
                       n_pixels=int(v.size), zero_fraction=float(zero_fraction),
                       median_nonzero=float(np.median(nz)) if nz.size else np.nan,
                       p99_nonzero=float(np.percentile(nz, 99)) if nz.size else np.nan,
                       saturated_fraction=float((v >= SATURATION).mean()) if v.size else np.nan,
                       channel_index=ci, entry_id=im.entry_id, level_index=li,
                       image_path=str(im.path), method_version=METHOD_VERSION,
                       anchor_source=ANCHOR_SOURCE)
            row.update(log_background(nz))
            if with_smo:
                op_mask = ~tissue_c | (a >= SATURATION)
                row.update({k: v_ for k, v_ in smo_background(a, op_mask, tissue_c & (a > 0)).items()
                            if k != "smo_n"})
            else:
                row.update(bg_smo_mean=np.nan, bg_smo_std=np.nan, bg_smo_median=np.nan,
                           smo_bg_fraction=np.nan)
            row["zero_fraction_l0"] = (_sample_level0_zero_fraction(s, ci, tissue, dx, dy)
                                       if with_l0 else np.nan)
            row["exposure_ms"] = meta.get("exposure_ms")
            row["xml_biomarker"] = meta.get("xml_biomarker")
            row["responsivity"] = meta.get("responsivity")
            row["filter_date"] = meta.get("filter_date")
            row["signal_units"] = meta.get("signal_units")
            xb = normalize_channel_name(meta.get("xml_biomarker"))
            row["biomarker_match"] = (None if xb is None else
                                      bool(xb == ch or xb.lower() in [x.lower() for x in CHANNEL_ALIASES.get(ch, ())]))
            if row["biomarker_match"] is False:
                log.warning("%s: project channel %s (index %d) but QPTIFF Biomarker=%r",
                            im.image, ch, ci, meta.get("xml_biomarker"))
            row["anchor"] = row["bg_halfmax"] if ANCHOR_SOURCE == "halfmax" else row["bg_logmode"]
            if not np.isfinite(row["anchor"]) or row["anchor"] <= 0:
                log.warning("%s %s: no valid anchor", im.image, ch)
            row["computed_at"] = _dt.datetime.now().astimezone().isoformat(timespec="seconds")
            rows.append(row)
            log.info("  %-4s idx=%-2d zero=%.3f smo=%.1f/%.1f logmode=%.1f halfmax=%.1f exp=%s ms",
                     ch, ci, zero_fraction, row["bg_smo_mean"], row["bg_smo_std"], row["bg_logmode"],
                     row["bg_halfmax"], row["exposure_ms"])
    log.info("%s done in %.1f s", im.image, time.time() - t0)
    return rows


# ----------------------------------------------------------------------------------------
# CSV
# ----------------------------------------------------------------------------------------

def upsert_csv(rows: list[dict], out: Path) -> "pandas.DataFrame":
    """Insert or replace rows keyed by (image, channel); atomic write; stable sort."""
    import pandas as pd
    new = pd.DataFrame(rows)
    for c in COLUMNS:
        if c not in new.columns:
            new[c] = np.nan
    if out.exists() and out.stat().st_size > 0:
        old = pd.read_csv(out)
        keys = set(zip(new["image"], new["channel"]))
        old = old[[(i, c) not in keys for i, c in zip(old["image"], old["channel"])]]
        df = pd.concat([old, new], ignore_index=True)
    else:
        df = new
    extra = [c for c in df.columns if c not in COLUMNS]
    df = df[COLUMNS + extra].sort_values(["image", "channel"], kind="stable").reset_index(drop=True)
    out.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=out.name + ".", suffix=".tmp", dir=str(out.parent))
    os.close(fd)
    try:
        df.to_csv(tmp, index=False)
        os.replace(tmp, out)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)
    return df


def existing_images(out: Path, method_version: str = METHOD_VERSION) -> set[str]:
    import pandas as pd
    if not out.exists():
        return set()
    df = pd.read_csv(out)
    ok = df[df["method_version"].astype(str) == str(method_version)]
    counts = ok.groupby("image")["channel"].nunique()
    return set(counts[counts >= len(CHANNELS)].index)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--images", nargs="+", help="substrings of project image names")
    g.add_argument("--all", action="store_true", help="every project image")
    ap.add_argument("--project", type=Path, default=DEFAULT_PROJECT)
    ap.add_argument("--manual-dir", type=Path, default=DEFAULT_MANUAL_DIR,
                    help="folder of <tag>.geojson Tissue exports")
    ap.add_argument("--tissue-geojson", type=Path,
                    help="explicit Tissue GeoJSON for a single image (e.g. 6457, whose .qpdata is "
                         "empty: results/Geojson/Tissue/Tissue__6457_...geojson)")
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--channels", default=",".join(CHANNELS))
    ap.add_argument("--target-um", type=float, default=TARGET_UM)
    ap.add_argument("--skip-existing", action="store_true",
                    help="skip images that already have all channels at this method_version")
    ap.add_argument("--no-smo", action="store_true", help="skip the SMO cross-check (faster)")
    ap.add_argument("--no-l0", action="store_true", help="skip level-0 zero-fraction sampling")
    ap.add_argument("-v", "--verbose", action="store_true")
    a = ap.parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if a.verbose else logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s", datefmt="%H:%M:%S")
    channels = tuple(c.strip() for c in a.channels.split(",") if c.strip())
    bad = [c for c in channels if c in FORBIDDEN_CHANNELS]
    if bad:
        raise SystemExit(f"refusing forbidden channels: {bad}")
    images = select_images(read_project(a.project), None if a.all else a.images)
    if a.tissue_geojson and len(images) != 1:
        raise SystemExit("--tissue-geojson needs exactly one selected image")
    done = existing_images(a.out) if a.skip_existing else set()
    n_ok = n_fail = 0
    for im in images:
        if im.image in done:
            log.info("skip %s (already at method_version %s)", im.image, METHOD_VERSION)
            continue
        try:
            rows = anchors_for_image(im, channels, a.manual_dir, a.target_um,
                                     tissue_geojson=a.tissue_geojson,
                                     with_smo=not a.no_smo, with_l0=not a.no_l0)
        except Exception as exc:  # keep going; report at the end
            log.error("FAILED %s: %s", im.image, exc)
            n_fail += 1
            continue
        upsert_csv(rows, a.out)   # after every image, so an interrupted --all keeps its work
        n_ok += 1
    log.info("anchors: %d image(s) written to %s, %d failed", n_ok, a.out, n_fail)
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
