#!/usr/bin/env python3
"""islet_anchors.py -- per-image, per-channel background anchors and clamp models for the v3 islet detector.

The v3 detector (qupath/1_panc_objects_v3.groovy, ISLET_METHOD='drotar_cd99') calls a pixel positive
for channel c when its smoothed value v satisfies

    v >= max(tau * anchor_c,  clamp_mu_c + k * clamp_sigma_c)          (k <= 0: the floor is off)

This module estimates, for CD99, INS, GCG and SST, the two per-image quantities in that rule and
writes them to results/v3/islet_anchors.csv, one row per (image, channel).

Governing rule (Phenocycler_Analysis docs/INTENSITY_SCALE.md): *standardize the unit, never the
distribution*. Every quantity is a property of the image's own tissue BACKGROUND (endocrine signal
is excluded, see "two passes"), estimated in the same pixel space the detector thresholds, and applied
as a linear scale or an additive floor. Nothing matches distributions between images, and nothing is
a quantile of the positive (endocrine) population.

Pixel space (v2.0): the detector's. The pyramid level nearest --resolution-um (1.0156 um = level 1)
is read and blurred with the detector's Gaussian (--sigma-um 1.5 um, truncate 4 sigma), so the anchor
and the detector see identical values. (v1.0 estimated at ~2 um without blur, so the effective tau
for INS differed by up to 21% between images.) Statistics use a seeded random subsample of at most
--sample-max mask pixels; the blur itself runs on the full crop.

Analysis mask (v2.0): union of the Tissue polygons MINUS the union of the excluded classes (pLN,
Artifact), clipped EDGE_MARGIN_UM inside the image edge -- the same area the detector analyses.
Polygons are rasterised one at a time and OR-ed (v1.0 filled all rings even-odd, so overlapping
Tissue polygons cancelled). A SHA-256 over the polygons used (annotation_hash: class and every vertex
at the export's 2 decimals, see annotation_signature) is stored with the anchors. The detector
recomputes it from its hierarchy and refuses to run (ANCHOR_CHECK='error', default) or warns when they
differ; --skip-existing recomputes an image whose hash or estimation parameters changed.

No DAPI mask (v2.1). Until 2026-10-01 a pixel counted as "occupied tissue" only where blurred DAPI
exceeded the valley between an empty-space and a tissue peak. The maintainer dropped it: the markers
cover far more than nucleated cells, and the redrawn Tissue annotations already separate the block
from the glass. Anchors and clamp models now use every analysis-mask pixel that carries data; only
"no data" pixels (all four detection channels exactly 0) are left out. On an image whose outline
encloses much empty space the level-0 clamp fit is therefore lower than v2.0's (6450panc, v2.0
outline: CD99 mu 5 vs 93 counts with the mask).

Two passes (v2.0): pass 1 estimates every channel on all data-carrying analysis pixels; the detector rule at
(EXCL_TAU_CD99 3, EXCL_TAU_H 5, EXCL_FLOOR_K 3) then marks endocrine pixels, which are dilated by
EXCL_DILATE_UM (25 um); pass 2 re-estimates on those pixels outside that mask. This removes the
dependence of the anchor on islet prevalence and on the signal halo around islets (24 images: pass 2 is 3%
lower at the median, >= 12% lower for a tenth of the (image, channel) pairs, up to 26%: 6534 SST).
anchor_pass1 and excluded_fraction (12-32% of the mask) are recorded.

Columns per channel (see COLUMNS):
  anchor, bg_halfmax      right half-maximum of the tissue-background peak of the log density of the
                          blurred values (the UPPER substantial peak when bimodal); scale-equivariant.
  bg_logmode              that peak's mode.
  clamp_fraction          fraction of exact zeros among background pixels at FULL resolution
                          (level 0), where the scanner clamps negative values after subtracting offset
                          and autofluorescence (INS: up to ~0.8 on T1D/AAB images).
  clamp_mu, clamp_sigma   the latent background N(mu, sigma^2) in LINEAR counts, fitted to level-0
                          integer values by maximum likelihood with zeros left-censored (value < 0.5)
                          and the likelihood truncated above the median of the positives (so the fit
                          uses the zero mass plus the lower half of the positives and is not driven by
                          signal or speckle in the upper tail). A linear model, not a log model: the
                          level-0 histograms are flat across 1-8 counts next to the zero spike (e.g.
                          6533 INS ~0.8%/count after 80% zeros), i.e. a normal whose mean lies BELOW
                          zero, which no log-normal can produce (a zero-censored log-normal fitted the
                          same way runs to its bounds with KS 0.035-0.22, vs <= 0.018 for this model).
                          Exactly scale-equivariant (clamp at 0 is invariant).
                          clamp_mu + k*clamp_sigma is the detector's positivity floor: it keeps the
                          threshold from collapsing when over-subtraction shrinks the half-max anchor
                          (Spearman(zero fraction, v1 INS anchor) = -0.92 over the 24 images).
  clamp_fit_ks            max |ECDF - model CDF| over [0, upper] (fit diagnostic).
  zero_fraction           exact zeros among mask pixels at the detector scale (after the blur; rare).
  zero_fraction_l0        exact zeros among all sampled level-0 mask pixels (pass 1).
  anchor_l2_unblurred     the v1.0 estimator (half-max at ~2 um, no blur) on the v2.0 mask, for comparison.
  bg_smo_mean/std/median  SMO (KINTSUGI signal/smo.py port) at ~2 um, recorded as a diagnostic only
                          (it is not scale-robust at this resolution; see ISLETS.md). --no-smo skips it.
  exposure_ms, responsivity, filter_date, signal_units, xml_biomarker, biomarker_match
                          from the level-0 QPTIFF page XML (ExposureTime is in microseconds).

Usage:
  islet_anchors.py --images 6450panc 6374_Scan1         # substring match on image names
  islet_anchors.py --all [--skip-existing]
  islet_anchors.py --all --out /scratch/.../anchors.csv # anywhere; default results/v3/islet_anchors.csv

Idempotent: rows are upserted by (image, channel); rerunning an image replaces its rows, and the
values are deterministic (only computed_at changes).
"""
from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import html
import json
import logging
import math
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
HORMONES = ("INS", "GCG", "SST")
# Never used for detection: IAPP failed QC; HLA-A is disease-associated (rises in T1D exocrine).
FORBIDDEN_CHANNELS = frozenset({"IAPP", "HLA-A"})
# Same alias map as 1_panc_objects*.groovy, so a channel resolves to the same index in both.
CHANNEL_ALIASES = {
    "SST": ("SST", "Somatostatin", "Stt"),
    "GCG": ("GCG", "Glucagon", "Gcg", "Opal 520"),
    "INS": ("INS", "Insulin", "Ins", "Opal 690"),
    "CD99": ("CD99",),
    "DAPI": ("DAPI",),
}

METHOD_VERSION = "2.1"          # 2.1 (2026-10-01): the DAPI occupied-tissue mask is gone
RESOLUTION_UM = 1.0156        # the detector's RESOLUTION_UM (pyramid level 1)
SIGMA_UM = 1.5                # the detector's SIGMA_UM
GAUSS_TRUNCATE = 4.0          # kernel radius in sigmas (QuPath/OpenCV use a similar 4-sigma kernel)
TISSUE_CLASSES = ("Tissue",)
EXCLUDE_CLASSES = ("pLN", "Artifact")   # the detector's EXCLUDE_CLASSES
EDGE_MARGIN_UM = 5.0          # the detector clips 5 um inside the image edge
HASH_DTYPES = {"annotation_hash": str, "params_hash": str, "method_version": str}   # read back as text, never numbers

# two-pass exclusion: the detector's provisional tau values; k = 3 is one step below the detector's FLOOR_K (4),
# so the exclusion is the more inclusive of the two (it also drops weaker speckle from the background sample)
EXCL_TAU_CD99 = 3.0
EXCL_TAU_H = 5.0
EXCL_FLOOR_K = 3.0
EXCL_DILATE_UM = 25.0
EXCL_BLOCK = 4                # dilation is computed on 4x4 blocks (+-4 um), then expanded

SAMPLE_MAX = 4_000_000        # mask pixels per channel for the log density
SEED = 20260928

SMO_SIGMA = 1.0
SMO_SIZE = 7
SMO_THRESHOLD = 0.05          # probability under the SMO null (smo.SMO.bg_rv default)
SMO_NULL_SHAPE = (1024, 1024)  # null depends only on sigma/size; 1M samples fix its 5% quantile
SMO_TARGET_UM = 2.0

LOG_STEP = 0.01               # log-density grid (natural-log units)
LOG_SMOOTH = 0.10             # Gaussian smoothing of the log-density (log units)
PEAK_MIN_HEIGHT = 0.25        # a "substantial" peak is >= 25% of the tallest one ...
PEAK_MIN_PROMINENCE = 0.05    # ... and stands out by >= 5% of it
ANCHOR_SOURCE = "halfmax"     # which estimator becomes the anchor column

L0_TILES = 64                 # level-0 tiles sampled for the clamp model
L0_TILE = 512
L0_TILE_MIN_MASK = 0.9        # tile footprint inside the analysis mask
CLAMP_UPPER_Q = 0.5           # right-censor the positives above their median
CLAMP_MIN_N = 5000
CLAMP_MAX_FRACTION = 0.95     # above this share of zeros the latent normal is not identifiable
SATURATION = 65535

COLUMNS = [
    # used by the detector
    "image", "tag", "channel", "anchor", "clamp_mu", "clamp_sigma", "clamp_fraction",
    "annotation_hash", "annotation_classes", "params_hash", "method_version", "computed_at",
    # estimation context
    "level_um", "level_index", "sigma_um", "n_pixels", "n_occupied", "n_sampled", "n_sampled_pass2",
    "occupied_fraction", "nodata_fraction", "excluded_fraction", "zero_fraction",
    "anchor_source", "anchor_pass1", "bg_halfmax", "bg_logmode", "logmode_global", "logmode_lower",
    "lower_mass_fraction", "n_peaks", "median", "p99",
    "clamp_n", "clamp_upper", "clamp_fit_ks", "clamp_converged", "clamp_mu_pass1", "clamp_sigma_pass1",
    "clamp_fraction_pass1", "zero_fraction_l0",
    "anchor_l2_unblurred", "bg_smo_mean", "bg_smo_std", "bg_smo_median", "smo_bg_fraction",
    "excl_tau_cd99", "excl_tau_h", "excl_floor_k", "excl_dilate_um",
    "exposure_ms", "responsivity", "filter_date", "signal_units", "xml_biomarker", "biomarker_match",
    "saturated_fraction", "channel_index", "entry_id", "image_path", "annotation_source",
    "n_tissue_polygons", "n_excluded_polygons", "excluded_area_mm2", "analysis_area_mm2",
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
                 pixel_um: float, target_um: float = RESOLUTION_UM) -> tuple[int, float, float, float]:
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
# Annotations: polygons, rasterisation (union of polygons), hash
# ----------------------------------------------------------------------------------------

def read_annotation_polygons(geojson: Path, classes) -> list[tuple[str, list[np.ndarray]]]:
    """(class, [exterior, hole, ...]) for every polygon of every feature whose class is in `classes`.

    Coordinates are QuPath full-resolution pixels, exactly as exported (QuPath writes 2 decimals).
    A MultiPolygon contributes one entry per part. Rings are kept as exported (the annotation
    signature needs them all); rasterisation and areas skip rings with fewer than 3 vertices.
    """
    with open(geojson, encoding="utf-8") as fh:
        gj = json.load(fh)
    feats = gj["features"] if gj.get("type") == "FeatureCollection" else [gj]
    out = []
    classes = tuple(classes)
    for ft in feats:
        cls = ((ft.get("properties") or {}).get("classification") or {}).get("name")
        if cls not in classes:
            continue
        g = ft.get("geometry") or {}
        if g.get("type") == "Polygon":
            polys = [g["coordinates"]]
        elif g.get("type") == "MultiPolygon":
            polys = g["coordinates"]
        else:
            continue
        for poly in polys:
            rings = [np.asarray(r, dtype=np.float64).reshape(-1, 2) if len(r) else np.zeros((0, 2)) for r in poly]
            if rings and len(rings[0]) >= 3:
                out.append((cls, rings))
    return out


def read_annotation_rings(geojson: Path, classes=("Tissue",)) -> list[np.ndarray]:
    """All rings (exteriors and holes) of the matching features (v1.0 interface, kept for scripts)."""
    return [r for _, rings in read_annotation_polygons(geojson, classes) for r in rings]


def rasterize_even_odd(rings: list[np.ndarray], shape: tuple[int, int],
                       dx: float = 1.0, dy: float = 1.0, origin: tuple[int, int] = (0, 0)) -> np.ndarray:
    """Boolean mask of `shape`; pixel (r, c) is inside if its centre ((c+c0+.5)dx, (r+r0+.5)dy) is.

    Even-odd over the given rings: exact for ONE valid polygon (exterior + holes). `rings` are
    closed coordinate arrays in full-resolution pixels; dx, dy are the level downsamples; origin
    (r0, c0) is the level pixel at the window's top-left. Scanline, O(perimeter + pixels).
    """
    H, W = shape
    r0w, c0w = origin
    cnt = np.zeros((H, W + 1), dtype=np.int32)
    for ring in rings:
        if len(ring) < 3:
            continue
        x = ring[:, 0] / dx - c0w
        y = ring[:, 1] / dy - r0w
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


def rasterize_union(polys: list[tuple[str, list[np.ndarray]]], shape: tuple[int, int],
                    dx: float = 1.0, dy: float = 1.0) -> np.ndarray:
    """Union of polygons: each polygon is filled even-odd (exterior minus holes) in its own
    bounding box and OR-ed into the result, so overlapping polygons never cancel."""
    H, W = shape
    out = np.zeros((H, W), dtype=bool)
    for _, rings in polys:
        allc = np.concatenate(rings)
        c0 = max(int(np.floor(allc[:, 0].min() / dx)) - 1, 0)
        c1 = min(int(np.ceil(allc[:, 0].max() / dx)) + 1, W)
        r0 = max(int(np.floor(allc[:, 1].min() / dy)) - 1, 0)
        r1 = min(int(np.ceil(allc[:, 1].max() / dy)) + 1, H)
        if r1 <= r0 or c1 <= c0:
            continue
        out[r0:r1, c0:c1] |= rasterize_even_odd(rings, (r1 - r0, c1 - c0), dx, dy, origin=(r0, c0))
    return out


def analysis_mask(tissue: list, excluded: list, shape: tuple[int, int], dx: float, dy: float,
                  full_shape: tuple[int, int] | None = None, edge_margin_px0: float = 0.0) -> np.ndarray:
    """union(Tissue) - union(excluded classes), optionally clipped edge_margin_px0 full-resolution
    pixels inside the image edge (pixel centres, as the detector's inner rectangle)."""
    m = rasterize_union(tissue, shape, dx, dy)
    if excluded:
        m &= ~rasterize_union(excluded, shape, dx, dy)
    if full_shape is not None and edge_margin_px0 > 0:
        H0, W0 = full_shape
        yc = (np.arange(shape[0]) + 0.5) * dy
        xc = (np.arange(shape[1]) + 0.5) * dx
        m &= ((yc >= edge_margin_px0) & (yc <= H0 - edge_margin_px0))[:, None]
        m &= ((xc >= edge_margin_px0) & (xc <= W0 - edge_margin_px0))[None, :]
    return m


def polygons_area_um2(polys: list, pixel_um: float) -> float:
    """Area of the union of the polygons in um^2 (overlaps counted once)."""
    if not polys:
        return 0.0
    from shapely.geometry import Polygon
    from shapely.ops import unary_union
    u = unary_union([Polygon(r[0], [h for h in r[1:] if len(h) >= 3]).buffer(0) for _, r in polys])
    return float(u.area) * pixel_um ** 2


def annotation_signature(polys: list[tuple[str, list[np.ndarray]]]) -> list[str]:
    """One canonical line per polygon, 'class|x,y;x,y;.../x,y;...' (rings separated by '/'), every
    vertex printed with 2 decimals -- the precision of QuPath's GeoJSON export. Lines are sorted, so
    feature order does not matter.

    The detector (1_panc_objects_v3.groovy) serialises its current Tissue/pLN/Artifact annotations with
    the same exporter (PathIO.exportObjectsAsGeoJSON), parses that text and builds the same lines. Both
    sides therefore format identical doubles and no floating-point arithmetic is involved (areas or
    centroids computed by GEOS and JTS differ in the last bits, which rounding does not always hide)."""
    lines = []
    for cls, rings in polys:
        body = "/".join(";".join(f"{x:.2f},{y:.2f}" for x, y in np.asarray(r, dtype=np.float64).reshape(-1, 2))
                        for r in rings)
        lines.append(f"{cls}|{body}")
    return sorted(lines)


def annotation_hash(polys: list[tuple[str, list[np.ndarray]]], classes) -> str:
    """First 16 hex digits of SHA-256 over 'classes=<sorted classes>' and the sorted signature lines."""
    text = "classes=" + ",".join(sorted(classes)) + "\n" + "\n".join(annotation_signature(polys))
    return hashlib.sha256(text.encode("utf-8")).hexdigest()[:16]


def estimation_params(resolution_um: float = RESOLUTION_UM, sigma_um: float = SIGMA_UM,
                      exclude_classes=EXCLUDE_CLASSES, two_pass: bool = True,
                      excl_tau_cd99: float = EXCL_TAU_CD99, excl_tau_h: float = EXCL_TAU_H,
                      excl_floor_k: float = EXCL_FLOOR_K, excl_dilate_um: float = EXCL_DILATE_UM,
                      sample_max: int = SAMPLE_MAX, n_l0_tiles: int = L0_TILES, l0_tile: int = L0_TILE) -> dict:
    """Everything that changes the estimates (hashed into params_hash)."""
    return dict(method_version=METHOD_VERSION, resolution_um=float(resolution_um), sigma_um=float(sigma_um),
                exclude_classes=list(exclude_classes), two_pass=bool(two_pass),
                excl_tau_cd99=float(excl_tau_cd99), excl_tau_h=float(excl_tau_h),
                excl_floor_k=float(excl_floor_k), excl_dilate_um=float(excl_dilate_um),
                sample_max=int(sample_max), n_l0_tiles=int(n_l0_tiles), edge_margin_um=EDGE_MARGIN_UM,
                clamp_upper_q=CLAMP_UPPER_Q, l0_tile=int(l0_tile), seed=SEED)


def params_hash(params: dict) -> str:
    """Hash of the estimation parameters (so --skip-existing recomputes when they change)."""
    text = json.dumps({k: params[k] for k in sorted(params)}, sort_keys=True, default=str)
    return hashlib.sha256(text.encode("utf-8")).hexdigest()[:12]


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
# Diagnostic only.
# ----------------------------------------------------------------------------------------

def _masked_filter(filt, data: np.ndarray, mask: np.ndarray | None, **kw):
    """scipy.ndimage filter respecting a mask exactly as smo._filter does (mode='mirror')."""
    if mask is not None and mask.any():
        out = filt(np.where(mask, 0.0, data), mode="mirror", **kw)
        m = ~filt(~mask, mode="mirror", **kw)
        return out, m
    return filt(data, mode="mirror", **kw), mask


def smo_operator(image: np.ndarray, mask: np.ndarray | None = None,
                 sigma: float = SMO_SIGMA, size: int = SMO_SIZE):
    """SMO image: norm of the window-averaged unit gradient of the Gaussian-smoothed image."""
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
    """SMO background statistics, KINTSUGI-style: SMO(sigma, size) -> bg_rv().mean()/.std()."""
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
    """Smoothed density of log(intensity) for positive values.

    Integer counts: each count k carries its mass uniformly over [k-0.5, k+0.5] (exact, deterministic
    dequantisation), so the density at low counts is not a comb. Non-integer values (blurred data):
    a plain histogram of log values. Returns (bin centres in log units, density normalised to sum 1).
    """
    from scipy.ndimage import gaussian_filter1d
    v = np.asarray(values)
    v = v[v > 0]
    if v.size == 0:
        return np.array([]), np.array([])
    if np.issubdtype(v.dtype, np.integer) or np.all(np.mod(v, 1) == 0):
        counts = np.bincount(v.astype(np.int64))
    else:
        lv = np.log(v.astype(np.float64))
        lo = max(float(lv.min()), float(np.log(1e-3)))   # sub-1e-3 values: numerical dust after the blur
        lv = lv[lv >= lo]
        edges = np.arange(lo, float(lv.max()) + 2 * step, step)
        h, _ = np.histogram(lv, bins=edges)
        hs = gaussian_filter1d(h.astype(float), smooth / step)
        return 0.5 * (edges[1:] + edges[:-1]), hs / max(hs.sum(), 1e-300)
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
               logmode_lower=np.nan, lower_mass_fraction=0.0, n_peaks=0, valley=np.nan)
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
        out["valley"] = float(np.exp(centres[valley]))
    return out


def log_background(values: np.ndarray) -> dict:
    """bg_logmode / bg_halfmax etc. for positive values (see background_peaks)."""
    c, d = log_density(values)
    return background_peaks(c, d)


# ----------------------------------------------------------------------------------------
# Clamp model: latent N(mu, sigma^2) in linear counts, zeros left-censored
# ----------------------------------------------------------------------------------------

def _log_prob_between(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """log(Phi(b) - Phi(a)) for a < b, accurate in both tails."""
    from scipy.special import ndtr
    upper = a > 0
    p = np.where(upper, ndtr(-a) - ndtr(-b), ndtr(b) - ndtr(a))
    return np.log(np.maximum(p, 1e-300))


def _fit_truncated(n0: int, ck: np.ndarray, U: int, n: int):
    """MLE of (mu, sigma) from the zero count n0 and bin counts ck[k-1] (k = 1..U), conditioned on
    X <= U (right truncation: values above U carry no information, so a bright tail cannot bias it)."""
    from scipy.optimize import minimize
    from scipy.special import log_ndtr, ndtri
    ks = np.arange(1, U + 1, dtype=np.float64)
    m = n0 + float(ck.sum())
    eps = 1.0 / (n + 1)
    F0 = min(max(n0 / n, eps), 1 - eps)
    FU = min(max(m / n, F0 + eps), 1 - eps)
    z0, zU = float(ndtri(F0)), float(ndtri(FU))
    S = float(max(U, 10))
    s0 = min(max(U / max(zU - z0, 1e-3), 0.5), 40 * S)
    m0 = min(max(0.5 - z0 * s0, -40 * S), 40 * S)

    def nll(theta):
        mu, ls = theta
        s = math.exp(ls)
        ll = n0 * log_ndtr((0.5 - mu) / s)
        ll += float((ck * _log_prob_between((ks - 0.5 - mu) / s, (ks + 0.5 - mu) / s)).sum())
        ll -= m * log_ndtr((U + 0.5 - mu) / s)
        v = -ll / m
        return v if np.isfinite(v) else 1e30

    res = minimize(nll, np.array([m0, math.log(s0)]), method="Nelder-Mead",
                   bounds=[(-50 * S, 50 * S), (math.log(0.3), math.log(50 * S))],
                   options=dict(xatol=1e-5, fatol=1e-12, maxiter=4000))
    return float(res.x[0]), float(math.exp(res.x[1])), bool(res.success)


def fit_clamped_normal(values: np.ndarray, upper_q: float = CLAMP_UPPER_Q, min_n: int = CLAMP_MIN_N) -> dict:
    """Latent background Y ~ N(mu, sigma^2) behind integer counts X = max(0, round(Y)).

    Maximum likelihood on the binned counts, conditioned on X <= U: P(X = 0) = Phi((0.5 - mu)/sigma)
    (the clamped zeros are left-censored), P(X = k) = Phi((k + .5 - mu)/sigma) - Phi((k - .5 - mu)/sigma)
    for 1 <= k <= U, all divided by Phi((U + .5 - mu)/sigma). U is the `upper_q` quantile of the
    positive values, so the fit sees the zero mass and the lower part of the positives only; values
    above U (signal, speckle) do not enter at all (right truncation). A fit that hits a bound, does
    not converge, or leaves |ECDF - model| > 0.02 on [0, U] is retried with U = max(8, 10% quantile);
    if that fails too, or more than CLAMP_MAX_FRACTION of the values are zero, mu and sigma are NaN
    (the detector then uses tau * anchor alone for that channel and warns).
    Returns clamp_mu, clamp_sigma (counts), clamp_fraction (share of zeros), clamp_n, clamp_upper,
    clamp_fit_ks, clamp_converged. Scale-equivariant: fit(a*X) = a*fit(X) up to quantisation.
    """
    from scipy.special import ndtr
    v = np.asarray(values).ravel()
    out = dict(clamp_mu=np.nan, clamp_sigma=np.nan, clamp_fraction=np.nan, clamp_n=int(v.size),
               clamp_upper=np.nan, clamp_fit_ks=np.nan, clamp_converged=False)
    if v.size < min_n:
        return out
    vi = np.rint(v).astype(np.int64)
    n = vi.size
    n0 = int((vi <= 0).sum())
    out["clamp_fraction"] = n0 / n
    pos = vi[vi > 0]
    if pos.size < 50 or n0 / n > CLAMP_MAX_FRACTION:   # unidentifiable: almost all mass censored
        return out
    for U in (int(max(2, np.quantile(pos, upper_q))), int(max(8, np.quantile(pos, 0.1)))):
        cnt = np.bincount(np.minimum(pos, U + 1), minlength=U + 2)
        ck = cnt[1:U + 1].astype(np.float64)
        if ck.sum() < 50:
            continue
        mu, s, ok = _fit_truncated(n0, ck, U, n)
        S = float(max(U, 10))
        grid = np.arange(0, U + 1, dtype=np.float64)
        m = n0 + ck.sum()
        ecdf = (n0 + np.concatenate([[0.0], np.cumsum(ck)])) / m
        model = ndtr((grid + 0.5 - mu) / s) / max(ndtr((U + 0.5 - mu) / s), 1e-300)
        ks = float(np.max(np.abs(ecdf - model)))
        sane = ok and 0.31 < s < 49 * S and -49 * S < mu < 49 * S and ks <= 0.02
        out.update(clamp_upper=U, clamp_fit_ks=ks, clamp_converged=bool(sane))
        if sane:
            out.update(clamp_mu=mu, clamp_sigma=s)
            return out
    return out


def no_clamp_model() -> dict:
    """The clamp columns when no level-0 values were sampled (--no-l0): the floor is off."""
    return dict(clamp_mu=np.nan, clamp_sigma=np.nan, clamp_fraction=np.nan, clamp_n=0, clamp_upper=np.nan,
                clamp_fit_ks=np.nan, clamp_converged=False)


def positivity_threshold(anchor: float, tau: float, mu: float = np.nan, sigma: float = np.nan,
                         k: float = 0.0) -> float:
    """The detector's per-channel rule: value >= max(tau * anchor, mu + k * sigma); k <= 0 or a
    missing clamp model turns the floor off. Mirrors 1_panc_objects_v3.groovy."""
    t = tau * anchor
    if k > 0 and np.isfinite(mu) and np.isfinite(sigma):
        t = max(t, mu + k * sigma)
    return float(t)


# ----------------------------------------------------------------------------------------
# Masks: dilation in um, sampling
# ----------------------------------------------------------------------------------------

def dilate_um(mask: np.ndarray, radius_um: float, px_um: float, block: int = EXCL_BLOCK) -> np.ndarray:
    """mask dilated by radius_um (Euclidean), computed on block x block cells: conservative by at
    most block * px_um (every pixel within radius_um of a mask pixel is covered)."""
    from scipy.ndimage import distance_transform_edt
    if radius_um <= 0 or not mask.any():
        return mask.copy()
    H, W = mask.shape
    Hb, Wb = -(-H // block), -(-W // block)
    pad = np.zeros((Hb * block, Wb * block), dtype=bool)
    pad[:H, :W] = mask
    coarse = pad.reshape(Hb, block, Wb, block).any(axis=(1, 3))
    d = distance_transform_edt(~coarse) * (block * px_um)
    dil = d <= radius_um + block * px_um
    return np.repeat(np.repeat(dil, block, axis=0), block, axis=1)[:H, :W]


def sample_mask_indices(mask: np.ndarray, n_max: int, seed: int = SEED) -> tuple[np.ndarray, np.ndarray]:
    """Up to n_max distinct (row, col) positions inside `mask`, uniformly at random (seeded).

    Rejection sampling on the bounding array, so no index array of the whole mask is built."""
    rng = np.random.default_rng(seed)
    H, W = mask.shape
    n_in = int(mask.sum())
    if n_in == 0:
        return np.zeros(0, np.int64), np.zeros(0, np.int64)
    if n_in <= n_max:
        r, c = np.nonzero(mask)
        return r.astype(np.int64), c.astype(np.int64)
    frac = n_in / mask.size
    fm = np.ascontiguousarray(mask).ravel()
    got = np.zeros(0, np.int64)
    while got.size < n_max:
        m = int((n_max - got.size) / frac * 1.3) + 1024
        flat = rng.integers(0, mask.size, size=m)
        flat = flat[fm[flat]]
        got = np.unique(np.concatenate([got, flat]))
    got = rng.permutation(got)[:n_max]
    got.sort()
    return got // W, got % W


# ----------------------------------------------------------------------------------------
# Core estimation on arrays (testable without a TIFF)
# ----------------------------------------------------------------------------------------

@dataclass
class L0Tile:
    """Full-resolution values of one sampled tile and the level pixel of each tile row/column."""
    values: dict           # channel -> 2D integer array (level 0)
    rows: np.ndarray       # level row index of each tile row
    cols: np.ndarray       # level column index of each tile column


def endocrine_mask(blurred: dict, stats: dict, tau_cd99: float, tau_h: float, floor_k: float,
                   region: np.ndarray | None = None) -> np.ndarray:
    """The detector's union rule: CD99 >= thr(CD99, tau_cd99) or any hormone >= thr(h, tau_h)."""
    out = None
    for ch, arr in blurred.items():
        st = stats[ch]
        tau = tau_cd99 if ch == "CD99" else tau_h
        thr = positivity_threshold(st["anchor"], tau, st.get("clamp_mu", np.nan),
                                   st.get("clamp_sigma", np.nan), floor_k)
        if not np.isfinite(thr):
            continue
        p = arr >= thr
        out = p if out is None else (out | p)
    if out is None:
        out = np.zeros(next(iter(blurred.values())).shape, dtype=bool)
    if region is not None:
        out &= region
    return out


def _channel_stats(values: np.ndarray) -> dict:
    v = np.asarray(values, dtype=np.float64)
    bp = log_background(v)
    return dict(anchor=bp["bg_halfmax"], bg_halfmax=bp["bg_halfmax"], bg_logmode=bp["bg_logmode"],
                logmode_global=bp["logmode_global"], logmode_lower=bp["logmode_lower"],
                lower_mass_fraction=bp["lower_mass_fraction"], n_peaks=bp["n_peaks"],
                median=float(np.median(v)) if v.size else np.nan,
                p99=float(np.percentile(v, 99)) if v.size else np.nan,
                zero_fraction=float(np.mean(v == 0)) if v.size else np.nan)


def _l0_values(tiles: list[L0Tile], ch: str, level_mask: np.ndarray) -> np.ndarray:
    parts = []
    for t in tiles:
        sel = level_mask[np.ix_(t.rows, t.cols)]
        parts.append(np.asarray(t.values[ch])[sel])
    return np.concatenate(parts) if parts else np.zeros(0, np.int64)


def compute_anchors(blurred: dict, mask: np.ndarray, occupied: np.ndarray, px_um: float,
                    l0_tiles: list[L0Tile] | None = None, two_pass: bool = True,
                    excl_tau_cd99: float = EXCL_TAU_CD99, excl_tau_h: float = EXCL_TAU_H,
                    excl_floor_k: float = EXCL_FLOOR_K, excl_dilate_um: float = EXCL_DILATE_UM,
                    sample_max: int = SAMPLE_MAX, seed: int = SEED) -> dict:
    """Two-pass anchors and clamp models from blurred level arrays (channel -> float array).

    mask: analysis mask; occupied: the pixels that carry data (same shape; AND-ed with mask).
    Returns {channel: stats} plus '_meta' (fractions and counts). Pass 1 uses every occupied
    pixel; pass 2 excludes the endocrine mask of pass 1 dilated by excl_dilate_um.
    """
    l0_tiles = l0_tiles or []
    occ = mask & occupied
    rr, cc = sample_mask_indices(occ, sample_max, seed)
    stats1 = {}
    for ch, arr in blurred.items():
        st = _channel_stats(arr[rr, cc])
        st.update(fit_clamped_normal(_l0_values(l0_tiles, ch, occ)) if l0_tiles else no_clamp_model())
        stats1[ch] = st
    meta = dict(n_pixels=int(mask.sum()), n_occupied=int(occ.sum()), n_sampled=int(rr.size),
                excluded_fraction=0.0, n_sampled_pass2=int(rr.size))
    if not two_pass:
        for ch, st in stats1.items():
            st.update(anchor_pass1=st["anchor"], clamp_mu_pass1=st["clamp_mu"],
                      clamp_sigma_pass1=st["clamp_sigma"], clamp_fraction_pass1=st["clamp_fraction"])
        stats1["_meta"] = meta
        return stats1
    endo = endocrine_mask(blurred, stats1, excl_tau_cd99, excl_tau_h, excl_floor_k, region=mask)
    excl = dilate_um(endo, excl_dilate_um, px_um)
    keep = ~excl[rr, cc]
    occ2 = occ & ~excl
    meta.update(excluded_fraction=float(1.0 - occ2.sum() / max(occ.sum(), 1)),
                n_sampled_pass2=int(keep.sum()), endocrine_fraction=float(endo.sum() / max(mask.sum(), 1)))
    out = {}
    for ch, arr in blurred.items():
        st = _channel_stats(arr[rr[keep], cc[keep]])
        st.update(fit_clamped_normal(_l0_values(l0_tiles, ch, occ2)) if l0_tiles else no_clamp_model())
        s1 = stats1[ch]
        st.update(anchor_pass1=s1["anchor"], clamp_mu_pass1=s1["clamp_mu"],
                  clamp_sigma_pass1=s1["clamp_sigma"], clamp_fraction_pass1=s1["clamp_fraction"])
        out[ch] = st
    out["_meta"] = meta
    return out


# ----------------------------------------------------------------------------------------
# Per image (I/O)
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


def gaussian_blur(a: np.ndarray, sigma_px: float) -> np.ndarray:
    from scipy.ndimage import gaussian_filter
    return gaussian_filter(np.asarray(a, dtype=np.float32), sigma_px, mode="reflect", truncate=GAUSS_TRUNCATE)


def _sample_l0_tiles(series, channel_idx: dict, occ: np.ndarray, rs: slice, cs: slice,
                     dx: float, dy: float, n_tiles: int = L0_TILES, tile: int = L0_TILE,
                     seed: int = SEED) -> tuple[list[L0Tile], np.ndarray]:
    """Random level-0 tiles whose level footprint lies >= L0_TILE_MIN_MASK inside `occ` (crop frame).

    Returns the tiles (rows/cols in CROP level coordinates) and an int array of the tile origins."""
    import zarr
    store = series.aszarr(level=0)
    try:
        z = zarr.open(store, mode="r")
        H0, W0 = z.shape[-2:]
        rng = np.random.default_rng(seed + 1)
        rr, cc = sample_mask_indices(occ, 200_000, seed + 2)
        tiles, origins = [], []
        tries = 0
        while len(tiles) < n_tiles and tries < n_tiles * 60 and rr.size:
            tries += 1
            i = int(rng.integers(rr.size))
            r, c = int(rr[i]) + rs.start, int(cc[i]) + cs.start
            y0 = int((r + 0.5) * dy) - tile // 2
            x0 = int((c + 0.5) * dx) - tile // 2
            if y0 < 0 or x0 < 0 or y0 + tile > H0 or x0 + tile > W0:
                continue
            rows = np.minimum((np.arange(y0, y0 + tile) / dy).astype(np.int64) - rs.start, occ.shape[0] - 1)
            cols = np.minimum((np.arange(x0, x0 + tile) / dx).astype(np.int64) - cs.start, occ.shape[1] - 1)
            if rows.min() < 0 or cols.min() < 0:
                continue
            if occ[np.ix_(rows[::8], cols[::8])].mean() < L0_TILE_MIN_MASK:
                continue
            if any(abs(y0 - oy) < tile and abs(x0 - ox) < tile for oy, ox in origins):
                continue
            vals = {ch: np.asarray(z[ci, y0:y0 + tile, x0:x0 + tile]) for ch, ci in channel_idx.items()}
            tiles.append(L0Tile(values=vals, rows=rows, cols=cols))
            origins.append((y0, x0))
        return tiles, np.asarray(origins)
    finally:
        try:
            store.close()
        except Exception:
            pass


def anchors_for_image(im: ProjectImage, channels=CHANNELS, manual_dir: Path = DEFAULT_MANUAL_DIR,
                      resolution_um: float = RESOLUTION_UM, sigma_um: float = SIGMA_UM,
                      tissue_geojson: Path | None = None, exclude_classes=EXCLUDE_CLASSES,
                      two_pass: bool = True, excl_tau_cd99: float = EXCL_TAU_CD99,
                      excl_tau_h: float = EXCL_TAU_H, excl_floor_k: float = EXCL_FLOOR_K,
                      excl_dilate_um: float = EXCL_DILATE_UM, with_smo: bool = True,
                      with_l0: bool = True, sample_max: int = SAMPLE_MAX,
                      n_l0_tiles: int = L0_TILES, l0_tile: int = L0_TILE) -> list[dict]:
    """One row per channel for one project image (see the module docstring for the columns)."""
    import tifffile
    t0 = time.time()
    gj = tissue_geojson or (manual_dir / f"{im.tag}.geojson")
    if not gj.exists():
        raise FileNotFoundError(f"no Tissue annotation export for {im.image}: {gj}")
    tissue_polys = read_annotation_polygons(gj, TISSUE_CLASSES)
    if not tissue_polys:
        raise ValueError(f"{gj} has no Tissue polygon")
    excl_polys = read_annotation_polygons(gj, exclude_classes) if exclude_classes else []
    hash_classes = tuple(TISSUE_CLASSES) + tuple(exclude_classes)
    ahash = annotation_hash(tissue_polys + excl_polys, hash_classes)
    phash = params_hash(estimation_params(resolution_um, sigma_um, exclude_classes, two_pass, excl_tau_cd99,
                                          excl_tau_h, excl_floor_k, excl_dilate_um, sample_max, n_l0_tiles,
                                          l0_tile))
    rows = []
    with tifffile.TiffFile(im.path) as tf:
        s = tf.series[0]
        full = s.levels[0].shape[-2:]
        if (full[0], full[1]) != (im.height, im.width):
            raise ValueError(f"{im.image}: TIFF {full} != project {im.height}x{im.width}")
        shapes = [lv.shape[-2:] for lv in s.levels]
        li, level_um, dx, dy = choose_level(shapes, full, im.pixel_um, resolution_um)
        lv = s.levels[li]
        mask_full = analysis_mask(tissue_polys, excl_polys, lv.shape[-2:], dx, dy, full,
                                  EDGE_MARGIN_UM / im.pixel_um)
        sigma_px = sigma_um / level_um
        margin = int(math.ceil(GAUSS_TRUNCATE * sigma_px)) + 2
        rs, cs = _crop_to_mask(mask_full, margin)
        mask = np.ascontiguousarray(mask_full[rs, cs])
        del mask_full
        n_mask = int(mask.sum())
        if n_mask == 0:
            raise ValueError(f"{im.image}: empty analysis mask")
        excluded_area = polygons_area_um2(excl_polys, im.pixel_um) / 1e6
        log.info("%s: level %d (%.4f um/px, %s), analysis area %.1f mm2 (excluded %.2f mm2 of %s), hash %s",
                 im.image, li, level_um, "x".join(map(str, lv.shape[-2:])), n_mask * level_um ** 2 / 1e6,
                 excluded_area, ",".join(exclude_classes) or "-", ahash)
        n_pages = len(lv.pages)
        idx = {ch: resolve_channel(im.channels, ch) for ch in channels}
        for ch, ci in idx.items():
            if ci >= n_pages:
                raise ValueError(f"{im.image}: channel {ch} index {ci} >= {n_pages} pages")
        blurred, zero_all, sat = {}, None, {}
        for ch, ci in idx.items():
            a = lv.pages[ci].asarray()[rs, cs]
            z = a == 0
            zero_all = z if zero_all is None else (zero_all & z)
            sat[ch] = float((a[mask] >= SATURATION).mean())
            blurred[ch] = gaussian_blur(a, sigma_px)
            del a
        # no DAPI mask: every analysis pixel that carries data is background-eligible
        nodata = zero_all
        occupied = ~nodata
        # level-0 tiles for the clamp model (data-carrying analysis pixels, crop coordinates)
        l0_tiles, zero_l0 = [], {}
        if with_l0:
            l0_tiles, _ = _sample_l0_tiles(s, idx, mask & occupied, rs, cs, dx, dy, n_tiles=n_l0_tiles, tile=l0_tile)
            for ch in idx:
                v0 = _l0_values(l0_tiles, ch, mask)
                zero_l0[ch] = float(np.mean(v0 == 0)) if v0.size else np.nan
            if len(l0_tiles) < n_l0_tiles:
                log.warning("%s: only %d of %d level-0 tiles fit inside the analysis mask", im.image,
                            len(l0_tiles), n_l0_tiles)
        res = compute_anchors(blurred, mask, occupied, level_um, l0_tiles, two_pass=two_pass,
                              excl_tau_cd99=excl_tau_cd99, excl_tau_h=excl_tau_h, excl_floor_k=excl_floor_k,
                              excl_dilate_um=excl_dilate_um, sample_max=sample_max)
        meta = res.pop("_meta")
        del blurred
        # v1.0 estimator (half-max at ~2 um, no blur) and SMO on the same mask, for comparison
        l2 = {}
        li2, um2, dx2, dy2 = choose_level(shapes, full, im.pixel_um, SMO_TARGET_UM)
        lv2 = s.levels[li2]
        m2_full = analysis_mask(tissue_polys, excl_polys, lv2.shape[-2:], dx2, dy2, full,
                                EDGE_MARGIN_UM / im.pixel_um)
        rs2, cs2 = _crop_to_mask(m2_full, 4 * SMO_SIZE)
        m2 = np.ascontiguousarray(m2_full[rs2, cs2])
        del m2_full
        for ch, ci in idx.items():
            a2 = lv2.pages[ci].asarray()[rs2, cs2]
            v2 = a2[m2]
            d = dict(anchor_l2_unblurred=log_background(v2[v2 > 0])["bg_halfmax"])
            if with_smo:
                d.update({k: v_ for k, v_ in smo_background(a2, ~m2 | (a2 >= SATURATION), m2 & (a2 > 0)).items()
                          if k != "smo_n"})
            else:
                d.update(bg_smo_mean=np.nan, bg_smo_std=np.nan, bg_smo_median=np.nan, smo_bg_fraction=np.nan)
            l2[ch] = d
        now = _dt.datetime.now().astimezone().isoformat(timespec="seconds")
        for ch, ci in idx.items():
            pg = s.levels[0].pages[ci]
            pg = pg.aspage() if hasattr(pg, "aspage") else pg      # OME-TIFF pages may be lazy TiffFrames
            meta_xml = parse_qpi_description(getattr(pg, "description", None))
            st = res[ch]
            row = dict(image=im.image, tag=im.tag, channel=ch, method_version=METHOD_VERSION,
                       annotation_hash=ahash, annotation_classes=",".join(hash_classes), params_hash=phash,
                       level_um=round(level_um, 4), level_index=li, sigma_um=sigma_um,
                       n_pixels=meta["n_pixels"], n_occupied=meta["n_occupied"], n_sampled=meta["n_sampled"],
                       n_sampled_pass2=meta["n_sampled_pass2"],
                       occupied_fraction=meta["n_occupied"] / max(meta["n_pixels"], 1),
                       nodata_fraction=float((nodata & mask).sum() / max(n_mask, 1)),
                       excluded_fraction=meta["excluded_fraction"], anchor_source=ANCHOR_SOURCE,
                       zero_fraction_l0=zero_l0.get(ch, np.nan), saturated_fraction=sat[ch],
                       excl_tau_cd99=excl_tau_cd99 if two_pass else np.nan,
                       excl_tau_h=excl_tau_h if two_pass else np.nan,
                       excl_floor_k=excl_floor_k if two_pass else np.nan,
                       excl_dilate_um=excl_dilate_um if two_pass else np.nan,
                       channel_index=ci, entry_id=im.entry_id,
                       image_path=str(im.path), annotation_source=str(gj),
                       n_tissue_polygons=len(tissue_polys), n_excluded_polygons=len(excl_polys),
                       excluded_area_mm2=excluded_area, analysis_area_mm2=n_mask * level_um ** 2 / 1e6,
                       computed_at=now)
            row.update({k: v_ for k, v_ in st.items() if k in COLUMNS})
            row.update(l2[ch])
            row.update(exposure_ms=meta_xml.get("exposure_ms"), xml_biomarker=meta_xml.get("xml_biomarker"),
                       responsivity=meta_xml.get("responsivity"), filter_date=meta_xml.get("filter_date"),
                       signal_units=meta_xml.get("signal_units"))
            xb = normalize_channel_name(meta_xml.get("xml_biomarker"))
            row["biomarker_match"] = (None if xb is None else
                                      bool(xb == ch or xb.lower() in [x.lower() for x in CHANNEL_ALIASES.get(ch, ())]))
            if row["biomarker_match"] is False:
                log.warning("%s: project channel %s (index %d) but QPTIFF Biomarker=%r",
                            im.image, ch, ci, meta_xml.get("xml_biomarker"))
            if not np.isfinite(row["anchor"]) or row["anchor"] <= 0:
                log.warning("%s %s: no valid anchor", im.image, ch)
            rows.append(row)
            log.info("  %-4s idx=%-2d anchor=%.2f (pass1 %.2f, v1-style %.2f) clamp: zero=%.3f mu=%.1f sigma=%.1f "
                     "ks=%.3f floor(k=%g)=%.1f",
                     ch, ci, row["anchor"], row["anchor_pass1"], row["anchor_l2_unblurred"],
                     row["clamp_fraction"], row["clamp_mu"], row["clamp_sigma"], row["clamp_fit_ks"],
                     excl_floor_k, row["clamp_mu"] + excl_floor_k * row["clamp_sigma"])
    log.info("%s done in %.1f s (with data %.2f, excluded %.3f of the mask; %d level-0 tiles)", im.image,
             time.time() - t0, meta["n_occupied"] / max(n_mask, 1), meta["excluded_fraction"], len(l0_tiles))
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
        old = pd.read_csv(out, dtype=HASH_DTYPES)   # hashes are text: "0123..." or "12e45..." must not become numbers
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
        os.chmod(tmp, 0o664)
        os.replace(tmp, out)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)
    return df


def existing_images(out: Path, method_version: str = METHOD_VERSION, channels=CHANNELS,
                    current_hash: dict | None = None, current_params: str | None = None) -> set[str]:
    """Images already in `out` with every channel at this method_version and, when given, the same
    annotation hash (current_hash: image -> hash) and estimation-parameter hash."""
    import pandas as pd
    if not out.exists() or out.stat().st_size == 0:
        return set()
    df = pd.read_csv(out, dtype=HASH_DTYPES)
    ok = df[df["method_version"].astype(str) == str(method_version)]
    if current_params is not None and "params_hash" in ok:
        ok = ok[ok["params_hash"].astype(str) == current_params]
    if current_hash is not None and "annotation_hash" in ok:
        ok = ok[[current_hash.get(i) == str(h) for i, h in zip(ok["image"], ok["annotation_hash"])]]
    counts = ok.groupby("image")["channel"].nunique()
    return set(counts[counts >= len(channels)].index)


def current_annotation_hash(im: ProjectImage, manual_dir: Path, tissue_geojson: Path | None,
                            exclude_classes=EXCLUDE_CLASSES) -> str | None:
    gj = tissue_geojson or (manual_dir / f"{im.tag}.geojson")
    if not gj.exists():
        return None
    polys = read_annotation_polygons(gj, tuple(TISSUE_CLASSES) + tuple(exclude_classes))
    return annotation_hash(polys, tuple(TISSUE_CLASSES) + tuple(exclude_classes))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--images", nargs="+", help="substrings of project image names")
    g.add_argument("--all", action="store_true", help="every project image")
    ap.add_argument("--project", type=Path, default=DEFAULT_PROJECT)
    ap.add_argument("--manual-dir", type=Path, default=DEFAULT_MANUAL_DIR,
                    help="folder of <tag>.geojson annotation exports (Tissue, pLN, Artifact)")
    ap.add_argument("--tissue-geojson", type=Path,
                    help="explicit annotation GeoJSON for a single image (e.g. 6457, whose .qpdata is "
                         "empty: results/Geojson/Tissue/Tissue__6457_...geojson)")
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--channels", default=",".join(CHANNELS))
    ap.add_argument("--resolution-um", type=float, default=RESOLUTION_UM,
                    help="the detector's RESOLUTION_UM (nearest pyramid level is used)")
    ap.add_argument("--sigma-um", type=float, default=SIGMA_UM, help="the detector's SIGMA_UM")
    ap.add_argument("--exclude-classes", default=",".join(EXCLUDE_CLASSES),
                    help="annotation classes cut from the analysis area (the detector's EXCLUDE_CLASSES)")
    ap.add_argument("--single-pass", action="store_true", help="skip the endocrine exclusion (pass 2)")
    ap.add_argument("--excl-tau-cd99", type=float, default=EXCL_TAU_CD99)
    ap.add_argument("--excl-tau-h", type=float, default=EXCL_TAU_H)
    ap.add_argument("--excl-floor-k", type=float, default=EXCL_FLOOR_K)
    ap.add_argument("--excl-dilate-um", type=float, default=EXCL_DILATE_UM)
    ap.add_argument("--sample-max", type=int, default=SAMPLE_MAX)
    ap.add_argument("--l0-tiles", type=int, default=L0_TILES)
    ap.add_argument("--skip-existing", action="store_true",
                    help="skip images whose rows are at this method_version with the same annotation "
                         "hash and estimation parameters (a changed annotation export is recomputed)")
    ap.add_argument("--no-smo", action="store_true", help="skip the SMO diagnostic (faster)")
    ap.add_argument("--no-l0", action="store_true", help="skip level-0 sampling (no clamp model)")
    ap.add_argument("-v", "--verbose", action="store_true")
    a = ap.parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if a.verbose else logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s", datefmt="%H:%M:%S")
    channels = tuple(c.strip() for c in a.channels.split(",") if c.strip())
    bad = [c for c in channels if c in FORBIDDEN_CHANNELS]
    if bad:
        raise SystemExit(f"refusing forbidden channels: {bad}")
    exclude = tuple(c.strip() for c in a.exclude_classes.split(",") if c.strip())
    images = select_images(read_project(a.project), None if a.all else a.images)
    if a.tissue_geojson and len(images) != 1:
        raise SystemExit("--tissue-geojson needs exactly one selected image")
    kw = dict(resolution_um=a.resolution_um, sigma_um=a.sigma_um, exclude_classes=exclude,
              two_pass=not a.single_pass, excl_tau_cd99=a.excl_tau_cd99, excl_tau_h=a.excl_tau_h,
              excl_floor_k=a.excl_floor_k, excl_dilate_um=a.excl_dilate_um, sample_max=a.sample_max,
              n_l0_tiles=a.l0_tiles)
    done = set()
    if a.skip_existing:
        cur = {im.image: current_annotation_hash(im, a.manual_dir, a.tissue_geojson, exclude) for im in images}
        phash = params_hash(estimation_params(**kw))
        done = existing_images(a.out, METHOD_VERSION, channels, cur, phash)
    n_ok = n_fail = 0
    for im in images:
        if im.image in done:
            log.info("skip %s (method_version %s, same annotation hash and parameters)", im.image, METHOD_VERSION)
            continue
        try:
            rows = anchors_for_image(im, channels, a.manual_dir, tissue_geojson=a.tissue_geojson,
                                     with_smo=not a.no_smo, with_l0=not a.no_l0, **kw)
        except Exception as exc:  # keep going; report at the end
            log.exception("FAILED %s: %s", im.image, exc)
            n_fail += 1
            continue
        upsert_csv(rows, a.out)   # after every image, so an interrupted --all keeps its work
        n_ok += 1
    log.info("anchors: %d image(s) written to %s, %d failed", n_ok, a.out, n_fail)
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
