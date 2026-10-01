#!/usr/bin/env python3
"""islet_classifier_fit.py -- one islet pixel classifier per image, set from that image's own cells.

WHY (maintainer, 2026-10-01): "All I wanted was the pixel classifier to be more accurate and not filter out so
many cells while adapting to batch effects that changed the classifier for each image. I could go through each
image and do that manually, but I asked you to do it. Then the rest of the workflow works without major changes."

So this script does, per image, what would otherwise be done by hand per image: it sets the islet pixel
classifier. Nothing else in the workflow changes. qupath/1_panc_objects_perimage.groovy is the maintainer's
1_panc_objects.groovy loading "Islets_<tag>" where it loaded "Islets2"; its size-adaptive erode / dilate, its
50 um^2 minimum and every later script are the maintainer's.

WHAT IS WRITTEN. <project>/classifiers/pixel_classifiers/Islets_<tag>.json -- an ordinary QuPath pixel classifier,
at Islets2's resolution (1.0156 um) and on Islets2's four channels, loadable in the GUI (Classify > Pixel
classification > Load pixel classifier) and by loadPixelClassifier():

    CD99, INS, GCG, SST  ->  each divided by ITS OWN threshold for this image
                         ->  Gaussian, sigma 1 um
                         ->  HORMONE pixels: INS, GCG or SST above its threshold
                             ANY pixels:     CD99 or a hormone above its threshold
                         ->  islet = an ANY pixel within SEED_PX of a HORMONE pixel
                                                               (CD99 widens an islet a hormone has found; CD99 alone
                                                                never makes one)
                         ->  closing, radius 3 px              (covers the nucleus of a hormone-positive cell: a cell
                                                                belongs to the object its CENTROID is in)
                         ->  opening, radius 1 px              (a speck or a thread one pixel wide is not an islet)
                         ->  margin, radius 1 px               (a fixed 1 um, the same for a dim and a bright cell, so
                                                                that a single small cell survives script 1's 50 um^2
                                                                minimum and its 1 um erosion)

Only the four thresholds differ between images. Islets2 is one set of raw-count weights for every image, so it
moves with the acquisition batch; and its sigma-4 features put its outline 3-5 um outside the cells, which is the
ring of neighbours the maintainer rejected ("Too many false positives handling segmentations around
endocrine-positive cells").

WHY CD99 NEEDS A HORMONE NEARBY. A first version let any of the four channels make an islet. On the control donors
(6482, HDL112) it outlined CD99-bright membranes in exocrine tissue with no hormone anywhere near -- thousands of
objects that are not islets. Islets2 has the same property in its own way: CD99 alone must be very bright, CD99
beside a hormone counts.

WHY NOT QUPATH'S OWN closing / minimum / maximum FILTERS. They share one kernel between threads and it is freed
when the tile that created it finishes (ImageOps.Filters.MorphOp, QuPath 0.7.0): on 2026-10-01 a classifier
written with them stopped script 1 on 6533 with "Pointer address of argument 3 is NULL" inside morphologyEx. Every
shaping step here is QuPath's `sum` filter (a disk from a synchronized, retained cache) followed by a threshold: a
dilation is "any pixel of the disk is on", an erosion "every pixel of the disk is on".

HOW THE FOUR THRESHOLDS ARE SET. From the image's own cells, used ONLY here (no cell test is a workflow step):

  endocrine  cells this image itself confirms: a hormone (INS / GCG / SST) is the cell's highest identity marker
             AND its CD99 is at or above the image's CD99 reference level ("confirmed"); or, within TERRITORY_UM
             of a confirmed cell, a cell whose CD99 reaches the 10th percentile of the confirmed cells' CD99 and
             whose highest marker is not another lineage's.
  other      no hormone highest, and either CD99 below the reference level or the cell is ordinary tissue (farther
             than TERRITORY_UM from every hormone-highest cell -- the population the reference level is drawn
             from, whatever its CD99). Those enclosed by endocrine cells (an islet's own capillaries) are not
             counted against the classifier -- they belong inside.
  left out   everything in between (a hormone highest without CD99; CD99 between the two levels, or on a vessel or
             immune cell, beside an islet): they neither pull the thresholds down nor push them up.

  The thresholds are the four numbers that minimise
        (other cells whose centroid the classifier covers) + MISS_WEIGHT x (endocrine cells it does not cover).

  Values are the source workflow's REDSEA-corrected cell means, as the marker registry declares for these four.
  The CD99 reference level is the typing lane's: the 99th percentile over ordinary tissue -- cells with a nucleus,
  in Tissue, in no Nerve / Blood vessel / Lymphatic object, and farther than TERRITORY_UM from every cell whose
  highest marker is a hormone. Fitted on this image alone, never pooled.

WHAT IS SHARED BY EVERY IMAGE. The shape of the chain (sigma, the seed distance, closing, opening, margin) and
MISS_WEIGHT. They were chosen on the non-diabetic / control donors only (maintainer, 2026-10-01: "use all
non-diabetic/controls for testing classifier since the endocrine differences with disease interfere"), scored after
script 1's own steps. The chain was applied by QuPath 0.7 itself to a real crop and compared with this file's own
arithmetic: 0 pixels differ.

REFUSALS. A channel that is not in the image is unavailable, not zero. A level or a threshold at the floor is
refused (it would admit every cell while looking fitted). An image with fewer than MIN_CONFIRMED confirmed
endocrine cells, or one where no hormone ends up with a threshold, gets no classifier: that is reported, not
papered over.

WHAT THE REPORT SAYS. The counts above, the thresholds, and what 1_panc_objects.groovy will make of the classifier
(its 50 um^2 minimum and size-adaptive erode / dilate are emulated on the classifier's own grid): cells inside,
how many of them are `other` (at the edge / enclosed), how many endocrine cells are left outside. These are
in-sample numbers about the fit. The check that counts is looking at the image.

    islet_classifier_fit.py --tag TAG --image IMG.qptiff --cells cells/TAG.parquet --redsea redsea/donor_id=D/data_0.parquet \\
        --structures <export>/structures --classifier <project>/classifiers/pixel_classifiers/Islets_TAG.json \\
        --report out/TAG.json [--mask out/TAG.mask.npy] [--objects out/TAG.objects.npy]
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import time
from pathlib import Path

import cv2
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy import ndimage as ndi
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parent))
import islet_cell_test as ict  # noqa: E402  (one implementation of the level, the highest marker, the object readers)

CHANNELS = ("CD99", "INS", "GCG", "SST")
HORMONES = ict.HORMONES
CLASSIFIER_UM = 1.0156          # Islets2's resolution: level 1 of the pyramid (2 x 0.5078)
SIGMA_UM = 1.0
SEED_PX = 6                     # CD99 counts within this distance of a hormone pixel
CLOSE_PX = 3
OPEN_PX = 1
MARGIN_PX = 1
MISS_WEIGHT = 8.0               # an endocrine cell left outside counts as eight neighbours let in
TERRITORY_UM = 30.0             # cells this close to a hormone-highest cell are not "ordinary tissue"
INTERIOR_UM = 20.0              # the neighbourhood that says whether an `other` cell is enclosed by endocrine cells
CONFIRMED_QUANTILE = 10.0       # of the confirmed cells' CD99: the level a hormone-less cell must reach
MIN_CONFIRMED = 50
OTHER_LINEAGES = ("Vascular", "T_NK", "B_Plasma", "Myeloid", "Granulocyte", "Neural")
MIN_AREA_UM2 = 50.0             # 1_panc_objects.groovy: createAnnotationsFromPixelClassifier(..., 50.0, -1.0, SPLIT)
EROSION_UM = ((2000.0, 1.0), (10000.0, 2.0), (np.inf, 3.0))   # 1_panc_objects.groovy: area below -> erosion
DILATION_FACTOR = 0.85
GRID = 240


# ----------------------------------------------------------------------------------------------- the image

def read_planes(image: Path, channels=CHANNELS) -> dict:
    """The classifier's input: level 1 of the pyramid (what QuPath reads at 1.0156 um), one float32 plane per channel."""
    import tifffile  # noqa: PLC0415
    with tifffile.TiffFile(image) as tf:
        series = tf.series[0]
        names = []
        for p in series.levels[0].pages:
            p = p.aspage() if hasattr(p, "aspage") else p
            m = re.search(r"<Biomarker>(.*?)</Biomarker>", p.description or "")
            names.append(m.group(1) if m else "?")
        missing = [c for c in channels if c not in names]
        if missing:
            raise SystemExit(f"{image.name}: channel(s) {missing} are not in this image ({len(names)} channels) -- unavailable, not zero")
        level = series.levels[1]
        scale = series.levels[0].shape[-1] / level.shape[-1]
        if abs(scale - 2.0) > 0.01:
            raise SystemExit(f"{image.name}: pyramid level 1 is {scale:.3f}x, not 2x -- the classifier's resolution would not be a stored level")
        return {c: level.pages[names.index(c)].asarray().astype(np.float32) for c in channels}


def smooth(plane: np.ndarray, sigma_px: float) -> np.ndarray:
    return cv2.GaussianBlur(plane, (0, 0), sigmaX=sigma_px, sigmaY=sigma_px, borderType=cv2.BORDER_REFLECT)


def disk(radius: int) -> np.ndarray:
    """The kernel of QuPath's `sum` filter (OpenCVTools.createDisk; printed by QuPath 0.7 itself for radius 1, 2, 3):
    every pixel within `radius` of the centre -- 5, 13, 29, 49 pixels for radius 1, 2, 3, 4."""
    yy, xx = np.ogrid[-radius:radius + 1, -radius:radius + 1]
    return (yy * yy + xx * xx <= radius * radius).astype(np.uint8)


def dilate(img: np.ndarray, radius: int) -> np.ndarray:
    return cv2.dilate(img, disk(radius), borderType=cv2.BORDER_REPLICATE) if radius > 0 else img


def erode(img: np.ndarray, radius: int) -> np.ndarray:
    return cv2.erode(img, disk(radius), borderType=cv2.BORDER_REPLICATE) if radius > 0 else img


def shaped(img: np.ndarray, close_px: int, open_px: int, margin_px: int) -> np.ndarray:
    """closing -> opening -> margin with disks, as the classifier file does them.

    Works on a 0/1 mask and on a grey image alike: a flat dilation / erosion and a threshold commute, so thresholding
    the shaped grey image at T is "threshold at T, then shape" -- which is what lets one pass serve every T in the fit.
    """
    out = erode(dilate(img, close_px), close_px)
    out = dilate(erode(out, open_px), open_px)
    return dilate(out, margin_px)


def classifier_mask(smoothed: dict, thresholds: dict, *, seed_px: int = SEED_PX, close_px: int = CLOSE_PX,
                    open_px: int = OPEN_PX, margin_px: int = MARGIN_PX) -> np.ndarray:
    """What the classifier file computes: True where the pixel is Islet. `smoothed` = Gaussian of each channel."""
    def above(channels):
        comb = None
        for c in channels:
            scaled = smoothed[c] * np.float32(1.0 / thresholds[c])
            comb = scaled if comb is None else np.maximum(comb, scaled, out=comb)
        return (comb > 1.0).astype(np.uint8)

    hormones = [c for c in smoothed if c in HORMONES and thresholds.get(c) is not None]
    if not hormones:
        raise SystemExit("no hormone carries a threshold: there is no classifier to apply")
    mask = above(hormones)
    if "CD99" in smoothed and thresholds.get("CD99") is not None:
        mask = above(["CD99", *hormones]) & dilate(mask, seed_px)
    return shaped(mask, close_px, open_px, margin_px).astype(bool)


# ----------------------------------------------------------------------------------------------- the cells

def enclosed_by(xy: np.ndarray, candidates: np.ndarray, endocrine: np.ndarray, radius_um: float = INTERIOR_UM) -> np.ndarray:
    """True for a candidate cell that endocrine cells surround: it lies inside the convex hull of the endocrine cells
    within `radius_um` (they are not all on one side of it). An islet's own capillary is; the neighbour outside is not."""
    out = np.zeros(len(xy), dtype=bool)
    if not candidates.any() or int(endocrine.sum()) < 3:
        return out
    end_xy = xy[endocrine]
    idx = np.flatnonzero(candidates)
    near = cKDTree(end_xy).query_ball_point(xy[idx], radius_um)
    for i, hits in zip(idx, near):
        if len(hits) < 3:
            continue
        d = end_xy[hits] - xy[i]
        angles = np.sort(np.arctan2(d[:, 1], d[:, 0]))
        gaps = np.diff(np.append(angles, angles[0] + 2 * np.pi))
        out[i] = gaps.max() < np.pi
    return out


def within(xy: np.ndarray, of: np.ndarray, radius_um: float) -> np.ndarray:
    """True for every row of `xy` within `radius_um` of a row flagged in `of`."""
    if not of.any():
        return np.zeros(len(xy), dtype=bool)
    dist, _ = cKDTree(xy[of]).query(xy, k=1, distance_upper_bound=radius_um)
    return np.isfinite(dist)


def label_cells(values: pd.DataFrame, x_um: np.ndarray, y_um: np.ndarray, eligible: np.ndarray, in_structure: np.ndarray,
                other_lineage_markers: set, quantile: float) -> dict:
    """endocrine / other / left out, and the levels they rest on. `values`: one row per cell, identity markers + CD99."""
    top = ict.top_marker(values)
    cd99 = values["CD99"].to_numpy(dtype=float)
    measured = np.isfinite(cd99)
    has_top = top != ""
    hormone_top = eligible & np.isin(top, HORMONES)
    xy = np.column_stack([x_um, y_um])
    territory = within(xy, hormone_top, TERRITORY_UM)
    reference = eligible & ~in_structure & ~territory
    level = ict.reference_level(cd99[reference], quantile)
    confirmed = hormone_top & measured & (cd99 >= level)
    if int(confirmed.sum()) < MIN_CONFIRMED:
        raise SystemExit(f"only {int(confirmed.sum())} cells have a hormone highest AND CD99 at or above the reference level ({level:.1f}); "
                         f"fewer than {MIN_CONFIRMED} cannot set a classifier -- this image needs to be looked at")
    by_cd99_level = max(level, float(np.percentile(cd99[confirmed], CONFIRMED_QUANTILE)))
    by_cd99 = (eligible & ~hormone_top & has_top & measured & (cd99 >= by_cd99_level) & ~np.isin(top, list(other_lineage_markers))
               & within(xy, confirmed, TERRITORY_UM))
    endocrine = confirmed | by_cd99
    # ordinary tissue is `other` whatever its CD99: it is the population the level was drawn from
    other = eligible & ~hormone_top & has_top & ~by_cd99 & ((measured & (cd99 < level)) | ~territory)
    # an `other` cell enclosed by endocrine cells is the islet's own (capillary, pericyte, immune cell): not counted against
    enclosed = enclosed_by(xy, other, endocrine)
    return {"top": top, "hormone_top": hormone_top, "reference": reference, "confirmed": confirmed, "by_cd99": by_cd99,
            "endocrine": endocrine, "other": other & ~enclosed, "enclosed": enclosed,
            "left_out": eligible & ~endocrine & ~other, "cd99_reference_level": level, "cd99_confirmed_level": by_cd99_level}


def centroid_values(smoothed: dict, cx: np.ndarray, cy: np.ndarray, *, seed_px: int = SEED_PX, close_px: int = CLOSE_PX,
                    open_px: int = OPEN_PX, margin_px: int = MARGIN_PX) -> tuple[np.ndarray, np.ndarray]:
    """Per cell and channel, two numbers a threshold is compared with:

    direct  the highest threshold at which that channel ALONE would still cover the cell's centroid
    seed    (hormones) the highest threshold at which the channel still has a pixel within `seed_px` of the centroid --
            what lets CD99 count there. NaN in the CD99 column: CD99 seeds nothing.
    """
    direct = np.column_stack([shaped(smoothed[c], close_px, open_px, margin_px)[cy, cx] for c in CHANNELS])
    seed = np.column_stack([np.full(len(cx), np.nan) if c == "CD99" else dilate(smoothed[c], seed_px)[cy, cx] for c in CHANNELS])
    return direct, seed


def covered(direct: np.ndarray, seed: np.ndarray, t: np.ndarray) -> np.ndarray:
    """A centroid is covered when a hormone covers it, or CD99 does within reach of a hormone pixel."""
    hormone = (direct[:, 1:] > t[None, 1:]).any(axis=1)
    seeded = (seed[:, 1:] > t[None, 1:]).any(axis=1)
    return hormone | ((direct[:, 0] > t[0]) & seeded)


def fit_thresholds(direct: np.ndarray, seed: np.ndarray, endocrine: np.ndarray, other: np.ndarray,
                   miss_weight: float = MISS_WEIGHT, grid: int = GRID) -> tuple[np.ndarray, float]:
    """Thresholds for (CD99, hormones...) -- np.inf = channel not used -- minimising
    #other covered + miss_weight * #endocrine not covered, with `covered` as above.

    Coordinate descent from "no channel used", over each channel's own grid; the best of several channel orders is kept
    (deterministic). Across a run of equal costs the middle threshold is taken, not its edge.
    """
    n_ch = direct.shape[1]
    keep = endocrine | other
    direct, seed, endocrine, other = direct[keep], seed[keep], endocrine[keep], other[keep]
    grids = []
    for j in range(n_ch):
        col = direct[:, j]
        lo = float(np.percentile(col[other], 50)) if other.any() else float(np.nanmin(col))
        hi = float(np.percentile(col[endocrine], 99.5)) if endocrine.any() else float(np.nanmax(col))
        lo = max(lo, 1e-6)
        if not np.isfinite(hi) or hi <= lo:
            grids.append(np.array([np.inf]))
            continue
        grids.append(np.append(np.exp(np.linspace(np.log(lo), np.log(hi), grid)), np.inf))

    def loss(t: np.ndarray) -> float:
        c = covered(direct, seed, t)
        return float((c & other).sum() + miss_weight * (~c & endocrine).sum())

    hormones = list(range(1, n_ch))
    best_t, best_loss = None, np.inf
    for order in (hormones + [0], [0] + hormones, hormones[::-1] + [0]):
        t = np.full(n_ch, np.inf)
        for _ in range(12):
            changed = False
            for j in order:
                if j == 0:      # CD99: counts only where a hormone seeds it; cells a hormone covers are settled
                    settled = (direct[:, 1:] > t[None, 1:]).any(axis=1)
                    seeded = (seed[:, 1:] > t[None, 1:]).any(axis=1)
                    critical = np.where(seeded, direct[:, 0], -np.inf)
                else:           # a hormone: covers directly, and seeds the cells CD99 already reaches
                    others = [k for k in hormones if k != j]
                    cd99_on = direct[:, 0] > t[0]
                    settled = (direct[:, others] > t[None, others]).any(axis=1) | (cd99_on & (seed[:, others] > t[None, others]).any(axis=1))
                    critical = np.where(cd99_on, np.maximum(direct[:, j], seed[:, j]), direct[:, j])
                free = ~settled
                vo = np.sort(critical[free & other]); ve = np.sort(critical[free & endocrine])
                # for threshold g: other covered = #(vo > g); endocrine missed = #(ve <= g)
                cost = (len(vo) - np.searchsorted(vo, grids[j], side="right")) + miss_weight * np.searchsorted(ve, grids[j], side="right")
                first = int(np.argmin(cost))
                if cost[-1] <= cost[first]:                   # using the channel gains nothing: it is not used
                    g = np.inf
                else:
                    last = first
                    while last + 1 < len(cost) and cost[last + 1] == cost[first] and np.isfinite(grids[j][last + 1]):
                        last += 1
                    g = grids[j][(first + last) // 2]
                if g != t[j]:
                    t[j] = g; changed = True
            if not changed:
                break
        value = loss(t)
        if value < best_loss:
            best_t, best_loss = t.copy(), value
    for j in range(n_ch):
        if np.isfinite(best_t[j]) and best_t[j] <= grids[j][0] * (1 + 1e-9):
            raise SystemExit(f"channel {j}: the fitted threshold sits at the floor of its range ({best_t[j]:.3g}) -- refusing to write a classifier")
    if not np.isfinite(best_t[1:]).any():
        raise SystemExit("no hormone ends up with a threshold: the classifier would have nothing to start an islet from -- refusing to write it")
    return best_t, best_loss


# ----------------------------------------------------------------------------------------------- script 1, emulated

def script1_objects(mask: np.ndarray, px_um: float = CLASSIFIER_UM) -> np.ndarray:
    """What 1_panc_objects.groovy makes of a classifier's mask, on the classifier's own grid (a label image).

    50 um^2 minimum at creation; erosion by 1 / 2 / 3 um by the object's area, then dilation by 0.85 x that;
    fragments under 50 um^2 removed. QuPath does this on vector outlines; on a 1 um grid it is an approximation
    good to about half a pixel, which is what the counts in the report are worth.
    """
    lab, n = ndi.label(mask)
    out = np.zeros(mask.shape, dtype=np.int32)
    if n == 0:
        return out
    areas = np.bincount(lab.ravel(), minlength=n + 1) * px_um * px_um
    last = 0
    for i, sl in enumerate(ndi.find_objects(lab), start=1):
        if sl is None or areas[i] < MIN_AREA_UM2:
            continue
        erosion = next(e for limit, e in EROSION_UM if areas[i] < limit)
        pad = int(np.ceil(erosion / px_um)) + 2
        y0, y1 = max(sl[0].start - pad, 0), min(sl[0].stop + pad, mask.shape[0])
        x0, x1 = max(sl[1].start - pad, 0), min(sl[1].stop + pad, mask.shape[1])
        sub = lab[y0:y1, x0:x1] == i
        inside = (ndi.distance_transform_edt(sub) - 0.5) * px_um          # um from the outline, inside
        eroded = sub & (inside > erosion)
        if not eroded.any():
            continue
        # an object the erosion cuts in two stays two objects; each part is dilated on its own
        parts, k = ndi.label(eroded)
        dist, (iy, ix) = ndi.distance_transform_edt(~eroded, return_indices=True)
        grown = (dist - 0.5) * px_um <= DILATION_FACTOR * erosion         # um from the eroded outline, outside (0 inside)
        part_of = parts[iy, ix]
        sizes = np.bincount(part_of[grown], minlength=k + 1) * px_um * px_um
        keep = sizes >= MIN_AREA_UM2                                      # RefineAnnotations: fragments under 50 um^2 go
        keep[0] = False
        new_id = np.zeros(k + 1, dtype=np.int32)
        new_id[keep] = last + np.arange(1, int(keep.sum()) + 1)
        last += int(keep.sum())
        target = out[y0:y1, x0:x1]
        write = grown & keep[part_of]
        target[write] = new_id[part_of[write]]
    return out


# ----------------------------------------------------------------------------------------------- the classifier file

def classifier_json(thresholds: dict, *, sigma_um: float = SIGMA_UM, seed_px: int = SEED_PX, close_px: int = CLOSE_PX,
                    open_px: int = OPEN_PX, margin_px: int = MARGIN_PX, px_um: float = CLASSIFIER_UM) -> dict:
    """The QuPath pixel classifier (see the module docstring for why the shaping is `sum` + threshold)."""
    hormones = [c for c in HORMONES if thresholds.get(c) is not None]
    if not hormones:
        raise SystemExit("no hormone carries a threshold: refusing to write a classifier")
    cd99 = thresholds.get("CD99") is not None
    used = (["CD99"] if cd99 else []) + hormones
    sigma_px = sigma_um / px_um

    def dilation(radius: int) -> list:        # on if any pixel of the disk is on
        return [{"type": "op.filters.sum", "radius": int(radius)}, {"type": "op.threshold.constant", "thresholds": [0.5]}]

    def erosion(radius: int) -> list:         # on if every pixel of the disk is on
        return [{"type": "op.filters.sum", "radius": int(radius)},
                {"type": "op.threshold.constant", "thresholds": [float(disk(radius).sum()) - 0.5]}]

    def any_above(n_channels: int) -> list:   # 0 / 1: some channel is above its threshold
        return ([{"type": "op.channels.maximum"}] if n_channels > 1 else []) + [{"type": "op.threshold.constant", "thresholds": [1.0]}]

    ops = [{"type": "op.core.multiply", "values": [1.0 / float(thresholds[c]) for c in used]},
           {"type": "op.filters.gaussian", "sigmaX": sigma_px, "sigmaY": sigma_px}]
    if cd99:
        hormone_branch = [{"type": "op.channels.extract-channels", "channels": list(range(1, len(used)))}, *any_above(len(hormones)),
                          *dilation(seed_px)]
        ops += [{"type": "op.core.split-merge", "ops": [{"type": "op.core.sequential", "ops": hormone_branch},
                                                        {"type": "op.core.sequential", "ops": any_above(len(used))}]},
                {"type": "op.channels.sum"}, {"type": "op.threshold.constant", "thresholds": [1.5]}]      # both: ANY and near a HORMONE
    else:
        ops += any_above(len(used))
    if close_px > 0:
        ops += dilation(close_px) + erosion(close_px)
    if open_px > 0:
        ops += erosion(open_px) + dilation(open_px)
    if margin_px > 0:
        ops += dilation(margin_px)
    ops.append({"type": "op.core.convert", "pixelType": "UINT8"})
    size = {"value": px_um, "unit": "µm"}
    return {
        "pixel_classifier_type": "OpenCVPixelClassifier",
        "metadata": {
            "inputPadding": 0,
            "inputResolution": {"pixelWidth": size, "pixelHeight": size, "zSpacing": {"value": 1.0, "unit": "z-slice"},
                                "timeUnit": "SECONDS", "timepoints": []},
            "inputWidth": 512, "inputHeight": 512, "inputNumChannels": 3,
            "outputType": "CLASSIFICATION", "outputChannels": [],
            "classificationLabels": {"0": {"name": "Ignore*", "color": [180, 180, 180]}, "1": {"name": "Islet", "color": [255, 0, 0]}},
        },
        "op": {"type": "data.op.channels", "colorTransforms": [{"channelName": c} for c in used],
               "op": {"type": "op.core.sequential", "ops": ops}},
    }


# ----------------------------------------------------------------------------------------------- one image

def strata(counts: np.ndarray) -> dict:
    return {"1": int((counts == 1).sum()), "2-3": int(((counts >= 2) & (counts <= 3)).sum()),
            "4-10": int(((counts >= 4) & (counts <= 10)).sum()), ">10": int((counts > 10).sum())}


def single_cell_objects(which: np.ndarray, inside: np.ndarray, per: np.ndarray, lab: dict) -> dict:
    """What the one cell of a one-cell object is."""
    lone = inside & (np.append(0, per)[which] == 1)
    return {"endocrine": int((lone & lab["endocrine"]).sum()), "other": int((lone & (lab["other"] | lab["enclosed"])).sum()),
            "left_out": int((lone & lab["left_out"] & ~lab["enclosed"]).sum())}


def score(objects: np.ndarray, cx: np.ndarray, cy: np.ndarray, lab: dict, eligible: np.ndarray, px_um: float = CLASSIFIER_UM) -> dict:
    which = objects[cy, cx]
    inside = eligible & (which > 0)
    per = np.bincount(which[inside], minlength=int(objects.max()) + 1)[1:]
    pct = lambda a, b: round(100.0 * a / b, 2) if b else None  # noqa: E731
    n_in = int(inside.sum())
    return {
        "objects": int(objects.max()),
        "objects_with_cells": int((per > 0).sum()), "objects_by_cells": strata(per[per > 0]),
        "area_mm2": round(float((objects > 0).sum() * px_um * px_um / 1e6), 4),
        "cells_inside": n_in,
        "inside_endocrine": int((inside & lab["endocrine"]).sum()),
        "inside_other_at_the_edge": int((inside & lab["other"]).sum()),
        "inside_other_enclosed": int((inside & lab["enclosed"]).sum()),
        "inside_left_out": int((inside & lab["left_out"]).sum()),
        "inside_other_at_the_edge_percent": pct(int((inside & lab["other"]).sum()), n_in),
        "endocrine_outside": int((~inside & lab["endocrine"]).sum()),
        "endocrine_outside_percent": pct(int((~inside & lab["endocrine"]).sum()), int(lab["endocrine"].sum())),
        "confirmed_outside": int((~inside & lab["confirmed"]).sum()),
        "confirmed_outside_percent": pct(int((~inside & lab["confirmed"]).sum()), int(lab["confirmed"].sum())),
        "hormone_highest_without_cd99_inside": int((inside & lab["hormone_top"] & ~lab["confirmed"]).sum()),
        "hormone_highest_without_cd99_outside": int((eligible & ~inside & lab["hormone_top"] & ~lab["confirmed"]).sum()),
        "single_cell_objects": single_cell_objects(which, inside, per, lab),
    }


def load_cells(cells: Path, redsea: Path, structures: Path | None, tag: str | None, senior_typing: Path):
    bcf, identity, quantile = ict.senior(senior_typing)
    import lineages  # noqa: PLC0415  (on the path after ict.senior)
    other_markers = {m for name in OTHER_LINEAGES for m in lineages.LINEAGE_MARKERS[name]}
    base = pq.read_table(cells, columns=["Object ID", "Parent", "Centroid X µm", "Centroid Y µm", "Nucleus: Area µm^2"]).to_pandas()
    if base["Object ID"].duplicated().any():
        raise SystemExit(f"{cells}: duplicate Object ID")
    base = base.set_index("Object ID")
    red = pq.read_table(redsea).to_pandas().set_index("object_id")
    wanted = [*identity, "CD99"]
    have = [m for m in wanted if f"Cell__{m}" in red.columns and red[f"Cell__{m}"].notna().any()]
    if "CD99" not in have:
        raise SystemExit(f"{redsea}: CD99 carries no value -- unavailable, the classifier cannot be set")
    values = red[[f"Cell__{m}" for m in have]].rename(columns=lambda c: c.split("__", 1)[1]).reindex(base.index)
    x_um = base["Centroid X µm"].to_numpy(dtype=float); y_um = base["Centroid Y µm"].to_numpy(dtype=float)
    eligible = base["Nucleus: Area µm^2"].notna().to_numpy() & (base["Parent"] == ict.TISSUE_PARENT).to_numpy()
    in_structure = np.zeros(len(base), dtype=bool)
    if structures is not None:
        px0 = CLASSIFIER_UM / 2.0
        for name in ict.STRUCTURES:
            f = structures / name / f"{name}__{tag}.geojson"
            if not f.exists():
                raise SystemExit(f"missing structure objects: {f}")
            in_structure |= ict.centroid_in_any(x_um / px0, y_um / px0, ict.read_objects(f), bcf.subdivide)
    return base, values, x_um, y_um, eligible, in_structure, other_markers, quantile, [m for m in wanted if m not in have]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--image", type=Path, required=True)
    ap.add_argument("--cells", type=Path, required=True)
    ap.add_argument("--redsea", type=Path, required=True, help="data/redsea/donor_id=<d>/data_0.parquet of the source workflow")
    ap.add_argument("--structures", type=Path, default=None, help="<export>/structures")
    ap.add_argument("--classifier", type=Path, required=True, help="<project>/classifiers/pixel_classifiers/Islets_<tag>.json")
    ap.add_argument("--report", type=Path, required=True)
    ap.add_argument("--mask", type=Path, default=None, help="also save the classifier's mask (npy, the classifier's grid)")
    ap.add_argument("--objects", type=Path, default=None, help="also save the emulated script-1 objects (npy label image)")
    ap.add_argument("--sigma-um", type=float, default=SIGMA_UM)
    ap.add_argument("--seed-px", type=int, default=SEED_PX)
    ap.add_argument("--close-px", type=int, default=CLOSE_PX)
    ap.add_argument("--open-px", type=int, default=OPEN_PX)
    ap.add_argument("--margin-px", type=int, default=MARGIN_PX)
    ap.add_argument("--miss-weight", type=float, default=MISS_WEIGHT)
    ap.add_argument("--senior-typing", type=Path, default=ict.SENIOR_TYPING)
    args = ap.parse_args()
    t0 = time.time()
    chain = {"seed_px": args.seed_px, "close_px": args.close_px, "open_px": args.open_px, "margin_px": args.margin_px}

    base, values, x_um, y_um, eligible, in_structure, other_markers, quantile, not_in_panel = load_cells(
        args.cells, args.redsea, args.structures, args.tag, args.senior_typing)
    lab = label_cells(values, x_um, y_um, eligible, in_structure, other_markers, quantile)
    planes = read_planes(args.image)
    h, w = next(iter(planes.values())).shape
    smoothed = {c: smooth(planes[c], args.sigma_um / CLASSIFIER_UM) for c in CHANNELS}
    del planes
    cx = np.clip((x_um / CLASSIFIER_UM).astype(int), 0, w - 1); cy = np.clip((y_um / CLASSIFIER_UM).astype(int), 0, h - 1)
    direct, seed = centroid_values(smoothed, cx, cy, **chain)
    t, fit_loss = fit_thresholds(direct, seed, lab["endocrine"], lab["other"], args.miss_weight)
    thresholds = {c: (float(t[j]) if np.isfinite(t[j]) else None) for j, c in enumerate(CHANNELS)}
    mask = classifier_mask(smoothed, thresholds, **chain)
    objects = script1_objects(mask)
    result = score(objects, cx, cy, lab, eligible)

    args.classifier.parent.mkdir(parents=True, exist_ok=True)
    args.classifier.write_text(json.dumps(classifier_json(thresholds, sigma_um=args.sigma_um, **chain), indent=2))
    report = {
        "tag": args.tag, "classifier": str(args.classifier), "image": str(args.image), "cells_file": str(args.cells),
        "redsea_file": str(args.redsea), "thresholds": thresholds,
        "channels_not_used": [c for c in CHANNELS if thresholds[c] is None],
        "sigma_um": args.sigma_um, **chain, "miss_weight": args.miss_weight, "resolution_um": CLASSIFIER_UM,
        "cd99_reference_level": lab["cd99_reference_level"], "cd99_confirmed_level": lab["cd99_confirmed_level"],
        "reference_quantile": quantile, "markers_not_in_this_panel": not_in_panel,
        "cells": {"with_nucleus_in_tissue": int(eligible.sum()), "reference_population": int(lab["reference"].sum()),
                  "hormone_highest": int(lab["hormone_top"].sum()), "confirmed": int(lab["confirmed"].sum()),
                  "endocrine_by_cd99": int(lab["by_cd99"].sum()), "endocrine": int(lab["endocrine"].sum()),
                  "other": int(lab["other"].sum()), "other_enclosed": int(lab["enclosed"].sum()), "left_out": int(lab["left_out"].sum())},
        "fit_loss": fit_loss, "after_script_1": result, "seconds": round(time.time() - t0, 1),
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2))
    if args.mask is not None:
        np.save(args.mask, mask)
    if args.objects is not None:
        np.save(args.objects, objects.astype(np.int32))
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
