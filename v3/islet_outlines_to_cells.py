#!/usr/bin/env python3
"""islet_outlines_to_cells.py -- what qupath/1b_islet_outlines_to_cells.groovy does, in Python, for checking it.

The maintainer, 2026-10-02, after looking at the per-image classifier's outlines:

    "The islet outlines overall look very good, but make the following changes: All outlines need to erode while
    stopping at nearest cell segmentation outline. Distinction between hormone single-cell and new class "Hormone
    Fragment" by cell segmentation outline: if one whole cell fits in islet outline, use cell segmentation; if not,
    label as hormone fragment, the object erodes to marker outline and markers do not contribute to cell type."

In steps (geometry only -- no marker of any cell is read):

  1. A detection FITS an islet outline when no part of it reaches farther than FIT_TOLERANCE_UM beyond the outline
     (one pixel of the islet classifier: the outline is drawn on a 1.0156 um grid, the cells on a 0.5 um one) and
     most of it lies inside. Every other detection the outline touches is CUT BY it.
  2. The outline erodes to the cells: the detections it cuts are taken out of it, so its edge now runs along their
     outlines; the strips left between two neighbouring cut cells (one or two pixels wide) are removed; the
     detections that fit are kept whole.
  3. What is left is split into its connected pieces. A piece holding two or more whole cells (a cell = a detection
     with a nucleus) is an Islet with that piece as its outline. A piece holding exactly one is an Islet whose
     outline IS that cell's segmentation outline. A piece holding none is dropped.
  4. An outline in which no whole cell fits is a "Hormone Fragment": its outline becomes the hormone marker's own
     outline inside it (the image's IsletHormone_<tag> classifier: a hormone above its threshold, nothing else), and
     it is never an islet -- so the hormone in it makes no cell endocrine.

The Groovy step is what runs; this module exists so that its result can be checked against an independent
implementation (islet_classifier_check.py) and tested on synthetic shapes. It reads; it writes nothing into QuPath.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import shapely
from shapely import STRtree

FULL_RES_UM = 0.5078
FIT_TOLERANCE_UM = 1.0156        # one pixel of the islet classifier
# Strips narrower than twice this are removed: the gaps left between two neighbouring cut cells. InstanSeg's outlines
# step in its own pixels (0.985 of an image pixel), and neighbouring cells are one or two of them apart, so those
# strips are 0.985 or 1.97 px wide. 1.25 px removes both and sits half a pixel from either a 2- or a 3-pixel strip,
# so QuPath (JTS) and this module (GEOS) cannot come out on different sides of it. (2.0 px did: 15 of 2,407 islets
# on HDL115 were split by one and not the other.)
WHISKER_UM = 0.63475             # 1.25 full-resolution pixels
CLASS_ISLET, CLASS_FRAGMENT = "Islet", "Hormone Fragment"


def polygons_of(geom) -> list:
    """The polygons of a geometry of any kind (lines and points left by an overlay are not area)."""
    if geom is None or geom.is_empty:
        return []
    if geom.geom_type == "Polygon":
        return [geom]
    if geom.geom_type in ("MultiPolygon", "GeometryCollection"):
        return [p for g in geom.geoms for p in polygons_of(g)]
    return []


def area_only(geom):
    """The polygons of an overlay's result as one geometry: two areas that also share a stretch of boundary intersect
    in polygons PLUS lines, and such a mixed collection cannot be the input of the next overlay."""
    if geom.is_empty or geom.geom_type in ("Polygon", "MultiPolygon"):
        return geom
    parts = polygons_of(geom)
    return shapely.MultiPolygon(parts) if parts else shapely.Polygon()


def opened(geom, width_px: float):
    """Remove what is narrower than 2 x width: shrink by it and grow back, with round joins.

    Round, not mitre: with mitre joins QuPath's geometry library (JTS 1.20) returns almost nothing for the negative
    buffer of some multi-part outlines, so the Groovy step uses round joins and this module does what it does.
    """
    if width_px <= 0 or geom.is_empty:
        return geom
    # an opening never adds anything; the arcs of the round joins are polygons, though, and can reach a few
    # thousandths of a pixel beyond the outline, so the result is cut back to it
    return area_only(shapely.intersection(shapely.buffer(shapely.buffer(geom, -width_px), width_px), geom))


def fits(objs: np.ndarray, cells: np.ndarray, tol_px: float):
    """For every detection: the outline it fits in (-1 = none), and every (detection, outline) pair that overlaps.

    A detection fits ONE outline at most: where two outlines both hold it, the one covering more of it.
    """
    home = np.full(len(cells), -1, dtype=np.int64)
    if len(objs) == 0 or len(cells) == 0:
        return home, np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64)
    ci, oi = STRtree(objs).query(cells, predicate="intersects")
    inside = shapely.area(shapely.intersection(objs[oi], cells[ci]))
    overlap = inside > 0                                    # touching along an edge is not overlapping
    ci, oi, inside = ci[overlap], oi[overlap], inside[overlap]
    grown = np.empty(len(objs), dtype=object)
    used = np.unique(oi)
    grown[used] = shapely.buffer(objs[used], tol_px) if tol_px > 0 else objs[used]
    ok = shapely.covers(grown[oi], cells[ci]) & (inside >= 0.5 * shapely.area(cells[ci]))
    order = np.argsort(inside[ok], kind="stable")           # the largest overlap is written last
    home[ci[ok][order]] = oi[ok][order]
    return home, ci, oi


def erode_to_cells(objs, cells, has_nucleus, *, tol_um: float = FIT_TOLERANCE_UM, whisker_um: float = WHISKER_UM,
                   px_um: float = FULL_RES_UM) -> dict:
    """The eroded outlines. Returns
         islets     list of {geometry, source (index into objs), cells (indices of the whole cells), detections}
         fragments  indices of the outlines no whole cell fits in (their new outline is the marker's; not made here)
         home       per detection, the outline it fits in (-1 = none)
         islet_of   per detection, the index into `islets` of the islet that holds it (-1 = none)
    """
    objs = np.asarray(objs, dtype=object); cells = np.asarray(cells, dtype=object); has_nucleus = np.asarray(has_nucleus, dtype=bool)
    home, ci, oi = fits(objs, cells, tol_um / px_um)
    whisker_px = whisker_um / px_um
    islets, fragments = [], []
    islet_of = np.full(len(cells), -1, dtype=np.int64)
    pairs_by_obj = {}
    for c, o in zip(ci, oi):
        pairs_by_obj.setdefault(int(o), []).append(int(c))
    for o in range(len(objs)):
        touching = np.array(pairs_by_obj.get(o, []), dtype=np.int64)
        whole = touching[home[touching] == o] if len(touching) else touching
        if not has_nucleus[whole].any():
            fragments.append(o)
            continue
        cut = touching[home[touching] != o]
        region = objs[o]
        if len(cut):
            region = area_only(shapely.difference(region, shapely.union_all(cells[cut])))
        region = opened(region, whisker_px)
        region = shapely.union(region, shapely.union_all(cells[whole]))
        pieces = polygons_of(region)
        points = shapely.point_on_surface(cells[whole])
        kept = False
        for piece in pieces:
            members = whole[shapely.contains(piece, points)]
            nucleated = members[has_nucleus[members]]
            if len(nucleated) == 0:
                continue
            geometry = cells[nucleated[0]] if len(nucleated) == 1 else piece
            if len(nucleated) == 1:
                members = nucleated
            islet_of[members] = len(islets)
            islets.append({"geometry": geometry, "source": o, "cells": nucleated, "detections": members})
            kept = True
        if not kept:                                        # cannot happen while a whole cell lies in its own outline
            fragments.append(o)
    return {"islets": islets, "fragments": fragments, "home": home, "islet_of": islet_of}


def cells_near(geojson: Path, objs, margin_px: float = 64.0):
    """Detections of a QuPath cell export that lie near an outline: ids, cell polygons, whether each has a nucleus.

    The export is one line of several GB; it is streamed, and a detection far from every outline is never built.
    """
    import ijson  # noqa: PLC0415
    from shapely.geometry import shape  # noqa: PLC0415
    grid = margin_px
    bounds = shapely.bounds(np.asarray(objs, dtype=object))
    width = int(bounds[:, 2].max() // grid) + 3; height = int(bounds[:, 3].max() // grid) + 3
    near = np.zeros((height, width), dtype=bool)
    for x0, y0, x1, y1 in bounds:
        near[max(int(y0 // grid) - 1, 0):int(y1 // grid) + 2, max(int(x0 // grid) - 1, 0):int(x1 // grid) + 2] = True
    ids, polys, nucleus = [], [], []
    with open(geojson, "rb") as fh:
        for f in ijson.items(fh, "features.item", use_float=True):
            g = f["geometry"]
            ring = g["coordinates"][0] if g["type"] == "Polygon" else g["coordinates"][0][0]
            x, y = ring[0]
            gy, gx = int(y // grid), int(x // grid)
            if gy < 0 or gx < 0 or gy >= height or gx >= width or not near[gy, gx]:
                continue
            ids.append(f["id"]); polys.append(shape(g)); nucleus.append(f.get("nucleusGeometry") is not None)
    polys = np.array(polys, dtype=object)
    bad = ~shapely.is_valid(polys)
    if bad.any():
        polys[bad] = shapely.make_valid(polys[bad])
    return np.array(ids), polys, np.array(nucleus, dtype=bool)
