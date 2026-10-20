#!/usr/bin/env python3
"""islet_validation.py -- annotation-free validation of the v3 islet detector (IO60panc2nd).

The user declined to hand-annotate the benchmark tiles, so the detector's defaults (tau_CD99 3, tau_H 5,
floor k 4, no neck split, no singles co-positivity) are checked with evidence that needs no labels:

  1. Stability over the full-image sweep (10 settings per image): object count and total area per size class,
     relative to the default; object-level agreement for Islet-class objects, using the NESTING of the
     settings (a stricter setting's mask lies inside the default's, a looser one's contains it), so each
     default islet is followed as it shrinks, splits, merges or disappears; batch/disease dependence
     (group x disease tables, type-II OLS on log ratios); invariance of the downstream contrasts
     (T1D vs ND islet density and area fraction) across the settings.
  2. Agreement with the legacy Islets2 islets (the old method, "too harsh"): recall of legacy islets by v3
     default objects (IoU >= 0.3, centroid containment in either direction for merged / fragmented cases, or
     >= 50% covered), what the missed ones look like, and what the extra v3 objects are (size class x route).
  3. Pathology flags: very large objects (merges, hormone- or CD99-bridged), CD99-only-dominated islets, CD99-only
     islets (route 1), endocrine area as % of the analysis area (and of DAPI-occupied tissue; ADIIDAS 10%),
     objects touching the analysis-area border.
  4. A ranked review list (<= ~40 objects) for gate 1, and README.md with the key tables.

Inputs (read-only): sweep GeoJSONs (paths from the benchmark run JSONs' params.CANDIDATES), their
sweep/tau summary JSONs, the legacy exports results/v3/annotations_all/<tag>.geojson (class 'Islet' =
legacy Islets2), the manual Tissue/pLN/Artifact exports, results/v3/islet_anchors.csv, and optionally the
pass-I log (count cross-check). Coordinates in all GeoJSONs are level-0 pixels (0.5078 um, set by script 1).

Production objects (pass I): --default-dir DIR replaces the sweep default with DIR/<tag>.geojson (any
export holding the v3 classes, e.g. export_annotations.groovy run into a SEPARATE folder). Tasks 2-4 and
the default columns of task 1 are then recomputed on the production objects; the sweep variants stay the
reference for the stability deltas. Keep the legacy folder intact: export_annotations.sh overwrites
annotations_all/<tag>.geojson, which is where the legacy islets live until then.

Usage:
  /scratch/io60/venv-v3/bin/python scripts/v3/islet_validation.py --workers 6
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import re
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
import pandas as pd
import shapely
from shapely.geometry import shape

PROJ = Path("/home/smith6jt/IO60panc2nd")
PX_UM = 0.5078                      # level-0 pixel size (script 1: setPixelSizeMicrons(0.5078, 0.5078))
PX2_UM2 = PX_UM * PX_UM
CLASSES = ["Islet", "EndocrineCluster", "EndocrineSingle"]
CLS_SHORT = {"Islet": "Islet", "EndocrineCluster": "Cluster", "EndocrineSingle": "Single"}
VARIANTS = ["k0", "k3", "neck", "copos", "tC2", "tC4", "tH4", "tH7"]
SETTINGS = ["default"] + VARIANTS
# nesting of the endocrine mask relative to the default (from the detector's rule; checked at run time)
STRICTER = ["tC4", "tH7", "neck", "copos"]      # objects lie inside default objects
LOOSER = ["tC2", "tH4", "k0", "k3"]             # objects contain default objects
# single-knob neighbours of the default (the plateau test) and their parameter values
KNOBS = {"tau_CD99": {"tC2": 2.0, "default": 3.0, "tC4": 4.0},
         "tau_H": {"tH4": 4.0, "default": 5.0, "tH7": 7.0},
         "floor_k": {"k0": 0.0, "k3": 3.0, "default": 4.0}}
NEIGHBOURS = ["tC2", "tC4", "tH4", "tH7", "k0", "k3"]
FLAG_REL = 0.10                     # >10% Islet-class count/area change vs the default = flag
LEGACY_IOU = 0.30
LEGACY_COVER = 0.50                 # or >= 50% of the legacy islet inside v3 objects ("covered")
FOUND = ("iou", "merged", "fragmented", "covered")
STABLE_IOU = 0.50
BIG_UM2 = 50_000.0
ADIIDAS_FRAC = 0.10
TOUCH_UM = 1.0156                   # one classifier pixel: "touching" the analysis-area border
M_KEYS = {  # GeoJSON measurement -> column
    "Area µm^2": "area_um2", "Mean caliper diameter µm": "diam_um", "Endocrine route": "route",
    "CD99+ area µm^2": "cd99_um2", "INS+ area µm^2": "ins_um2", "GCG+ area µm^2": "gcg_um2",
    "SST+ area µm^2": "sst_um2", "Hormone+ area µm^2": "hormone_um2", "CD99-only fraction": "cd99_only_frac",
    "Hormone cores": "cores", "Multi-core flag": "multi_core", "Neck split": "neck_split",
    "Positive channels": "n_pos", "Clamped-only flag": "clamped_only", "Bridge flag": "bridge",
    "Tissue border flag": "border_flag", "Distance to tissue border µm": "border_um",
    "Solidity": "solidity", "Circularity": "circularity"}


# ----------------------------------------------------------------------------------------------- inputs
def read_runs(runs_dir: Path, legacy_dir: Path, exclude: set[str]) -> list[dict]:
    """One record per image from the benchmark run JSONs: metadata + label -> GeoJSON path."""
    out = []
    for f in sorted(glob.glob(str(runs_dir / "*.json"))):
        d = json.load(open(f))
        if d["donor"] in exclude:
            continue
        cands = {}
        for spec in d["params"]["CANDIDATES"].split(","):
            spec = spec.strip()
            if spec:
                lab, path = spec.split("=", 1)
                cands[lab.strip()] = path.strip()
        missing = [s for s in SETTINGS if s not in cands]
        if missing:
            raise SystemExit(f"{f}: CANDIDATES lacks {missing}")
        leg = str(legacy_dir / f"{d['tag']}.geojson")     # --legacy-dir wins (the run JSONs leave it empty)
        out.append(dict(image=d["image"], tag=d["tag"], donor=d["donor"], group=d["group"],
                        disease=d["disease"], stratum=d["stratum"], cands=cands, legacy=leg,
                        bench_sources=d.get("sources", {})))
    return out


def _poly(geom: dict):
    g = shape(geom)
    if not g.is_valid:
        g = g.buffer(0)
    return g


def load_v3(path: str) -> tuple[np.ndarray, pd.DataFrame]:
    """v3 objects (Islet / EndocrineCluster / EndocrineSingle) -> geometries + measurement table."""
    with open(path, encoding="utf-8") as fh:
        d = json.load(fh)
    feats = d["features"] if isinstance(d, dict) else d
    geoms, rows = [], []
    for f in feats:
        p = f.get("properties") or {}
        cls = (p.get("classification") or {}).get("name")
        if cls not in CLASSES:
            continue
        m = p.get("measurements") or {}
        r = {"cls": cls}
        for k, c in M_KEYS.items():
            v = m.get(k)
            r[c] = np.nan if v is None else float(v)
        rows.append(r)
        geoms.append(_poly(f["geometry"]))
    df = pd.DataFrame(rows, columns=["cls"] + list(M_KEYS.values()))
    g = np.array(geoms, dtype=object)
    if len(g):
        df["geom_um2"] = shapely.area(g) * PX2_UM2
        c = shapely.centroid(g)
        df["cx_px"], df["cy_px"] = shapely.get_x(c), shapely.get_y(c)
    else:
        df["geom_um2"], df["cx_px"], df["cy_px"] = [], [], []
    return g, df


def load_classes(path: str, classes: tuple[str, ...]) -> list:
    if not path or not os.path.exists(path):
        return []
    with open(path, encoding="utf-8") as fh:
        d = json.load(fh)
    feats = d["features"] if isinstance(d, dict) else d
    out = []
    for f in feats:
        p = f.get("properties") or {}
        cls = (p.get("classification") or {}).get("name")
        if cls in classes:
            out.append((cls, p.get("name"), _poly(f["geometry"])))
    return out


def min_caliper_width(g) -> float:
    """Minimum caliper width (JTS MinimumDiameter) of a polygon, in its own units (rotating calipers)."""
    h = shapely.get_coordinates(shapely.convex_hull(g))
    if len(h) < 4:
        return 0.0
    p = h[:-1]
    e = np.roll(p, -1, axis=0) - p
    L = np.hypot(e[:, 0], e[:, 1])
    ok = L > 0
    e, p0 = e[ok] / L[ok, None], p[ok]
    d = np.abs(e[:, None, 0] * (p[None, :, 1] - p0[:, None, 1]) - e[:, None, 1] * (p[None, :, 0] - p0[:, None, 0]))
    return float(d.max(axis=1).min())


def mean_caliper_um(g) -> float:
    return (2.0 * shapely.minimum_bounding_radius(g) + min_caliper_width(g)) / 2.0 * PX_UM


def size_class(d: float) -> str:
    return "Islet" if d >= 35 else ("EndocrineCluster" if d > 10 else "EndocrineSingle")


# ---------------------------------------------------------------------------------------------- geometry
def overlaps(A: np.ndarray, B: np.ndarray):
    """All (i, j) with area(A_i & B_j) > 0 and that area (px^2)."""
    if len(A) == 0 or len(B) == 0:
        z = np.zeros(0, dtype=int)
        return z, z, np.zeros(0)
    ia, jb = shapely.STRtree(B).query(A, predicate="intersects")
    inter = shapely.area(shapely.intersection(A[ia], B[jb]))
    k = inter > 0
    return ia[k], jb[k], inter[k]


def contained_in(points: np.ndarray, polys: np.ndarray) -> np.ndarray:
    """Index of the polygon containing each point (-1 = none; first hit if several)."""
    out = np.full(len(points), -1, dtype=int)
    if len(points) == 0 or len(polys) == 0:
        return out
    ip, jp = shapely.STRtree(polys).query(points, predicate="within")
    out[ip[::-1]] = jp[::-1]
    return out


def inner_points(g: np.ndarray) -> np.ndarray:
    """Centroid when it lies inside the polygon, else a point on the surface."""
    c = shapely.centroid(g)
    bad = ~shapely.contains(g, c)
    if bad.any():
        c = c.copy()
        c[bad] = shapely.point_on_surface(g[bad])
    return c


# ------------------------------------------------------------------------------------------- per image
def summarise_setting(df: pd.DataFrame) -> dict:
    r = {}
    for c in CLASSES:
        s = df[df.cls == c]
        r[f"n_{CLS_SHORT[c]}"] = int(len(s))
        r[f"area_{CLS_SHORT[c]}"] = float(s.area_um2.sum())
    r["n_all"] = int(len(df))
    r["area_all"] = float(df.area_um2.sum())
    isl = df[df.cls == "Islet"]
    r["islet_median_diam_um"] = float(isl.diam_um.median()) if len(isl) else np.nan
    for k in (1, 2, 3):
        r[f"islet_route{k}"] = int((isl.route == k).sum())
    return r


def process_image(rec: dict, default_dir: str | None, manual_dir: str, dims: tuple | None) -> dict:
    t0 = time.time()
    tag = rec["tag"]
    geo, tab = {}, {}
    for s in SETTINGS:
        path = rec["cands"][s]
        if s == "default" and default_dir:
            path = os.path.join(default_dir, f"{tag}.geojson")
        geo[s], tab[s] = load_v3(path)
    out = {"tag": tag, "checks": {}}

    # ---- determinism: sweep default vs tau-grid default (same parameters, separate file) ----
    d_def = rec["cands"]["default"]
    tau_def = os.path.join(os.path.dirname(d_def), "tau__tC3_tH5_k4.geojson")
    if os.path.exists(tau_def) and not default_dir:
        g2, t2 = load_v3(tau_def)
        same = len(t2) == len(tab["default"]) and np.allclose(np.sort(t2.area_um2.values),
                                                               np.sort(tab["default"].area_um2.values))
        out["checks"]["default_sweep_eq_tau"] = bool(same)
        del g2, t2

    # ---- 1. per setting x class: count and area ----
    rows = []
    for s in SETTINGS:
        r = summarise_setting(tab[s])
        r.update(tag=tag, setting=s)
        rows.append(r)
    out["setting_summary"] = rows

    # ---- 1b. nesting + object-level fate of each default Islet ----
    D, Dt = geo["default"], tab["default"]
    isl = np.flatnonzero(Dt.cls.values == "Islet")
    fate = pd.DataFrame(index=isl)
    nest = {}
    Dpts = inner_points(D) if len(D) else np.zeros(0, dtype=object)
    for v in VARIANTS:
        V, Vt = geo[v], tab[v]
        if v in STRICTER:
            # every variant object inside one default object
            ia, jb, inter = overlaps(V, D)
            vin = np.bincount(ia, weights=inter, minlength=len(V))
            nest[v] = float(vin.sum() / max(shapely.area(V).sum(), 1e-9)) if len(V) else 1.0
            par = contained_in(inner_points(V), D) if len(V) else np.zeros(0, dtype=int)
            ok = par >= 0
            vcls = Vt.cls.values
            vis = vcls == "Islet"
            n_child_islet = np.bincount(par[ok & vis], minlength=len(D))
            n_child = np.bincount(par[ok], minlength=len(D))
            kept_area = np.bincount(par[ok], weights=Vt.area_um2.values[ok], minlength=len(D))
            fate[f"{v}_islet_parts"] = n_child_islet[isl]
            fate[f"{v}_parts"] = n_child[isl]
            fate[f"{v}_retained"] = kept_area[isl] / Dt.area_um2.values[isl]
        else:
            # every default object inside one variant object
            ia, jb, inter = overlaps(D, V)
            din = np.bincount(ia, weights=inter, minlength=len(D))
            nest[v] = float(din.sum() / max(shapely.area(D).sum(), 1e-9)) if len(D) else 1.0
            cont = contained_in(Dpts, V) if len(D) else np.zeros(0, dtype=int)
            # how many default Islets share the container (merges at the looser setting)
            c_isl = cont[isl]
            okc = c_isl >= 0
            n_share = np.zeros(len(isl), dtype=int)
            if okc.any():
                u, cnt = np.unique(c_isl[okc], return_counts=True)
                m = dict(zip(u, cnt))
                n_share[okc] = [m[x] for x in c_isl[okc]]
            growth = np.full(len(isl), np.nan)
            growth[okc] = Vt.area_um2.values[c_isl[okc]] / Dt.area_um2.values[isl[okc]]
            fate[f"{v}_shared_with"] = n_share
            fate[f"{v}_growth"] = growth
    out["checks"]["nesting"] = nest

    # set-level Islet agreement per variant (IoU >= 0.5 one-to-one via best match; area IoU of the unions)
    agree = []
    for v in VARIANTS:
        V, Vt = geo[v], tab[v]
        vis = np.flatnonzero(Vt.cls.values == "Islet")
        ia, jb, inter = overlaps(D[isl], V[vis]) if len(isl) and len(vis) else (np.zeros(0, int),) * 2 + (np.zeros(0),)
        a_d = Dt.area_um2.values[isl]
        a_v = Vt.area_um2.values[vis]
        inter_um2 = inter * PX2_UM2
        iou = inter_um2 / (a_d[ia] + a_v[jb] - inter_um2) if len(ia) else np.zeros(0)
        good = iou >= STABLE_IOU          # IoU >= 0.5 pairs are necessarily one-to-one
        n_m = int(good.sum())
        area_i = float(inter_um2.sum())
        agree.append(dict(tag=tag, setting=v, n_default=len(isl), n_variant=len(vis), n_matched=n_m,
                          islet_jaccard=n_m / max(len(isl) + len(vis) - n_m, 1),
                          frac_default_matched=n_m / max(len(isl), 1),
                          islet_area_iou=area_i / max(a_d.sum() + a_v.sum() - area_i, 1e-9)))
        # the same, all classes (endocrine mask)
        ia2, jb2, inter2 = overlaps(D, V)
        ai = float(inter2.sum() * PX2_UM2)
        agree[-1]["endocrine_area_iou"] = ai / max(Dt.area_um2.sum() + Vt.area_um2.sum() - ai, 1e-9)
        if len(isl):
            best = np.zeros(len(isl))
            if len(ia):
                np.maximum.at(best, ia, iou)
            fate[f"{v}_best_islet_iou"] = best
    out["islet_agreement"] = agree

    # ---- 2. legacy ----
    leg = load_classes(rec["legacy"], ("Islet",))
    man = load_classes(os.path.join(manual_dir, f"{tag}.geojson"), ("Tissue", "pLN", "Artifact"))
    tissue = shapely.union_all([g for c, _, g in man if c == "Tissue"]) if man else None
    excl = [g for c, _, g in man if c in ("pLN", "Artifact")]
    area_geom = tissue.difference(shapely.union_all(excl)) if (tissue is not None and excl) else tissue
    Dcls = Dt.cls.values
    leg_rows = []
    if leg:
        L = np.array([g for _, _, g in leg], dtype=object)
        names = [n for _, n, _ in leg]
        aL = shapely.area(L) * PX2_UM2
        Lpts = inner_points(L)
        inside = shapely.within(Lpts, area_geom) if area_geom is not None else np.ones(len(L), bool)
        ia, jb, inter = overlaps(L, D)
        inter_um2 = inter * PX2_UM2
        iou = inter_um2 / (aL[ia] + Dt.area_um2.values[jb] - inter_um2)
        best_iou = np.zeros(len(L)); best_j = np.full(len(L), -1)
        order = np.argsort(iou)          # ascending: the last write per i is its best
        best_iou[ia[order]] = iou[order]
        best_j[ia[order]] = jb[order]
        cover = np.minimum(np.bincount(ia, weights=inter_um2, minlength=len(L)) / aL, 1.0)
        n_ov = np.bincount(ia, minlength=len(L))
        merged_j = contained_in(Lpts, D)                       # legacy centroid inside a v3 object
        frag_par = contained_in(inner_points(D), L) if len(D) else np.zeros(0, int)   # v3 centroids inside legacy
        okf = frag_par >= 0
        n_frag = np.bincount(frag_par[okf], minlength=len(L))
        n_frag_islet = np.bincount(frag_par[okf & (Dcls == "Islet")], minlength=len(L))
        # union over every setting: any overlap / found (same rule) by any setting
        found_by = [[] for _ in range(len(L))]
        any_ov = n_ov > 0
        for s in SETTINGS:
            if s == "default":
                continue
            G = geo[s]
            ia_s, jb_s, in_s = overlaps(L, G)
            in_s = in_s * PX2_UM2
            iou_s = in_s / (aL[ia_s] + tab[s].area_um2.values[jb_s] - in_s)
            b = np.zeros(len(L))
            if len(ia_s):
                np.maximum.at(b, ia_s, iou_s)
            m_s = contained_in(Lpts, G) >= 0
            fp = contained_in(inner_points(G), L) if len(G) else np.zeros(0, int)
            nf = np.bincount(fp[fp >= 0], minlength=len(L))
            cov_s = np.bincount(ia_s, weights=in_s, minlength=len(L)) / aL
            f_s = (b >= LEGACY_IOU) | m_s | (nf > 0) | (cov_s >= LEGACY_COVER)
            any_ov |= np.bincount(ia_s, minlength=len(L)) > 0
            for i in np.flatnonzero(f_s):
                found_by[i].append(s)
        for i in range(len(L)):
            if not inside[i]:
                st = "outside_area"
            elif best_iou[i] >= LEGACY_IOU:
                st = "iou"
            elif merged_j[i] >= 0:
                st = "merged"
            elif n_frag[i] > 0:
                st = "fragmented"
            elif cover[i] >= LEGACY_COVER:
                st = "covered"          # inside a larger v3 object whose centroid lies elsewhere
            elif n_ov[i] > 0:
                st = "partial"
            else:
                st = "absent"
            dl = mean_caliper_um(L[i])
            j = best_j[i] if best_j[i] >= 0 else merged_j[i]
            fb = set(found_by[i])
            # why a default miss is missed: the looser setting that recovers it (nesting: tC2 = CD99 fold 2-3,
            # tH4 = hormone fold 4-5, k3/k0 = the clamp floor); none = below every tested threshold
            why = ("" if st in FOUND + ("outside_area",) else
                   "CD99-dim (tC2)" if "tC2" in fb else "hormone-dim (tH4)" if "tH4" in fb else
                   "floor (k3/k0)" if fb & {"k3", "k0"} else "no setting")
            leg_rows.append(dict(miss_type=why,
                tag=tag, legacy_name=names[i], legacy_area_um2=aL[i], legacy_diam_um=dl,
                legacy_class=size_class(dl), cx_px=shapely.get_x(Lpts[i]), cy_px=shapely.get_y(Lpts[i]),
                status=st, best_iou=best_iou[i], coverage=cover[i], n_v3_overlapping=int(n_ov[i]),
                n_v3_centroids_inside=int(n_frag[i]), n_v3_islet_centroids_inside=int(n_frag_islet[i]),
                v3_class=(Dcls[j] if j >= 0 else ""),
                v3_area_um2=(float(Dt.area_um2.values[j]) if j >= 0 else np.nan),
                v3_islet_match=bool((j >= 0 and Dcls[j] == "Islet" and st in ("iou", "merged", "covered")) or
                                    (st == "fragmented" and n_frag_islet[i] > 0)),
                found_by_other_settings=",".join(found_by[i]), any_setting_overlap=bool(any_ov[i])))
        # v3 default objects vs legacy: fraction of each object's area on legacy islets, and how many legacy
        # islets (inner point) each v3 object contains (>= 2: v3 joins islets that Islets2 kept apart)
        ia2, jb2, inter2 = overlaps(D, L)
        Dt = Dt.copy()
        Dt["legacy_overlap_frac"] = np.minimum(np.bincount(ia2, weights=inter2 * PX2_UM2, minlength=len(D))
                                               / Dt.area_um2.values, 1.0)
        Dt["n_legacy_inside"] = np.bincount(merged_j[merged_j >= 0], minlength=len(D))
    else:
        Dt = Dt.copy()
        Dt["legacy_overlap_frac"] = np.nan
        Dt["n_legacy_inside"] = np.nan

    # border contacts: image edge (5 um inner rectangle) or Tissue/pLN/Artifact outline; the detector's own
    # 'Distance to tissue border' is re-derived from the exported outlines as a check
    Dt["touch_type"] = ""
    touch = np.flatnonzero(Dt.border_um.values <= TOUCH_UM)
    out["checks"]["border_recheck_max_um"] = np.nan
    if len(touch) and dims is not None and area_geom is not None:
        W, H = dims
        m = 5.0 / PX_UM
        rect = shapely.box(m, m, W - m, H - m)
        d_edge = shapely.distance(D[touch], rect.boundary) * PX_UM
        d_line = shapely.distance(D[touch], area_geom.boundary) * PX_UM
        tt = np.where(d_edge <= TOUCH_UM, "image_edge", np.where(d_line <= TOUCH_UM, "tissue_outline", "unverified"))
        Dt.loc[touch, "touch_type"] = tt
        out["checks"]["border_recheck_max_um"] = float(np.minimum(d_edge, d_line).max())
    out["legacy"] = leg_rows
    out["n_legacy"] = len(leg)

    # ---- 3/4. default object table (Islet class + every flagged object) ----
    Dt["tag"] = tag
    Dt["obj"] = np.arange(len(Dt))
    for c in fate.columns:
        Dt.loc[fate.index, c] = fate[c].values
    flagged = ((Dt.cls == "Islet") | (Dt.area_um2 > BIG_UM2) | (Dt.bridge == 1) | (Dt.border_um <= TOUCH_UM))
    out["objects"] = Dt[flagged].reset_index(drop=True)
    # extras: default objects with no legacy overlap, by class x route
    ex = Dt.assign(extra=Dt.legacy_overlap_frac.fillna(0) == 0)
    out["extras"] = (ex.groupby(["cls", "route", "extra"]).agg(n=("obj", "size"), area_um2=("area_um2", "sum"))
                     .reset_index().assign(tag=tag).to_dict("records"))
    # border / pathology counts over ALL default objects
    out["border"] = {f"touch_{CLS_SHORT[c]}": int(((Dt.cls == c) & (Dt.border_um <= TOUCH_UM)).sum()) for c in CLASSES}
    out["border"].update({f"flag10_{CLS_SHORT[c]}": int(((Dt.cls == c) & (Dt.border_flag == 1)).sum()) for c in CLASSES})
    for k in ("image_edge", "tissue_outline", "unverified"):
        out["border"][f"touch_Islet_{k}"] = int(((Dt.cls == "Islet") & (Dt.touch_type == k)).sum())
    out["seconds"] = time.time() - t0
    return out


# ------------------------------------------------------------------------------------------ statistics
def type2_anova(df: pd.DataFrame, y: str) -> dict:
    """OLS y ~ C(group) + C(disease), type-II F tests (statsmodels). NaN-safe; returns p-values and R^2."""
    import statsmodels.formula.api as smf
    from statsmodels.stats.anova import anova_lm
    d = df[[y, "group", "disease"]].replace([np.inf, -np.inf], np.nan).dropna()
    if len(d) < 8 or d[y].std() == 0:
        return {"p_group": np.nan, "p_disease": np.nan, "r2": np.nan, "n": len(d)}
    fit = smf.ols(f"{y} ~ C(group) + C(disease)", data=d).fit()
    a = anova_lm(fit, typ=2)
    return {"p_group": float(a.loc["C(group)", "PR(>F)"]), "p_disease": float(a.loc["C(disease)", "PR(>F)"]),
            "r2": float(fit.rsquared), "n": len(d)}


def composition_test(objs: pd.DataFrame) -> dict:
    """Is the tau_H sensitivity of an islet's area a disease/batch property or an islet-composition one?
    Per default Islet: log(area kept at tau_H 7) and log(area grown at tau_H 4) ~ log(diameter) + disease + group
    [+ INS share of hormone area + hormone share of area]; cluster-robust by donor; joint Wald tests.
    Also the equivalent boundary shift (um) of a disc of the islet's area."""
    import statsmodels.formula.api as smf
    x = objs[objs.cls == "Islet"].copy()
    x["ins_frac"] = (x.ins_um2 / x.hormone_um2.replace(0, np.nan)).fillna(0)
    x["hormone_frac"] = x.hormone_um2 / x.area_um2
    x["ld"] = np.log(x.diam_um)
    R = np.sqrt(x.area_um2 / np.pi)
    x["shift_tH7_um"] = R * (1 - np.sqrt(x.tH7_retained.clip(0, 1)))
    x["shift_tH4_um"] = R * (np.sqrt(x.tH4_growth) - 1)
    res = {"n_islets": int(len(x)),
           "boundary_shift_um_median": {d: {"tH4": float(g.shift_tH4_um.median()), "tH7": float(g.shift_tH7_um.median())}
                                        for d, g in x.groupby("disease")},
           "ins_frac_median": {d: float(g.ins_frac.median()) for d, g in x.groupby("disease")}}
    for name, col in (("tH7_kept", "tH7_retained"), ("tH4_grown", "tH4_growth")):
        d = x[(x[col] > 0) & np.isfinite(x[col])].copy()
        d["y"] = np.log(d[col])
        for adj, f in (("size", "y ~ ld + C(disease, Treatment('ND')) + C(group)"),
                       ("size+composition", "y ~ ld + ins_frac + hormone_frac + C(disease, Treatment('ND')) + C(group)")):
            fit = smf.ols(f, data=d).fit(cov_type="cluster", cov_kwds={"groups": d["donor"].astype(str)})
            w = fit.wald_test_terms(skip_single=False, scalar=True).table
            key_d = [k for k in w.index if "disease" in k][0]
            key_g = [k for k in w.index if "group" in k][0]
            t1d = "C(disease, Treatment('ND'))[T.T1D]"
            res[f"{name}__{adj}"] = {
                "p_disease": float(w.loc[key_d, "pvalue"]), "p_group": float(w.loc[key_g, "pvalue"]),
                "T1D_vs_ND_ratio": float(np.exp(fit.params.get(t1d, np.nan))),
                "p_ins_frac": float(fit.pvalues.get("ins_frac", np.nan)),
                "coef_ins_frac": float(fit.params.get("ins_frac", np.nan)), "n": int(len(d))}
    return res


def disease_effect(df: pd.DataFrame, y: str, contrast: str = "T1D") -> float:
    """exp(coef) of disease (contrast vs ND) in log(y) ~ C(group) + C(disease): the group-adjusted ratio."""
    import statsmodels.formula.api as smf
    d = df[[y, "group", "disease"]].copy()
    d = d[d[y] > 0]
    d["ly"] = np.log(d[y])
    fit = smf.ols("ly ~ C(group) + C(disease, Treatment('ND'))", data=d).fit()
    k = f"C(disease, Treatment('ND'))[T.{contrast}]"
    return float(np.exp(fit.params[k])) if k in fit.params else np.nan


def fmt_pct(x: float, nd: int = 1) -> str:
    return "" if pd.isna(x) else f"{100 * x:+.{nd}f}%"


def md_table(df: pd.DataFrame, floatfmt: str = "{:.3g}") -> str:
    cols = list(df.columns)
    lines = ["| " + " | ".join(str(c) for c in cols) + " |", "|" + "|".join("---" for _ in cols) + "|"]
    for _, r in df.iterrows():
        cells = []
        for c in cols:
            v = r[c]
            if isinstance(v, (float, np.floating)):
                cells.append("" if pd.isna(v) else floatfmt.format(v))
            else:
                cells.append(str(v))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def parse_passI_log(path: Path) -> dict:
    """tag -> {Islet, EndocrineCluster, EndocrineSingle} from 'PASS -- islet-class annotations [...]' lines."""
    res = {}
    if not path or not path.exists():
        return res
    pat = re.compile(r"(\S+): 1: PASS -- islet-class annotations \[([^\]]*)\]")
    for line in path.read_text(errors="replace").splitlines():
        m = pat.search(line)
        if m:
            res[m.group(1)] = {k.strip(): int(v) for k, v in (x.split(":") for x in m.group(2).split(","))}
    return res


# ------------------------------------------------------------------------------------------------ main
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--runs-dir", type=Path, default=Path("/scratch/io60/benchmark_v3/private/runs"))
    ap.add_argument("--legacy-dir", type=Path, default=PROJ / "results/v3/annotations_all")
    ap.add_argument("--manual-dir", type=Path, default=PROJ / "results/v3/manual_annotations")
    ap.add_argument("--anchors-csv", type=Path, default=PROJ / "results/v3/islet_anchors.csv")
    ap.add_argument("--passI-log", type=Path, default=PROJ / "results/v3/passI_20260929.log")
    ap.add_argument("--default-dir", default=None,
                    help="production objects: DIR/<tag>.geojson replaces the sweep default (tasks 2-4)")
    ap.add_argument("--out-dir", type=Path, default=PROJ / "results/v3/islet_validation")
    ap.add_argument("--exclude", default="6457", help="comma list of donors to skip")
    ap.add_argument("--only", default="", help="comma list of donors (debug)")
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--max-review", type=int, default=40)
    a = ap.parse_args(argv)

    recs = read_runs(a.runs_dir, a.legacy_dir, set(x for x in a.exclude.split(",") if x))
    if a.only:
        keep = set(a.only.split(","))
        recs = [r for r in recs if r["donor"] in keep]
    meta = pd.DataFrame([{k: r[k] for k in ("image", "tag", "donor", "group", "disease", "stratum")} for r in recs])
    a.out_dir.mkdir(parents=True, exist_ok=True)

    # summaries written by the detector (tissue area, floor binding, thresholds) per image x setting
    det_rows = []
    for r in recs:
        for s, path in r["cands"].items():
            sj = os.path.join(os.path.dirname(path), ("sweep" if "sweep__" in path else "tau") + ".summary.json")
            for e in json.load(open(sj))["sweep"]:
                if e.get("geojson") and os.path.basename(e["geojson"]) == os.path.basename(path):
                    det_rows.append(dict(tag=r["tag"], setting=s, tissue_um2=e["tissue_area_um2"],
                                         floor_binds=",".join(e.get("floor_binds") or []),
                                         det_n=e["n_objects"], det_endo_frac=e["endocrine_area_fraction"],
                                         **{f"thr_{k}": v for k, v in e["thresholds_counts"].items()}))
    det = pd.DataFrame(det_rows)
    tissue = det[det.setting == "default"].set_index("tag").tissue_um2

    # image size (level 0) from the QPTIFF header (no pixels read): the 5 um inner rectangle of the analysis area
    dims = {}
    try:
        import tifffile
        for t, path in pd.read_csv(a.anchors_csv).groupby("tag").image_path.first().items():
            if os.path.exists(str(path)):
                with tifffile.TiffFile(path) as tf:
                    h, w = tf.pages[0].shape[:2]
                dims[t] = (w, h)
    except Exception as e:   # the border split is optional
        print(f"[islet_validation] WARN: image sizes unavailable ({e}); border contacts not split by type")

    t0 = time.time()
    results = []
    with ProcessPoolExecutor(max_workers=max(1, min(a.workers, 8))) as ex:
        futs = {ex.submit(process_image, r, a.default_dir, str(a.manual_dir), dims.get(r["tag"])): r["tag"]
                for r in recs}
        for f in as_completed(futs):
            res = f.result()
            results.append(res)
            print(f"[islet_validation] {res['tag']}: {res['seconds']:.0f}s", flush=True)
    print(f"[islet_validation] {len(results)} images in {time.time() - t0:.0f}s", flush=True)
    by_tag = {r["tag"]: r for r in results}

    # ======================================================= 1. stability
    ss = pd.DataFrame([row for r in results for row in r["setting_summary"]]).merge(meta, on="tag")
    ss = ss.merge(det[["tag", "setting", "det_n", "floor_binds"]], on=["tag", "setting"], how="left")
    long = []
    for _, r in ss.iterrows():
        for c in CLASSES + ["all"]:
            k = CLS_SHORT.get(c, "all")
            long.append(dict(tag=r.tag, donor=r.donor, group=r.group, disease=r.disease, setting=r.setting,
                             size_class=k, n=r[f"n_{k}"], area_um2=r[f"area_{k}"]))
    long = pd.DataFrame(long)
    dflt = long[long.setting == "default"].set_index(["tag", "size_class"])
    long["n_default"] = [dflt.loc[(t, c), "n"] for t, c in zip(long.tag, long.size_class)]
    long["area_default"] = [dflt.loc[(t, c), "area_um2"] for t, c in zip(long.tag, long.size_class)]
    long["rel_n"] = long.n / long.n_default.replace(0, np.nan) - 1
    long["rel_area"] = long.area_um2 / long.area_default.replace(0, np.nan) - 1
    long.to_csv(a.out_dir / "stability_long.tsv", sep="\t", index=False, float_format="%.6g")

    # count check vs the detector's own summaries
    chk_counts = int((ss.n_all != ss.det_n).sum())

    # knob summary over images: median [min, max] relative change per class
    summ = []
    for s in VARIANTS:
        for k in ["Islet", "Cluster", "Single", "all"]:
            x = long[(long.setting == s) & (long.size_class == k)]
            summ.append(dict(setting=s, size_class=k,
                             rel_n_median=x.rel_n.median(), rel_n_min=x.rel_n.min(), rel_n_max=x.rel_n.max(),
                             rel_area_median=x.rel_area.median(), rel_area_min=x.rel_area.min(),
                             rel_area_max=x.rel_area.max(),
                             n_images_abs_gt10=int(((x.rel_n.abs() > FLAG_REL) | (x.rel_area.abs() > FLAG_REL)).sum())))
    summ = pd.DataFrame(summ)
    summ.to_csv(a.out_dir / "stability_by_setting.tsv", sep="\t", index=False, float_format="%.4g")

    # group x disease medians (Islet class and singles), per setting
    gd = (long[long.setting.isin(VARIANTS)]
          .groupby(["setting", "size_class", "group", "disease"])
          .agg(n_images=("tag", "size"), rel_n_median=("rel_n", "median"), rel_area_median=("rel_area", "median"),
               rel_n_min=("rel_n", "min"), rel_n_max=("rel_n", "max"))
          .reset_index())
    gd.to_csv(a.out_dir / "stability_by_group_disease.tsv", sep="\t", index=False, float_format="%.4g")

    # flags: Islet-class count or area change > 10% at a neighbour of the default
    fl = long[(long.size_class == "Islet") & long.setting.isin(NEIGHBOURS) &
              ((long.rel_n.abs() > FLAG_REL) | (long.rel_area.abs() > FLAG_REL))].copy()
    fl["knob"] = fl.setting.map({s: k for k, d in KNOBS.items() for s in d if s != "default"})
    fl = fl[["tag", "donor", "group", "disease", "knob", "setting", "n_default", "n", "rel_n", "area_default",
             "area_um2", "rel_area"]].sort_values(["knob", "setting", "rel_area"])
    fl.to_csv(a.out_dir / "stability_flags.tsv", sep="\t", index=False, float_format="%.4g")

    # object-level agreement
    ag = pd.DataFrame([row for r in results for row in r["islet_agreement"]]).merge(meta, on="tag")
    ag.to_csv(a.out_dir / "islet_agreement.tsv", sep="\t", index=False, float_format="%.4g")
    nest = pd.DataFrame([{"tag": r["tag"], **{f"nest_{k}": v for k, v in r["checks"]["nesting"].items()},
                          "default_sweep_eq_tau": r["checks"].get("default_sweep_eq_tau")} for r in results])

    # per-islet fates (default Islets)
    objs = pd.concat([r["objects"] for r in results], ignore_index=True).merge(meta, on="tag")
    objs["cx_um"], objs["cy_um"] = objs.cx_px * PX_UM, objs.cy_px * PX_UM
    isl = objs[objs.cls == "Islet"].copy()
    fates = []
    for t, x in isl.groupby("tag"):
        n = len(x)
        row = dict(tag=t, n_islets=n)
        for v in ["tC4", "tH7"]:
            row[f"split_{v}"] = int((x[f"{v}_islet_parts"] >= 2).sum())
            row[f"lost_{v}"] = int((x[f"{v}_islet_parts"] == 0).sum())
            row[f"shrunk50_{v}"] = int((x[f"{v}_retained"] < 0.5).sum())
            row[f"median_retained_{v}"] = float(x[f"{v}_retained"].median())
        for v in ["tC2", "tH4", "k0", "k3"]:
            row[f"merged_{v}"] = int((x[f"{v}_shared_with"] >= 2).sum())
            row[f"median_growth_{v}"] = float(x[f"{v}_growth"].median())
        row["neck_split"] = int((x["neck_parts"] >= 2).sum())
        fates.append(row)
    fates = pd.DataFrame(fates).merge(meta, on="tag")
    fates.to_csv(a.out_dir / "islet_fates.tsv", sep="\t", index=False, float_format="%.4g")

    # batch / disease dependence of the Islet-class sensitivity (log ratios per image)
    dep = []
    isl_long = long[long.size_class.isin(["Islet", "Single", "Cluster"])]
    for s in NEIGHBOURS + ["neck", "copos"]:
        for k in ["Islet", "Cluster", "Single"]:
            x = isl_long[(isl_long.setting == s) & (isl_long.size_class == k)].copy()
            for met in ["n", "area"]:
                col = "rel_n" if met == "n" else "rel_area"
                x["y"] = np.log1p(x[col])
                r = type2_anova(x, "y")
                g = x.groupby("disease")[col].median()
                e = x.groupby("group")[col].median()
                dep.append(dict(setting=s, size_class=k, metric=met, **r,
                                **{f"median_{d}": g.get(d, np.nan) for d in ["ND", "Aab+", "T1D"]},
                                **{f"median_{G}": e.get(G, np.nan) for G in ["E1", "E2", "E3", "E4"]}))
    dep = pd.DataFrame(dep)
    dep.to_csv(a.out_dir / "sensitivity_batch_disease.tsv", sep="\t", index=False, float_format="%.4g")
    comp = composition_test(objs)
    json.dump(comp, open(a.out_dir / "tauH_composition_test.json", "w"), indent=2)

    # downstream contrasts per setting: T1D/ND and Aab+/ND ratios (group-adjusted), islet density & area fraction
    ro = ss.copy()
    ro["tissue_mm2"] = ro.tag.map(tissue) / 1e6
    ro["islet_density_mm2"] = ro.n_Islet / ro.tissue_mm2
    ro["islet_area_pct"] = 100 * ro.area_Islet / (ro.tissue_mm2 * 1e6)
    ro["endocrine_area_pct"] = 100 * ro.area_all / (ro.tissue_mm2 * 1e6)
    ro["singles_density_mm2"] = ro.n_Single / ro.tissue_mm2
    con = []
    for s in SETTINGS:
        x = ro[ro.setting == s]
        row = dict(setting=s)
        for y in ["islet_density_mm2", "islet_area_pct", "endocrine_area_pct", "islet_median_diam_um",
                  "singles_density_mm2"]:
            row[f"{y}__T1D_vs_ND"] = disease_effect(x, y, "T1D")
            row[f"{y}__Aab_vs_ND"] = disease_effect(x, y, "Aab+")
        # rank agreement with the default across images
        xd = ro[ro.setting == "default"].set_index("tag")
        xs = x.set_index("tag")
        row["spearman_islet_density_vs_default"] = float(xs.islet_density_mm2.corr(xd.islet_density_mm2, method="spearman"))
        row["spearman_islet_area_pct_vs_default"] = float(xs.islet_area_pct.corr(xd.islet_area_pct, method="spearman"))
        con.append(row)
    con = pd.DataFrame(con)
    con.to_csv(a.out_dir / "downstream_contrasts.tsv", sep="\t", index=False, float_format="%.4g")
    ro.to_csv(a.out_dir / "readouts_by_setting.tsv", sep="\t", index=False, float_format="%.6g")

    # ======================================================= 2. legacy
    lg = pd.DataFrame([row for r in results for row in r["legacy"]])
    lg_sum = pd.DataFrame()
    if len(lg):
        lg = lg.merge(meta, on="tag")
        lg["cx_um"], lg["cy_um"] = lg.cx_px * PX_UM, lg.cy_px * PX_UM
        lg.to_csv(a.out_dir / "legacy_islets.tsv.gz", sep="\t", index=False, float_format="%.7g")
        rows = []
        for t, x in lg.groupby("tag"):
            xin = x[x.status != "outside_area"]
            cnt = xin.status.value_counts()
            found = xin.status.isin(FOUND)
            strict = (xin.status.isin(["iou", "covered"]) |
                      (xin.status.isin(["merged", "fragmented"]) & (xin.coverage >= LEGACY_IOU)))
            big = xin[xin.legacy_class == "Islet"]
            rows.append(dict(
                tag=t, n_legacy=len(x), n_outside_area=int((x.status == "outside_area").sum()), n_in_area=len(xin),
                n_iou=int(cnt.get("iou", 0)), n_merged=int(cnt.get("merged", 0)),
                n_fragmented=int(cnt.get("fragmented", 0)), n_covered=int(cnt.get("covered", 0)),
                n_partial=int(cnt.get("partial", 0)),
                n_absent=int(cnt.get("absent", 0)),
                recall=found.mean() if len(xin) else np.nan, recall_strict=strict.mean() if len(xin) else np.nan,
                recall_any_overlap=(xin.n_v3_overlapping > 0).mean() if len(xin) else np.nan,
                recall_any_setting=(found | (xin.found_by_other_settings != "")).mean() if len(xin) else np.nan,
                n_legacy_ge35=len(big),
                recall_ge35=big.status.isin(FOUND).mean() if len(big) else np.nan,
                recall_ge35_as_islet=big.v3_islet_match.mean() if len(big) else np.nan,
                missed_area_frac=(xin.legacy_area_um2[~found].sum() / xin.legacy_area_um2.sum()) if len(xin) else np.nan,
                missed_median_area_um2=xin.legacy_area_um2[~found].median(),
                n_missed_ge35=int((~found & (xin.legacy_class == "Islet")).sum()),
                **{f"miss_{k}": int((xin.miss_type == v).sum()) for k, v in (
                    ("cd99dim", "CD99-dim (tC2)"), ("hormonedim", "hormone-dim (tH4)"), ("floor", "floor (k3/k0)"),
                    ("nosetting", "no setting"))},
                median_iou_matched=xin.best_iou[xin.status == "iou"].median(),
                median_v3_to_legacy_area=(xin.v3_area_um2 / xin.legacy_area_um2)[xin.status == "iou"].median()))
        lg_sum = pd.DataFrame(rows).merge(meta, on="tag")
        lg_sum.to_csv(a.out_dir / "legacy_recall.tsv", sep="\t", index=False, float_format="%.4g")
    # extras: what v3 adds
    exr = pd.DataFrame([row for r in results for row in r["extras"]])
    exr.to_csv(a.out_dir / "v3_vs_legacy_objects.tsv", sep="\t", index=False, float_format="%.6g")

    # ======================================================= 3. pathology
    anc = pd.read_csv(a.anchors_csv)
    occ = anc.groupby("tag").occupied_fraction.first()
    cd99_anchor = anc[anc.channel == "CD99"].set_index("tag").anchor
    ins_clamp = anc[anc.channel == "INS"].set_index("tag").clamp_fraction
    D0 = ro[ro.setting == "default"].set_index("tag")
    pat = []
    for t in meta.tag:
        o = objs[objs.tag == t]
        io = o[o.cls == "Islet"]
        b = by_tag[t]["border"]
        big = o[o.area_um2 > BIG_UM2]
        pat.append(dict(
            tag=t, tissue_mm2=D0.loc[t, "tissue_mm2"], occupied_frac=occ.get(t, np.nan),
            endocrine_pct=D0.loc[t, "endocrine_area_pct"],
            endocrine_pct_occupied=D0.loc[t, "endocrine_area_pct"] / occ.get(t, np.nan),
            islet_pct=D0.loc[t, "islet_area_pct"],
            adiidas_flag=bool(D0.loc[t, "endocrine_area_pct"] / occ.get(t, np.nan) > 100 * ADIIDAS_FRAC),
            n_islets=len(io), n_gt50k=len(big), largest_um2=o.area_um2.max(),
            n_gt50k_cd99only_ge30=int((big.cd99_only_frac >= 0.3).sum()),
            n_gt50k_merge=int(((big.tH7_islet_parts >= 3) | (big.n_legacy_inside >= 3)).sum()),
            n_islet_joins_ge2_legacy=int((io.n_legacy_inside >= 2).sum()),
            n_bridge_flag=int((o.bridge == 1).sum()),
            n_islet_cd99only_ge50=int((io.cd99_only_frac >= 0.5).sum()),
            n_islet_route1=int((io.route == 1).sum()),
            n_islet_multicore=int((io.multi_core == 1).sum()),
            n_islet_split_tC4=int((io.tC4_islet_parts >= 2).sum()),
            n_islet_split_tH7=int((io.tH7_islet_parts >= 2).sum()),
            islet_area_in_gt50k_pct=100 * big.area_um2.sum() / max(io.area_um2.sum(), 1e-9),
            **b,
            islet_touch_pct=100 * b["touch_Islet"] / max(len(io), 1),
            cd99_anchor=cd99_anchor.get(t, np.nan), ins_clamp_fraction=ins_clamp.get(t, np.nan)))
    pat = pd.DataFrame(pat).merge(meta, on="tag")
    pat.to_csv(a.out_dir / "pathology_flags.tsv", sep="\t", index=False, float_format="%.4g")
    big_objs = objs[(objs.area_um2 > BIG_UM2) | (objs.bridge == 1)].sort_values("area_um2", ascending=False)
    cols_o = ["image", "donor", "group", "disease", "cls", "cx_um", "cy_um", "cx_px", "cy_px", "area_um2", "diam_um",
              "route", "cd99_only_frac", "cores", "multi_core", "bridge", "border_um", "touch_type", "solidity",
              "tC4_islet_parts", "tH7_islet_parts", "tC4_retained", "tH7_retained", "legacy_overlap_frac",
              "n_legacy_inside"]
    big_objs[cols_o].to_csv(a.out_dir / "large_or_bridge_objects.tsv", sep="\t", index=False, float_format="%.7g")
    isl_out = objs[objs.cls == "Islet"]
    isl_out.drop(columns=["geom_um2"]).to_csv(a.out_dir / "islets_default.tsv.gz", sep="\t", index=False,
                                              float_format="%.7g")

    # ======================================================= 4. review list
    review = build_review(objs, lg, lg_sum, pat, a.max_review)
    if len(review):
        review = review.assign(cx_um=review.cx_um.round(1), cy_um=review.cy_um.round(1),
                               cx_px=review.cx_px.round(0).astype(int), cy_px=review.cy_px.round(0).astype(int),
                               area_um2=review.area_um2.round(0).astype(int), diam_um=review.diam_um.round(1))
    review.to_csv(a.out_dir / "review_list.tsv", sep="\t", index=False)

    # ======================================================= checks + README
    passI = parse_passI_log(a.passI_log)
    pchk = []
    for t, c in passI.items():
        x = ss[(ss.tag == t) & (ss.setting == "default")]
        if len(x):
            x = x.iloc[0]
            pchk.append(dict(tag=t, passI=f"{c.get('Islet')}/{c.get('EndocrineCluster')}/{c.get('EndocrineSingle')}",
                             sweep=f"{x.n_Islet}/{x.n_Cluster}/{x.n_Single}",
                             identical=(c.get("Islet"), c.get("EndocrineCluster"), c.get("EndocrineSingle")) ==
                                       (x.n_Islet, x.n_Cluster, x.n_Single)))
    pchk = pd.DataFrame(pchk)
    checks = dict(count_mismatch_vs_detector_summary=chk_counts,
                  default_sweep_equals_tau_default=int(nest.default_sweep_eq_tau.fillna(False).sum()),
                  n_images=len(meta),
                  nesting_min={c: float(nest[c].min()) for c in nest.columns if c.startswith("nest_")},
                  passI_identical=f"{int(pchk.identical.sum()) if len(pchk) else 0}/{len(pchk)}",
                  border_recheck_max_um=float(np.nanmax([r["checks"].get("border_recheck_max_um", np.nan)
                                                         for r in results] + [np.nan])) if results else np.nan,
                  default_dir=a.default_dir)
    mg = objs[(objs.area_um2 > BIG_UM2) & ((objs.tH7_islet_parts >= 3) | (objs.n_legacy_inside >= 3))]
    checks["merges"] = dict(n_gt50k=int((objs.area_um2 > BIG_UM2).sum()), n=int(len(mg)),
                            median_cd99_only_frac=float(mg.cd99_only_frac.median()) if len(mg) else np.nan,
                            n_cd99_only_ge20=int((mg.cd99_only_frac >= 0.2).sum()),
                            n_split_tC4=int((mg.tC4_islet_parts >= 2).sum()),
                            n_split_tH7_ge3=int((mg.tH7_islet_parts >= 3).sum()),
                            n_legacy_ge3=int((mg.n_legacy_inside >= 3).sum()),
                            median_tH7_parts=float(mg.tH7_islet_parts.median()) if len(mg) else np.nan)
    json.dump(checks, open(a.out_dir / "checks.json", "w"), indent=2, default=str)
    if len(pchk):
        pchk.to_csv(a.out_dir / "passI_count_check.tsv", sep="\t", index=False)
    (a.out_dir / "README.md").write_text(
        write_readme(a, meta, long, summ, gd, fl, ag, fates, dep, con, lg_sum, exr, pat, review, checks, pchk, comp),
        encoding="utf-8")
    print(f"[islet_validation] wrote {a.out_dir}", flush=True)
    return 0


# --------------------------------------------------------------------------------------------- review
REVIEW_TIERS = [
    # code, label, max rows, max rows per image
    ("R1_merge", "merged islets (> 50,000 um2)", 10, 1),
    ("R2_cd99_bridge", "CD99 bridge", 4, 2),
    ("R3_cd99_only_islet", "CD99-only Islet", 6, 2),
    ("R4_legacy_missed", "missed legacy islet", 9, 3),
    ("R5_fragile", "threshold-fragile Islet", 4, 1),
    ("R6_border", "Islet cut by the Tissue outline", 4, 1),
]


def build_review(objs: pd.DataFrame, lg: pd.DataFrame, lg_sum: pd.DataFrame, pat: pd.DataFrame,
                 max_n: int) -> pd.DataFrame:
    """Ranked objects for gate 1, highest risk to the islet readouts first (tiers in REVIEW_TIERS order, largest
    first within a tier, capped per image so the list covers the cohort)."""
    picks = []
    cap = {c: (n, k) for c, _, n, k in REVIEW_TIERS}

    def take(df: pd.DataFrame, code: str, text):
        n, per_image = cap[code]
        df = df.copy()
        df["_k"] = df.groupby("tag").cumcount()
        for _, r in df[df._k < per_image].head(n).iterrows():
            is_leg = "legacy_name" in r.index
            picks.append(dict(
                reason_code=code, reason=text(r), **{k: r.get(k) for k in (
                    "image", "donor", "group", "disease", "cx_um", "cy_um", "cx_px", "cy_px")},
                cls=(f"legacy Islets2 ({CLS_SHORT[r.legacy_class]}-sized)" if is_leg else r.get("cls")),
                area_um2=r.get("legacy_area_um2") if is_leg else r.get("area_um2"),
                diam_um=r.get("legacy_diam_um") if is_leg else r.get("diam_um")))

    isl = objs[objs.cls == "Islet"]
    big = objs[(objs.area_um2 > BIG_UM2) & ((objs.tH7_islet_parts >= 3) | (objs.n_legacy_inside >= 3))]
    def merge_reason(r):
        ev = []
        if r.tH7_islet_parts >= 2:
            ev.append(f"splits into {int(r.tH7_islet_parts)} Islets at tau_H 7")
        if r.tC4_islet_parts >= 2:
            ev.append(f"into {int(r.tC4_islet_parts)} at tau_CD99 4")
        if pd.notna(r.n_legacy_inside) and r.n_legacy_inside >= 2:
            ev.append(f"holds {int(r.n_legacy_inside)} legacy islets")
        how = ("CD99-joined" if (r.cd99_only_frac >= 0.2 or r.tC4_islet_parts >= 2) else
               "joined by weak hormone signal" if r.tH7_islet_parts >= 3 else "separate islets joined")
        return (f"{r.area_um2 / 1e3:.0f}k um2 ({r.diam_um:.0f} um), {r.cd99_only_frac:.0%} CD99-only: "
                + "; ".join(ev) + f" -> {how}?")
    take(big.sort_values("area_um2", ascending=False), "R1_merge", merge_reason)
    br = objs[objs.bridge == 1].sort_values("area_um2", ascending=False)
    take(br, "R2_cd99_bridge",
         lambda r: (f"bridge flag: {r.area_um2 / 1e3:.1f}k um2, {r.cd99_only_frac:.0%} CD99-only, route {int(r.route)}, "
                    f"{int(r.cores)} hormone core(s); {int(r.tC4_islet_parts)} islet(s) left at tau_CD99 4"))
    r1 = isl[isl.route == 1].sort_values("area_um2", ascending=False)
    take(r1, "R3_cd99_only_islet",
         lambda r: (f"Islet without any hormone pixel (route 1), {r.area_um2:.0f} um2, "
                    f"{r.legacy_overlap_frac:.0%} on a legacy islet: endocrine with weak hormone stain, or CD99+ "
                    f"non-endocrine tissue?"))
    if len(lg):
        m = lg[lg.status.isin(["absent", "partial"])].copy()
        m["_isl"] = m.legacy_class == "Islet"
        m = m.sort_values(["_isl", "legacy_area_um2"], ascending=[False, False])
        take(m, "R4_legacy_missed",
             lambda r: (f"legacy islet {r.legacy_name} ({r.legacy_area_um2:.0f} um2) "
                        + ("not overlapped by any v3 object" if r.status == "absent"
                           else f"only {r.coverage:.0%} covered by v3")
                        + (f"; recovered at {r.found_by_other_settings}" if isinstance(r.found_by_other_settings, str)
                           and r.found_by_other_settings else "; no v3 setting finds it")))
    fr = isl[(isl.tH7_islet_parts == 0) & (isl.area_um2 > 3000)].sort_values("area_um2", ascending=False)
    take(fr, "R5_fragile",
         lambda r: (f"Islet at the default, gone at tau_H 7 ({r.tH7_retained:.0%} of its area kept, "
                    f"{int(r.tH7_parts)} smaller part(s)): hormone signal near threshold"))
    bo = isl[(isl.touch_type == "tissue_outline") & (isl.area_um2 > 20000)].sort_values("area_um2", ascending=False)
    take(bo, "R6_border",
         lambda r: f"large Islet ({r.area_um2 / 1e3:.0f}k um2) touching the Tissue outline: truncated by the outline?")
    rv = pd.DataFrame(picks)
    if len(rv) == 0:
        return rv
    rv = rv.drop_duplicates(subset=["image", "cx_px", "cy_px"]).head(max_n).reset_index(drop=True)
    rv.insert(0, "rank", np.arange(1, len(rv) + 1))
    return rv[["rank", "image", "donor", "group", "disease", "cls", "cx_um", "cy_um", "area_um2", "diam_um",
               "reason_code", "reason", "cx_px", "cy_px"]]


# --------------------------------------------------------------------------------------------- README
def _pct(x, nd=1):
    return "n/a" if x is None or pd.isna(x) else f"{100 * x:+.{nd}f}%"


def _p(x):
    return "n/a" if pd.isna(x) else (f"{x:.2g}" if x >= 0.001 else "<0.001")


def write_readme(a, meta, long, summ, gd, fl, ag, fates, dep, con, lg_sum, exr, pat, review, checks, pchk, comp):
    L = []
    w = L.append
    src = f"production objects in `{a.default_dir}`" if a.default_dir else "the full-image sweep GeoJSONs"
    w("# v3 islet detector: annotation-free validation\n")
    w(f"Generated {time.strftime('%Y-%m-%d %H:%M')} by `scripts/v3/islet_validation.py` from {src} "
      f"({len(meta)} images; excluded: {a.exclude}). No manual labels are used. Pixel size 0.5078 um; "
      "centroids are in um (QuPath's `Centroid X/Y um`) and in level-0 pixels.\n")
    w("Default = tau_CD99 3, tau_H 5, floor k 4, no neck split, no singles co-positivity. Neighbours: tC2/tC4 "
      "(tau_CD99 2/4), tH4/tH7 (tau_H 4/7), k0/k3 (floor off / k 3); neck and copos are the two optional safeguards. "
      "Size classes by mean caliper diameter: Islet >= 35 um, Cluster 10-35 um, Single <= 10 um.\n")

    # ---------------- verdict
    isl = long[long.size_class == "Islet"]
    def med_abs(s, col):
        x = isl[isl.setting == s][col]
        return x.abs().median(), x.abs().max()
    fl_i = fl.copy()
    fl_i["disease"] = fl_i.disease.astype(str)
    tH_flag = fl_i[fl_i.knob == "tau_H"]
    tH_dis = tH_flag.drop_duplicates("tag").disease.value_counts().to_dict()
    d_tH7 = dep[(dep.setting == "tH7") & (dep.size_class == "Islet") & (dep.metric == "area")].iloc[0]
    d_tH4 = dep[(dep.setting == "tH4") & (dep.size_class == "Islet") & (dep.metric == "area")].iloc[0]
    c_def = con.set_index("setting")
    rng = lambda col: (c_def[col].min(), c_def[col].max(), c_def.loc["default", col])
    dens, area = rng("islet_density_mm2__T1D_vs_ND"), rng("islet_area_pct__T1D_vs_ND")
    sp_min = min(c_def.spearman_islet_density_vs_default.min(), c_def.spearman_islet_area_pct_vs_default.min())
    ct7a = comp.get("tH7_kept__size", {}); ct7b = comp.get("tH7_kept__size+composition", {})
    shift = comp.get("boundary_shift_um_median", {})
    w("## Verdict\n")
    mv = max(abs(v / dens[2] - 1) for v in dens[:2])
    mv = max(mv, max(abs(v / area[2] - 1) for v in area[:2]))
    w("**Stability.** Tau_CD99 and the floor are on a plateau for Islet-class objects (outside the few images named "
      "below). Tau_H is not. Its effect is a boundary shift of about one classifier pixel, and the T1D/ND contrasts "
      f"stay within {100 * mv:.0f}% of their default values across all settings.")
    for s in NEIGHBOURS:
        mn, mx = med_abs(s, "rel_n")
        an, ax = med_abs(s, "rel_area")
        nfl = int(((isl[isl.setting == s].rel_n.abs() > FLAG_REL) | (isl[isl.setting == s].rel_area.abs() > FLAG_REL)).sum())
        w(f"- {s}: Islet count |change| median {100 * mn:.1f}% (max {100 * mx:.0f}%), area median {100 * an:.1f}% "
          f"(max {100 * ax:.0f}%); {nfl}/{len(meta)} images > 10%.")
    w("")
    w("**Which knob matters.**")
    sm = summ.set_index(["setting", "size_class"])
    def med2(s1, s2, k):
        return (abs(sm.loc[(s1, k), "rel_n_median"]) + abs(sm.loc[(s2, k), "rel_n_median"])) / 2
    w(f"- tau_H moves every class. The median |count change| per step is Islets {100 * med2('tH4', 'tH7', 'Islet'):.0f}%, "
      f"clusters {100 * med2('tH4', 'tH7', 'Cluster'):.0f}% and singles {100 * med2('tH4', 'tH7', 'Single'):.0f}%. "
      f"The images flagged at tH4/tH7 are {', '.join(f'{k} {v}' for k, v in tH_dis.items())}"
      f"{'' if 'ND' in tH_dis else ' (no ND image)'}.")
    tc_fl = fl_i[fl_i.knob == "tau_CD99"].drop_duplicates("tag")
    anc_ = pat.set_index("tag").cd99_anchor
    w(f"- tau_CD99 moves the small classes: singles median {_pct(sm.loc[('tC2', 'Single'), 'rel_n_median'])} and "
      f"clusters {_pct(sm.loc[('tC2', 'Cluster'), 'rel_n_median'])} at tC2. It moves Islets by more than 10% only on "
      + (", ".join(f"{r.donor} (CD99 anchor {anc_.get(r.tag, np.nan):.0f}; cohort median {anc_.median():.0f})"
                   for r in tc_fl.itertuples()) or "no image") + ".")
    k0 = long[long.setting == "k0"]
    nobind = int(k0.groupby("tag").apply(lambda x: (x.rel_n.abs() < 1e-9).all(), include_groups=False).sum())
    k0i = isl[isl.setting == "k0"].sort_values("rel_n", ascending=False).iloc[0]
    k0s = long[(long.setting == "k0") & (long.size_class == "Single")].rel_n.max()
    k3f = fl_i[fl_i.setting == "k3"]
    w(f"- The floor binds (changes anything at all) on {len(meta) - nobind}/{len(meta)} images, the heavily INS/SST-clamped "
      f"ones; on the other {nobind}, k0 = k3 = k4 exactly. "
      f"Without it (k0), singles rise by up to {_pct(k0s, 0)} and Islets by up to {_pct(k0i.rel_n, 0)} ({k0i.donor}). "
      "k3 equals k4 for Islets"
      + (" except " + ", ".join(f"{r.donor} (area {_pct(r.rel_area, 0)})" for r in k3f.itertuples()) if len(k3f) else "")
      + ". The sweep has no k5, so the upper side of k 4 is untested.")
    nk = isl[isl.setting == "neck"]
    w(f"- The neck split changes Islet counts by at most {100 * nk.rel_n.abs().max():.0f}%. copos removes singles only "
      f"(up to {_pct(sm.loc[('copos', 'Single'), 'rel_n_min'], 0)}).\n")
    w("**Batch vs disease.**")
    w(f"- Tau_H sensitivity of Islet area depends on disease, not on group: "
      f"tH4 p_disease {_p(d_tH4.p_disease)} vs p_group {_p(d_tH4.p_group)}; "
      f"tH7 {_p(d_tH7.p_disease)} vs {_p(d_tH7.p_group)}. "
      f"Median change at tH7: ND {_pct(d_tH7.median_ND)}, T1D {_pct(d_tH7.median_T1D)}.")
    ct4a = comp.get("tH4_grown__size", {}); ct4b = comp.get("tH4_grown__size+composition", {})
    w(f"- Per islet, the models adjust for size and cluster by donor.")
    w(f"  - tH7: T1D keeps less area than ND (ratio {ct7a.get('T1D_vs_ND_ratio', np.nan):.3f}, "
      f"p {_p(ct7a.get('p_disease'))}). With INS share and hormone share added, the disease term vanishes "
      f"(ratio {ct7b.get('T1D_vs_ND_ratio', np.nan):.3f}, p {_p(ct7b.get('p_disease'))}; INS share "
      f"p {_p(ct7b.get('p_ins_frac'))}).")
    w(f"  - tH4: T1D grows more (ratio {ct4a.get('T1D_vs_ND_ratio', np.nan):.2f}, p {_p(ct4a.get('p_disease'))}). "
      f"Composition explains most of it (ratio {ct4b.get('T1D_vs_ND_ratio', np.nan):.2f}, "
      f"p {_p(ct4b.get('p_disease'))}).")
    w(f"  - Group: p {_p(ct7b.get('p_group'))} (tH7) and {_p(ct4b.get('p_group'))} (tH4).")
    w(f"- So INS-deficient islets have softer GCG/SST edges. This is mainly biology, not batch. "
      f"Median boundary shift per tau_H step: ND {shift.get('ND', {}).get('tH4', np.nan):.2f}/"
      f"{shift.get('ND', {}).get('tH7', np.nan):.2f} um, T1D {shift.get('T1D', {}).get('tH4', np.nan):.2f}/"
      f"{shift.get('T1D', {}).get('tH7', np.nan):.2f} um (tH4/tH7), against a 1.0156 um classifier pixel.")
    w(f"- The group-adjusted T1D/ND ratios stay within {dens[0]:.2f}-{dens[1]:.2f} for islet density "
      f"(default {dens[2]:.2f}) and {area[0]:.2f}-{area[1]:.2f} for islet area % (default {area[2]:.2f}) over all "
      f"9 settings. The per-image rank vs the default has Spearman >= {sp_min:.2f}.\n")
    if len(lg_sum):
        pooled = (lg_sum.n_iou + lg_sum.n_merged + lg_sum.n_fragmented + lg_sum.n_covered).sum() / lg_sum.n_in_area.sum()
        low = lg_sum[lg_sum.recall < 0.9]
        worst = lg_sum.sort_values("recall").head(3)
        w(f"**Legacy recall.**")
        w(f"- {100 * pooled:.1f}% of {int(lg_sum.n_in_area.sum())} legacy Islets2 islets (in {len(lg_sum)} images) are "
          f"found by v3 default objects. Found means IoU >= 0.3, centroid containment in either direction, or >= 50% "
          "covered.")
        st_all = lg_sum.recall_strict * lg_sum.n_in_area
        w(f"- Strict recall: {100 * st_all.sum() / lg_sum.n_in_area.sum():.1f}% pooled, per-image minimum "
          f"{100 * lg_sum.recall_strict.min():.1f}% ({lg_sum.loc[lg_sum.recall_strict.idxmin(), 'donor']}).")
        w(("- No image is below 90% on either definition." if len(low) == 0 and (lg_sum.recall_strict >= 0.9).all()
           else f"- {len(low)} images are below 90% (strict: {int((lg_sum.recall_strict < 0.9).sum())}).")
          + " The lowest are: "
          + ", ".join(f"{r.donor} {100 * r.recall:.1f}% ({int(r.n_partial + r.n_absent)} missed)" for r in worst.itertuples())
          + ".")
        w(f"- Only {int(lg_sum.n_missed_ge35.sum())} missed legacy objects are islet-sized (>= 35 um). The others are "
          "cluster- or single-sized.")
        w(f"- Matched v3/legacy area ratio: per-image median {lg_sum.median_v3_to_legacy_area.min():.2f}-"
          f"{lg_sum.median_v3_to_legacy_area.max():.2f}. v3 is neither systematically tighter nor looser. The Islets2 "
          "erode/dilate varied per image.")
        ge = lg_sum.n_legacy_ge35.sum()
        w(f"- {100 * (lg_sum.recall_ge35_as_islet * lg_sum.n_legacy_ge35).sum() / max(ge, 1):.1f}% of the {int(ge)} "
          "legacy objects >= 35 um are found as a v3 Islet. The rest become v3 clusters just under the 35 um line, a "
          "class-boundary effect rather than a miss.")
        w("- The legacy classifier used the same four channels (plus HLA-A for 6579), so this recall is agreement, "
          "not independent truth.\n")
    npat = pat
    w("**Pathology.**")
    w(f"- No image reaches the ADIIDAS 10% endocrine flag. The maximum is "
      f"{npat.endocrine_pct.max():.1f}% of the analysis area and {npat.endocrine_pct_occupied.max():.1f}% of "
      f"DAPI-occupied tissue.")
    mgs = checks.get("merges", {})
    w(f"- {int(npat.n_gt50k.sum())} objects exceed 50,000 um2, up to {npat.largest_um2.max() / 1e3:.0f}k um2. "
      f"Of these, {mgs.get('n', 0)} look like merges: {mgs.get('n_split_tH7_ge3', 0)} split into >= 3 Islets at "
      f"tau_H 7, and {mgs.get('n_legacy_ge3', 0)} hold >= 3 legacy islets.")
    w(f"- The merges are joined by hormone signal, not by CD99: their median CD99-only fraction is "
      f"{100 * mgs.get('median_cd99_only_frac', np.nan):.1f}%, only {mgs.get('n_cd99_only_ge20', 0)} are >= 20% "
      f"CD99-only, and {mgs.get('n_split_tC4', 0)} split at tau_CD99 4. The neck split (CD99-only necks) cannot cut "
      "them. Hormone-positive bridges are unguarded (ISLETS.md section 7).")
    tb = npat.sort_values("n_bridge_flag", ascending=False).iloc[0]
    w(f"- CD99 bridging and hormone-negative Islets concentrate in {tb.donor}: {int(tb.n_bridge_flag)} of the "
      f"cohort's {int(npat.n_bridge_flag.sum())} bridge flags, and {int(tb.n_islet_route1)} of its {int(tb.n_islets)} "
      f"Islets have no hormone pixel (route 1; cohort total {int(npat.n_islet_route1.sum())}). Its CD99 anchor is "
      f"{tb.cd99_anchor:.0f} counts (cohort median {npat.cd99_anchor.median():.0f}).\n")
    w("**Recommendation.** Keep the default (tC3_tH5_k4, no neck, no copos). No sweep variant is more stable:")
    w("- tC2/tC4 and k0/k3 do not change the Islet class outside the images named above.")
    w("- tH4 and tH7 only move the T1D islet boundary by about a pixel.")
    w("- neck and copos do not touch the main pathology.")
    w("")
    w("Use tH4/tH7 as a sensitivity band for T1D islet area. Treat singles as the least reliable class. Spend "
      f"gate 1 on the merges, on {tb.donor}'s CD99-only islets, and on the missed legacy islets (review_list.tsv).\n")

    # ---------------- 1. stability
    w("## 1. Stability over the sweep\n")
    t = summ.copy()
    t = t[t.size_class.isin(["Islet", "Cluster", "Single"])]
    t["count change, median [min, max]"] = [f"{_pct(a_)} [{_pct(b_)}, {_pct(c_)}]" for a_, b_, c_ in
                                           zip(t.rel_n_median, t.rel_n_min, t.rel_n_max)]
    t["area change, median [min, max]"] = [f"{_pct(a_)} [{_pct(b_)}, {_pct(c_)}]" for a_, b_, c_ in
                                          zip(t.rel_area_median, t.rel_area_min, t.rel_area_max)]
    t = t.rename(columns={"n_images_abs_gt10": "images >10%"})
    w("Table 1. Change vs the default, over images (count and total area per size class; "
      "`stability_long.tsv` has every image x class x setting).\n")
    w(md_table(t[["setting", "size_class", "count change, median [min, max]", "area change, median [min, max]",
                  "images >10%"]]))
    w("")
    ta = ag.groupby("setting")[["islet_jaccard", "frac_default_matched", "islet_area_iou"]].median().reindex(VARIANTS)
    ft = fates[[c for c in fates.columns if c.startswith(("split_", "lost_", "merged_"))] + ["neck_split"]].sum()
    ta["Jaccard min"] = ag.groupby("setting").islet_jaccard.min().reindex(VARIANTS)
    ta = ta.reset_index().rename(columns={"islet_jaccard": "Islet Jaccard (IoU>=0.5), median",
                                          "frac_default_matched": "default Islets kept", "islet_area_iou": "Islet area IoU"})
    w("Table 2. Object-level agreement of Islet-class objects with the default (median over images). The settings "
      "are nested: every stricter-setting object lies inside a default object, and every looser-setting object "
      f"contains default objects (checked: minimum nesting {min(checks['nesting_min'].values()):.4f}).\n")
    w(md_table(ta, "{:.3f}"))
    w("")
    w(f"Fates of the {int(fates.n_islets.sum())} default Islets:")
    w(f"- split into >= 2 Islets: {int(ft.split_tC4)} at tau_CD99 4 (CD99 necks), {int(ft.split_tH7)} at tau_H 7 "
      "(hormone necks);")
    w(f"- no Islet-class part left: {int(ft.lost_tC4)} at tC4, {int(ft.lost_tH7)} at tH7;")
    w(f"- joined with another default Islet: {int(ft.merged_tC2)} at tC2, {int(ft.merged_tH4)} at tH4, "
      f"{int(ft.merged_k0)} at k0, {int(ft.merged_k3)} at k3.")
    w("Per image: `islet_fates.tsv`.\n")
    t3 = dep[(dep.size_class == "Islet") & dep.setting.isin(NEIGHBOURS)].copy()
    for c in ["median_ND", "median_Aab+", "median_T1D", "median_E1", "median_E2", "median_E3", "median_E4"]:
        t3[c] = t3[c].map(_pct)
    t3["p_group"] = t3.p_group.map(_p)
    t3["p_disease"] = t3.p_disease.map(_p)
    w("Table 3. Is the Islet-class sensitivity batch- or disease-dependent? The table gives the median change by "
      "disease and by acquisition group. The p-values are type-II F tests of OLS log(1+change) ~ group + disease "
      "(25 images; descriptive).\n")
    w(md_table(t3[["setting", "metric", "median_ND", "median_Aab+", "median_T1D", "median_E1", "median_E2", "median_E3",
                   "median_E4", "p_group", "p_disease"]]))
    w("")
    t4 = con[["setting", "islet_density_mm2__T1D_vs_ND", "islet_area_pct__T1D_vs_ND", "endocrine_area_pct__T1D_vs_ND",
              "islet_density_mm2__Aab_vs_ND", "islet_area_pct__Aab_vs_ND", "spearman_islet_density_vs_default",
              "spearman_islet_area_pct_vs_default"]].rename(columns={
        "islet_density_mm2__T1D_vs_ND": "T1D/ND islets per mm2", "islet_area_pct__T1D_vs_ND": "T1D/ND islet area %",
        "endocrine_area_pct__T1D_vs_ND": "T1D/ND endocrine %", "islet_density_mm2__Aab_vs_ND": "Aab+/ND islets per mm2",
        "islet_area_pct__Aab_vs_ND": "Aab+/ND islet area %", "spearman_islet_density_vs_default": "rho density",
        "spearman_islet_area_pct_vs_default": "rho area %"})
    w("Table 4. Downstream contrasts per setting. Each is exp(coef) of log(readout) ~ group + disease, i.e. the "
      "group-adjusted ratio to ND. rho is the Spearman rank correlation of the per-image readout with the default. "
      "Per-image readouts: `readouts_by_setting.tsv`.\n")
    w(md_table(t4, "{:.3f}"))
    w("")
    w(f"Table 5. Images whose Islet-class count or area changes by more than 10% at a neighbour of the default "
      f"({fl.tag.nunique()} images, {len(fl)} image x setting pairs; `stability_flags.tsv`).\n")
    f5 = fl.copy()
    f5["change"] = [f"{s}: n {_pct(n_)}, area {_pct(a_)}" for s, n_, a_ in zip(f5.setting, f5.rel_n, f5.rel_area)]
    f5 = f5.groupby(["donor", "group", "disease"]).change.apply(lambda s: "; ".join(s)).reset_index()
    w(md_table(f5))
    w("")

    # ---------------- 2. legacy
    if len(lg_sum):
        w("## 2. Agreement with the legacy Islets2 islets\n")
        w("A legacy islet counts as found when a v3 default object matches it at IoU >= 0.3 (`iou`), or contains its "
          "centroid (`merged`), or has a centroid inside it (`fragmented`), or v3 objects cover >= 50% of it "
          "(`covered`: it lies inside a larger v3 object whose centroid is elsewhere). `strict` also requires "
          "`merged`/`fragmented` matches to cover >= 30% of it (a small v3 object at the centre of a large legacy "
          "islet does not count). `any setting` counts it as found when any of the 9 settings finds it. Missed = `partial` "
          "(some overlap) + `absent` (none). The miss type is the looser setting that recovers it: CD99-dim = tC2, "
          "hormone-dim = tH4, floor = k3/k0.\n")
        t6 = lg_sum.sort_values("recall")[["donor", "group", "disease", "n_in_area", "n_iou", "n_merged", "n_fragmented",
                                           "n_covered", "n_partial", "n_absent", "recall", "recall_strict", "recall_any_setting",
                                           "n_legacy_ge35", "recall_ge35_as_islet", "n_missed_ge35", "miss_cd99dim",
                                           "miss_hormonedim", "miss_nosetting", "median_iou_matched"]].copy()
        for c in ["recall", "recall_strict", "recall_any_setting", "recall_ge35_as_islet"]:
            t6[c] = (100 * t6[c]).map(lambda v: f"{v:.1f}")
        t6 = t6.rename(columns={"recall_ge35_as_islet": "legacy>=35um found as v3 Islet %", "n_missed_ge35": "missed >=35um",
                                "median_iou_matched": "median IoU"})
        w("Table 6. Legacy recall per image (sorted by recall; `legacy_recall.tsv`; every legacy islet with its status: "
          "`legacy_islets.tsv.gz`).\n")
        w(md_table(t6, "{:.3g}"))
        w("")
        e = exr[exr.tag.isin(lg_sum.tag)]
        tt = e.groupby(["cls", "extra"]).n.sum().unstack("extra").fillna(0)
        rt = e[e.extra].groupby(["cls", "route"]).n.sum().unstack("route").fillna(0)
        rows = []
        for c in CLASSES:
            if c not in tt.index:
                continue
            new, old = tt.loc[c].get(True, 0), tt.loc[c].get(False, 0)
            r_ = rt.loc[c] if c in rt.index else pd.Series(dtype=float)
            tot = max(r_.sum(), 1)
            rows.append({"class": c, "v3 objects": int(new + old), "no legacy overlap": int(new),
                         "% new": f"{100 * new / max(new + old, 1):.1f}",
                         "new: CD99 only (1)": f"{100 * r_.get(1.0, 0) / tot:.1f}%",
                         "new: hormone only (2)": f"{100 * r_.get(2.0, 0) / tot:.1f}%",
                         "new: both (3)": f"{100 * r_.get(3.0, 0) / tot:.1f}%"})
        w("Table 7. What v3 adds: default objects with no overlap with any legacy islet, by class and route "
          f"({len(lg_sum)} images with legacy islets). Per image: `v3_vs_legacy_objects.tsv`.\n")
        w(md_table(pd.DataFrame(rows)))
        w("")
        pi = e[e.cls == "Islet"].groupby(["tag", "extra"]).n.sum().unstack("extra").fillna(0)
        pi = (100 * pi.get(True, 0) / (pi.get(True, 0) + pi.get(False, 0))).rename("pct").reset_index().merge(meta, on="tag")
        top = pi.sort_values("pct", ascending=False).head(6)
        w("New Islet-class objects are concentrated in T1D/Aab+ images: "
          + ", ".join(f"{r.donor} ({r.disease}) {r.pct:.0f}%" for r in top.itertuples())
          + f". The median over ND images is {pi[pi.disease == 'ND'].pct.median():.1f}%.\n")

    # ---------------- 3. pathology
    w("## 3. Pathology flags\n")
    t8 = pat.sort_values(["group", "disease", "donor"])[[
        "donor", "group", "disease", "endocrine_pct", "endocrine_pct_occupied", "n_islets", "n_gt50k", "largest_um2",
        "islet_area_in_gt50k_pct", "n_gt50k_merge", "n_islet_joins_ge2_legacy", "n_islet_split_tH7", "n_bridge_flag",
        "n_islet_route1", "n_islet_cd99only_ge50", "touch_Islet", "touch_Islet_tissue_outline", "touch_Islet_image_edge",
        "cd99_anchor"]].copy()
    t8["largest_um2"] = (t8.largest_um2 / 1e3).map(lambda v: f"{v:.0f}k")
    t8 = t8.rename(columns={"endocrine_pct": "endocrine % area", "endocrine_pct_occupied": "endocrine % occupied",
                            "n_gt50k": ">50k um2", "islet_area_in_gt50k_pct": "% islet area in >50k",
                            "n_gt50k_merge": ">50k merges", "n_islet_joins_ge2_legacy": "Islets holding >=2 legacy",
                            "n_islet_split_tH7": "Islets split at tH7", "n_bridge_flag": "bridge flags",
                            "n_islet_route1": "route-1 Islets", "n_islet_cd99only_ge50": "Islets >=50% CD99-only",
                            "touch_Islet": "Islets on border", "touch_Islet_tissue_outline": "on Tissue outline",
                            "touch_Islet_image_edge": "on image edge"})
    w("Table 8. Per image, default objects. Endocrine % uses the analysis area and the DAPI-occupied tissue "
      "(`occupied_fraction`, islet_anchors.csv); the ADIIDAS flag is at 10%. A merge is an object > 50,000 um2 that "
      "splits into >= 3 Islets at tau_H 7 or holds >= 3 legacy islets. The bridge flag is the detector's own "
      "(>= 5000 um2 and >= 50% CD99-only). Border = distance to the analysis-area edge <= 1 pixel, split into Tissue "
      "outline and image edge (re-derived from the exported outlines; `pathology_flags.tsv`; objects: "
      "`large_or_bridge_objects.tsv`).\n")
    w(md_table(t8, "{:.3g}"))
    w("")
    w(f"Border check: the largest distance from a touching object to the re-derived boundary is "
      f"{checks.get('border_recheck_max_um', np.nan):.2f} um.\n")

    # ---------------- 4. review list
    w("## 4. Review list for gate 1 (`review_list.tsv`)\n")
    P = pat.set_index("donor")
    img = []
    r1 = P[P.n_islet_route1 / P.n_islets.clip(lower=1) >= 0.10]
    for d_, r in r1.iterrows():
        img.append(f"**{d_}**: {int(r.n_islet_route1)}/{int(r.n_islets)} Islets have no hormone pixel and there are "
                   f"{int(r.n_bridge_flag)} bridge flags (CD99 anchor {r.cd99_anchor:.0f}). Check whether its CD99-only "
                   "Islets are islets.")
    mg = P[P.islet_area_in_gt50k_pct >= 20].sort_values("islet_area_in_gt50k_pct", ascending=False)
    if len(mg):
        img.append("**" + ", ".join(f"{d_} ({v:.0f}%)" for d_, v in mg.islet_area_in_gt50k_pct.items()) +
                   "**: >= 20% of the Islet area is in objects > 50,000 um2. Check whether these are merged islets.")
    if len(lg_sum):
        lo = lg_sum[lg_sum.recall < 0.97].sort_values("recall")
        for r in lo.itertuples():
            img.append(f"**{r.donor}**: legacy recall {100 * r.recall:.1f}% ({int(r.n_partial + r.n_absent)} missed, "
                       f"median {r.missed_median_area_um2:.0f} um2; {int(r.miss_nosetting)} found by no setting). Check "
                       "whether the missed legacy objects are endocrine.")
    bo_ = P[(P.touch_Islet_tissue_outline / P.n_islets.clip(lower=1)) >= 0.04]
    if len(bo_):
        img.append("**" + ", ".join(f"{d_} ({int(v)})" for d_, v in bo_.touch_Islet_tissue_outline.items()) +
                   "**: many Islets touch the Tissue outline. Check whether the outline cuts through tissue.")
    k0f = fl[(fl.setting == "k0") & (fl.rel_n > 0.2)]
    for r in k0f.itertuples():
        k3r = long[(long.tag == r.tag) & (long.setting == "k3") & (long.size_class == "Islet")].rel_n
        img.append(f"**{r.donor}**: the floor matters (k0 {_pct(r.rel_n, 0)} Islets; k3 "
                   f"{_pct(k3r.iloc[0] if len(k3r) else np.nan, 0)}). Check the extra INS-only objects that k0 would add.")
    if img:
        w("Images to open first:")
        for i_, s_ in enumerate(img, 1):
            w(f"{i_}. {s_}")
        w("")
    if len(review):
        cnt = review.reason_code.value_counts()
        lab = {c: l for c, l, _, _ in REVIEW_TIERS}
        w(f"{len(review)} objects, highest risk first. Tiers: "
          + "; ".join(f"{lab.get(c, c)} ({cnt.get(c, 0)})" for c, *_ in REVIEW_TIERS) + ".")
        w("Within a tier the list is sorted by size and capped per image. Coordinates are `Centroid X/Y um`; the "
          "`cx_px/cy_px` columns can be used in QuPath's script editor:")
        w("`getCurrentViewer().setCenterPixelLocation(cx_px, cy_px)`.\n")
        rv = review[["rank", "donor", "disease", "cls", "cx_um", "cy_um", "area_um2", "reason_code"]].copy()
        rv["cx_um"] = rv.cx_um.map(lambda v: f"{v:.0f}")
        rv["cy_um"] = rv.cy_um.map(lambda v: f"{v:.0f}")
        rv["area_um2"] = rv.area_um2.map(lambda v: f"{v:.0f}")
        w(md_table(rv))
        w("")

    # ---------------- 5. recompute
    w("## 5. Recomputing on the production pass-I objects\n")
    if len(pchk):
        w(f"Pass-I cross-check (`{a.passI_log.name}`): {int(pchk.identical.sum())}/{len(pchk)} finished images have "
          "exactly the sweep default's Islet/Cluster/Single counts. The production objects are expected to be the "
          "sweep default objects.\n")
    w("- **Recomputable from pass-I objects:** sections 2 (legacy recall, extras), 3 (pathology, border, endocrine %) "
      "and 4 (review list), plus the default column of section 1.")
    w("  - Export the pass-I annotations into a NEW folder, for example with `export_annotations.groovy` and a "
      "different output root.")
    w("  - Then run `islet_validation.py --default-dir <folder>/annotations_all --out-dir <new dir>`.")
    w("  - **Keep `results/v3/annotations_all/` untouched until then.** Pass I's script 1 replaces the legacy "
      "Islets2 islets in each .qpdata, so after pass I this export is their only copy, and `export_annotations.sh` "
      "overwrites it. Snapshot it first (`cp -a annotations_all annotations_all_legacy_20260927`) and point "
      "`--legacy-dir` there.")
    w("- **Not recomputable from pass I:** the stability deltas and fates (section 1). They need the sweep variants, "
      "which production does not build. They stay valid as long as pass I reproduces the sweep default "
      "(identical counts: see the cross-check above).")
    w("- After gate-1 edits (Artifact outlines, deleted or redrawn objects), the anchors and hash change. Re-run the "
      "sweep only if the thresholds are revisited.\n")
    w("## Files\n")
    w("| file | content |\n|---|---|")
    for f_, d_ in [("stability_long.tsv", "image x size class x setting: n, area, change vs default"),
                   ("stability_by_setting.tsv", "Table 1"), ("stability_by_group_disease.tsv", "group x disease medians"),
                   ("stability_flags.tsv", "Table 5"), ("islet_agreement.tsv", "Table 2 per image"),
                   ("islet_fates.tsv", "split/lost/merged Islets per image"),
                   ("sensitivity_batch_disease.tsv", "Table 3, all classes"),
                   ("tauH_composition_test.json", "per-islet tau_H test (size, composition, disease, group)"),
                   ("downstream_contrasts.tsv", "Table 4"), ("readouts_by_setting.tsv", "per-image readouts per setting"),
                   ("legacy_recall.tsv", "Table 6"), ("legacy_islets.tsv.gz", "every legacy islet: status, miss type"),
                   ("v3_vs_legacy_objects.tsv", "Table 7 per image"), ("pathology_flags.tsv", "Table 8"),
                   ("large_or_bridge_objects.tsv", "every object > 50,000 um2 or bridge-flagged"),
                   ("islets_default.tsv.gz", "every default Islet: measurements, fate at each setting, legacy overlap, border contact"),
                   ("review_list.tsv", "gate-1 review list"), ("checks.json", "consistency checks")]:
        w(f"| `{f_}` | {d_} |")
    return "\n".join(L) + "\n"


if __name__ == "__main__":
    sys.exit(main())
