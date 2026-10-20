#!/usr/bin/env python3
"""Same-slide RNA check of REDSEA spillover correction (IO60panc2nd v3; donors 6533 and HDL112).

Question, per protein marker that has a gene on the Xenium panel: on the SAME physical section,
does REDSEA-corrected protein agree with the same cell's RNA better or worse than raw protein?
Secondary evidence for the per-marker REDSEA policy (primary = impossible co-positives).

Inputs (read-only)
  v3 PhenoCycler   results/v3/cells/<image>.parquet and results/v3/cell_geojson/<image>.geojson
                   (paths from results/v3/senior/qupath_manifest_v3.json)
  Senior tables    results/v3/senior/{data,pilot_data}/{cells,redsea,geometry_qc}/donor_id=<d>/
                   raw = QuPath whole-cell mean (cells table), corrected = redsea ``Cell__<m>``
  Xenium bundle    cells.parquet, nucleus_boundaries.parquet, cell_feature_matrix.h5
                   (gene_panel.json / metrics_summary.csv optional). HiPerGator only:
                   /blue/maigan/smith6jt/Xenium_Analysis/data/output-XETG00298__0080870__panc__20260910__211917 (6533)
                   /blue/maigan/smith6jt/Xenium_Analysis/data/output-XETG00298__0081366__panc__20260910__211918 (HDL112)
  Transform        optional: an existing Xenium->PhenoCycler Transform (Senior ``Transform.save``
                   format, e.g. 6533's ``refined_field2``). Without one it is refitted from the two
                   nucleus clouds.

Subcommands
  nuclei   v3 nucleus centroids (um) from the cell GeoJSON, plus the 25 um placement null
  protein  per-marker REDSEA zeroing (raw>0 -> 0) and which markers have a gene on the panel
  run      register Xenium nuclei onto v3 nuclei -> validate with the lane's nucleus-pair
           statistic -> mutual-nearest nucleus match -> per-marker protein-vs-RNA agreement

Matching rule: Xenium nucleus centroids are warped FORWARD into the v3 frame (the lane's
convention) and paired with v3 nucleus centroids by mutual nearest neighbour within --cap-um
(3.5 um = half a nucleus diameter, as in the lane's fan-out analysis).

RNA positivity is gDNA-aware (docs/XENIUM_PIPELINE.md sec. 1): lambda_c is the shrunken per-gene
genomic-control rate; a gene set of k genes is RNA+ when counts >= 2 AND P(K >= counts | k*lambda_c)
< 1e-3, RNA- when counts == 0, and ambiguous (excluded) otherwise. Cells need >= 50 gDNA-corrected
transcripts.

Resource use: run under ``taskset -c <8 cores>``; the Senior helpers call cKDTree(workers=-1).
Outputs: results/v3/redsea_policy/rna_check/<donor>/ (never under results/v3/senior).
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

V3 = Path(os.environ.get("IO60_V3", "/home/smith6jt/IO60panc2nd/results/v3"))
SENIOR_OUT = V3 / "senior"
MANIFEST = SENIOR_OUT / "qupath_manifest_v3.json"
OUT = V3 / "redsea_policy" / "rna_check"
SENIOR_REPO = Path(os.environ.get("SENIOR_REPO", "/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior"))
#: A local Prime 5K gene_panel.json (Y9QRJ6 add-on design). Its 'base' category is the stock
#: 5,001-gene panel the same-slide runs used; 'current' genes are the 100-gene T1D add-on.
REF_GENE_PANEL = Path("/data/QP_Xenium/output-XETG00298__0041435__TMA172_PANC__20250606__213055/"
                      "gene_panel.json")

DONORS = {
    "6533": {"serial": "0080870", "bundle": "output-XETG00298__0080870__panc__20260910__211917"},
    "0112": {"serial": "0081366", "bundle": "output-XETG00298__0081366__panc__20260910__211918"},
}
ALIASES = {"6533": "6533", "0112": "0112", "112": "0112", "HDL112": "0112"}

#: Channels REDSEA passes through uncorrected (full-run log, "[no-spillover]"). Kept as
#: negative controls: their delta-AUC measures raster-vs-QuPath measurement noise only.
PASS_THROUGH = {"DAPI", "FOXP3", "PCNA", "Ki67", "Collagen_IV", "TCF_1", "TOX", "SOX2", "TP63"}

#: marker -> list of (kind, genes). The first kind whose genes are on the bundle is used.
#: direct = the protein's own gene; proxy = co-expressed genes (weaker, labelled as such).
MARKER_GENES: dict[str, list[tuple[str, list[str]]]] = {
    "INS": [("direct", ["INS"]), ("proxy", ["PCSK1"])],
    "GCG": [("direct", ["GCG"]), ("proxy", ["ARX"])],
    "SST": [("direct", ["SST"]), ("proxy", ["HHEX"])],
    "CD3e": [("direct", ["CD3E"])],
    "CD20": [("direct", ["MS4A1"])],
    "CD4": [("direct", ["CD4"])],
    "CD8": [("direct", ["CD8A"])],
    "E_cadherin": [("direct", ["CDH1"])],
    "Vimentin": [("direct", ["VIM"])],
    "Caveolin": [("direct", ["CAV1"]), ("proxy", ["CAV2", "CAVIN2"])],
    "SMA": [("direct", ["ACTA2"]), ("proxy", ["MYLK", "SMTN", "NOTCH3", "SYNPO2", "KCNMB1"])],
    # pan-epithelial proxy only: no ductal/acinar claim is made (Senior CLAUDE.md)
    "Pan_Cytokeratin": [("direct", ["KRT8", "KRT18", "KRT19", "KRT7"]),
                        ("proxy", ["EPCAM", "CDH1", "CLDN7"])],
    "Keratin_5": [("direct", ["KRT5"])],
    "CD31": [("direct", ["PECAM1"])],
    "CD68": [("direct", ["CD68"])],
    "CD99": [("direct", ["CD99"])],
    "EpCAM": [("direct", ["EPCAM"])],
    "CD34": [("direct", ["CD34"])],
    "Podoplanin": [("direct", ["PDPN"])],
    "HLA_DR": [("direct", ["HLA-DRA"]), ("proxy", ["CIITA"])],
    "HLA_A": [("direct", ["HLA-A"])],
    "CD38": [("direct", ["CD38"])],
    "CD44": [("direct", ["CD44"])],
    "CD11c": [("direct", ["ITGAX"])],
    "CD11b": [("direct", ["ITGAM"])],
    "CD21": [("direct", ["CR2"])],
    "iNOS": [("direct", ["NOS2"])],
    "B3TUBB": [("direct", ["TUBB3"])],
    "CD57": [("direct", ["B3GAT1"])],
    "Granzyme_B": [("direct", ["GZMB"])],
    "VISTA": [("direct", ["VSIR"])],
    "PD_1": [("direct", ["PDCD1"])],
    "PD_L1": [("direct", ["CD274"])],
    "b_Catenin1": [("direct", ["CTNNB1"])],
    "ICOS": [("direct", ["ICOS"])],
    "LAG3": [("direct", ["LAG3"])],
    "CD163": [("direct", ["CD163"])],
    "CD206": [("direct", ["MRC1"])],
    "CD209": [("direct", ["CD209"])],
    "Iba1": [("direct", ["AIF1"]), ("proxy", ["CD14", "CSF1R"])],
    "M2Gal3": [("direct", ["LGALS3"])],
    "CD79a": [("direct", ["CD79A"])],
    "Bcl_2": [("direct", ["BCL2"])],
    "CD56": [("direct", ["NCAM1"])],
    "CD141": [("direct", ["THBD"])],
    "CD39": [("direct", ["ENTPD1"])],
    "CD66": [("direct", ["CEACAM1", "CEACAM8"])],
    "CD107a": [("direct", ["LAMP1"])],
    "IDO1": [("direct", ["IDO1"])],
    "MPO": [("direct", ["MPO"])],
    "Beta_actin": [("direct", ["ACTB"])],
    "ATP5A": [("direct", ["ATP5F1A"])],
    "GLUT1": [("direct", ["SLC2A1"])],
    "LDHA": [("direct", ["LDHA"])],
    "HK1": [("direct", ["HK1"])],
    "IDH2": [("direct", ["IDH2"])],
    "SDHA": [("direct", ["SDHA"])],
    "ASCT2": [("direct", ["SLC1A5"])],
    "G6PD": [("direct", ["G6PD"])],
    "CS": [("direct", ["CS"])],
    # pass-through controls
    "FOXP3": [("direct", ["FOXP3"])],
    "PCNA": [("direct", ["PCNA"])],
    "Ki67": [("direct", ["MKI67"])],
    "TCF_1": [("direct", ["TCF7"])],
    "TOX": [("direct", ["TOX"])],
    "SOX2": [("direct", ["SOX2"])],
    "TP63": [("direct", ["TP63"])],
    "Collagen_IV": [("direct", ["COL4A1", "COL4A2"])],
}

CAP_UM = 3.5            # mutual-NN cap: half a nucleus diameter (lane fan-out analysis)
CELL_GATE_UM = 2.0      # lane gate (SAME_SECTION_LINKAGE): held-out median nucleus offset
MIN_CORRECTED_COUNTS = 50
RNA_POS_P = 1e-3
N_FOLDS = 20
FOLD_BLOCK_UM = 1000.0


def log(msg: str) -> None:
    print(f"[rna_check {time.strftime('%H:%M:%S')}] {msg}", flush=True)


def canon(donor: str) -> str:
    try:
        return ALIASES[str(donor)]
    except KeyError:
        raise SystemExit(f"unknown donor {donor!r}; expected one of {sorted(ALIASES)}")


def senior_imports():
    """Senior's same-slide lane helpers, imported read-only from the checkout."""
    core = SENIOR_REPO / "core" / "phenocycler"
    for p in (core, core / "scripts" / "restain", core / "scripts" / "sameslide"):
        if str(p) not in sys.path:
            sys.path.insert(0, str(p))
    import refine as rf                                           # noqa: E402
    import register_points as rp                                  # noqa: E402
    from phenocycler.integration.rasterize import GridSpec        # noqa: E402
    from phenocycler.integration.transform import Transform       # noqa: E402
    return rp, rf, Transform, GridSpec


def manifest_entry(donor: str) -> dict:
    m = json.loads(MANIFEST.read_text())
    for r in m["images"]:
        if str(r["donor_id"]) == donor:
            return r
    raise SystemExit(f"donor {donor} not in {MANIFEST}")


def senior_root(donor: str) -> Path:
    """The Senior output root whose REDSEA partition exists for this donor (full run first)."""
    for root in (SENIOR_OUT / "data", SENIOR_OUT / "pilot_data"):
        if list((root / "redsea" / f"donor_id={donor}").glob("*.parquet")):
            return root
    raise SystemExit(f"no REDSEA output for donor {donor} under {SENIOR_OUT}/{{data,pilot_data}}")


def read_partition(root: Path, kind: str, donor: str, columns: list[str]) -> pd.DataFrame:
    import pyarrow.parquet as pq
    files = sorted((root / kind / f"donor_id={donor}").glob("*.parquet"))
    if not files:
        raise SystemExit(f"missing {root}/{kind}/donor_id={donor}")
    have = set(pq.ParquetFile(files[0]).schema.names)
    cols = [c for c in columns if c in have]
    return pd.concat([pq.read_table(f, columns=cols).to_pandas() for f in files], ignore_index=True)


# --------------------------------------------------------------------------------------------
# geometry helpers
# --------------------------------------------------------------------------------------------

def ring_stats(coords) -> tuple[float, float, float]:
    """Shoelace centroid and area of one ring (closed or open)."""
    a = np.asarray(coords, dtype=np.float64)
    if a.ndim != 2 or len(a) < 3:
        return np.nan, np.nan, 0.0
    x, y = a[:, 0], a[:, 1]
    x1, y1 = np.roll(x, -1), np.roll(y, -1)
    cr = x * y1 - x1 * y
    area = cr.sum() / 2.0
    if abs(area) < 1e-9:
        return float(x.mean()), float(y.mean()), 0.0
    return (float(((x + x1) * cr).sum() / (6 * area)),
            float(((y + y1) * cr).sum() / (6 * area)), float(abs(area)))


def geom_centroid(g) -> tuple[float, float, float]:
    """Centroid/area of a Polygon, or of the largest part of a MultiPolygon (exterior only)."""
    if not g:
        return np.nan, np.nan, 0.0
    t, c = g.get("type"), g.get("coordinates")
    if not c:
        return np.nan, np.nan, 0.0
    if t == "Polygon":
        return ring_stats(c[0])
    if t == "MultiPolygon":
        parts = [ring_stats(p[0]) for p in c if p]
        return max(parts, key=lambda s: s[2]) if parts else (np.nan, np.nan, 0.0)
    return np.nan, np.nan, 0.0


def grouped_ring_centroids(key: np.ndarray, x: np.ndarray, y: np.ndarray):
    """Vectorised shoelace over contiguous runs of ``key`` (one ring per run)."""
    n = len(key)
    starts = np.flatnonzero(np.r_[True, key[1:] != key[:-1]])
    ends = np.r_[starts[1:], n] - 1
    nxt = np.arange(1, n + 1)
    nxt[ends] = starts
    cr = x * y[nxt] - x[nxt] * y
    a2 = np.add.reduceat(cr, starts)
    cx = np.add.reduceat((x + x[nxt]) * cr, starts)
    cy = np.add.reduceat((y + y[nxt]) * cr, starts)
    with np.errstate(invalid="ignore", divide="ignore"):
        gx = np.where(np.abs(a2) > 1e-9, cx / (3 * a2), np.nan)
        gy = np.where(np.abs(a2) > 1e-9, cy / (3 * a2), np.nan)
    cnt = np.diff(np.r_[starts, n])
    mx = np.add.reduceat(x, starts) / cnt
    my = np.add.reduceat(y, starts) / cnt
    gx = np.where(np.isfinite(gx), gx, mx)
    gy = np.where(np.isfinite(gy), gy, my)
    return key[starts], gx, gy, np.abs(a2) / 2


# --------------------------------------------------------------------------------------------
# nuclei: v3 PhenoCycler nucleus centroids
# --------------------------------------------------------------------------------------------

def v3_nuclei_path(donor: str) -> Path:
    return OUT / "cache" / f"v3_nuclei_{donor}.parquet"


def extract_v3_nuclei(donor: str) -> pd.DataFrame:
    import ijson
    import pyarrow.parquet as pq

    ent = manifest_entry(donor)
    gj = Path(ent["cell_geojson"]["path"])
    pqf = Path(ent["measurement_csv"]["path"])
    log(f"{donor}: streaming {gj.name} ({gj.stat().st_size/1e9:.2f} GB)")
    t0 = time.time()
    ids, cx, cy, nx, ny, na = [], [], [], [], [], []
    with open(gj, "rb") as fh:
        for k, f in enumerate(ijson.items(fh, "features.item", use_float=True)):
            ids.append(str(f.get("id", "")))
            a, b, _ = geom_centroid(f.get("geometry"))
            c, d, e = geom_centroid(f.get("nucleusGeometry"))
            cx.append(a); cy.append(b); nx.append(c); ny.append(d); na.append(e)
            if (k + 1) % 250_000 == 0:
                log(f"  {k+1:,} features ({time.time()-t0:.0f}s)")
    geo = pd.DataFrame({"object_id": ids, "cell_px_x": cx, "cell_px_y": cy,
                        "nuc_px_x": nx, "nuc_px_y": ny, "nuc_area_px2": na})
    log(f"{donor}: {len(geo):,} features, {np.isfinite(geo.nuc_px_x).sum():,} with a nucleus "
        f"({time.time()-t0:.0f}s)")

    truth = pq.read_table(pqf, columns=["Object ID", "Centroid X µm", "Centroid Y µm",
                                        "Parent"]).to_pandas()
    truth.columns = ["object_id", "x_um", "y_um", "parent"]
    j = geo.merge(truth, on="object_id", how="inner")
    fit = {}
    for ax in ("x", "y"):
        ok = np.isfinite(j[f"cell_px_{ax}"]) & np.isfinite(j[f"{ax}_um"])
        s, b = np.polyfit(j.loc[ok, f"cell_px_{ax}"], j.loc[ok, f"{ax}_um"], 1)
        resid = np.abs(j.loc[ok, f"{ax}_um"] - (s * j.loc[ok, f"cell_px_{ax}"] + b))
        fit[ax] = {"um_per_px": float(s), "intercept_um": float(b),
                   "median_abs_resid_um": float(np.median(resid)),
                   "p99_abs_resid_um": float(np.percentile(resid, 99))}
        if not 0.4 < s < 0.6 or fit[ax]["median_abs_resid_um"] > 1.0:
            raise SystemExit(f"{donor}: px->um fit on {ax} is implausible: {fit[ax]}")
    sx, bx, sy, by = fit["x"]["um_per_px"], fit["x"]["intercept_um"], fit["y"]["um_per_px"], fit["y"]["intercept_um"]
    out = pd.DataFrame({
        "object_id": j["object_id"],
        "nuc_x_um": sx * j["nuc_px_x"] + bx, "nuc_y_um": sy * j["nuc_px_y"] + by,
        "nuc_area_um2": j["nuc_area_px2"] * sx * sy,
        "cell_x_um": j["x_um"], "cell_y_um": j["y_um"],
        "parent": j["parent"].astype(str)})
    meta = {"donor": donor, "geojson": str(gj), "parquet": str(pqf), "n_features": len(geo),
            "n_joined": len(j), "n_nucleus": int(np.isfinite(out.nuc_x_um).sum()),
            "px_to_um_fit": fit, "seconds": round(time.time() - t0, 1)}
    path = v3_nuclei_path(donor)
    path.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(path, index=False)
    path.with_suffix(".json").write_text(json.dumps(meta, indent=2) + "\n")
    log(f"{donor}: wrote {path} ({meta['n_nucleus']:,} nuclei); fit {fit}")
    return out


def load_v3_nuclei(donor: str) -> pd.DataFrame:
    p = v3_nuclei_path(donor)
    return pd.read_parquet(p) if p.exists() else extract_v3_nuclei(donor)


def placement_null(xy: np.ndarray, radii=(1.0, 2.0, 3.0, 3.5, 4.5, 5.0), offset_um=25.0) -> dict:
    """What a registration that is off by 25 um would score, from the v3 cloud itself.

    Approximates ``refine.score_transform``'s ``matched_within_null`` before any Xenium data
    exist: a displaced point finds a v3 nucleus within r at this rate by density alone.
    """
    from scipy.spatial import cKDTree
    tree = cKDTree(xy)
    d2, _ = tree.query(xy, k=2, workers=8)
    rows = []
    for dx, dy in ((offset_um, 0), (0, offset_um), (-offset_um, 0), (0, -offset_um)):
        d, _ = tree.query(xy + np.array([dx, dy]), k=1, workers=8)
        rows.append([(d <= r).mean() for r in radii])
    null = np.mean(rows, axis=0)
    nn = d2[:, 1]
    return {"offset_um": offset_um,
            "matched_within_null": {f"{r:g}um": round(float(v), 4) for r, v in zip(radii, null)},
            "self_nn_spacing_um": {"p10": float(np.percentile(nn, 10)),
                                   "p50": float(np.median(nn)), "mean": float(nn.mean())}}


def cmd_nuclei(args) -> int:
    for d in args.donor:
        donor = canon(d)
        nuc = extract_v3_nuclei(donor) if args.force or not v3_nuclei_path(donor).exists() \
            else load_v3_nuclei(donor)
        ok = np.isfinite(nuc.nuc_x_um)
        xy = nuc.loc[ok, ["nuc_x_um", "nuc_y_um"]].to_numpy()
        pn = placement_null(xy)
        pn["n_nuclei"] = int(ok.sum())
        pn["extent_um"] = {"x": [float(xy[:, 0].min()), float(xy[:, 0].max())],
                           "y": [float(xy[:, 1].min()), float(xy[:, 1].max())]}
        pn["by_parent"] = nuc.loc[ok, "parent"].value_counts().to_dict()
        (OUT / "cache" / f"placement_null_{donor}.json").write_text(json.dumps(pn, indent=2) + "\n")
        log(f"{donor}: placement null {pn['matched_within_null']}  nn {pn['self_nn_spacing_um']}")
    return 0


# --------------------------------------------------------------------------------------------
# protein: REDSEA vs raw
# --------------------------------------------------------------------------------------------

def load_protein(donor: str, markers: list[str]) -> tuple[pd.DataFrame, list[str], str]:
    """object_id-joined raw (QuPath whole-cell mean) and REDSEA ``Cell__`` values + QC."""
    root = senior_root(donor)
    cells = read_partition(root, "cells", donor, ["object_id", "cell_region"] + markers)
    red = read_partition(root, "redsea", donor, ["object_id"] + [f"Cell__{m}" for m in markers])
    qc = read_partition(root, "geometry_qc", donor, ["object_id", "analysis_eligible"])
    have = [m for m in markers if m in cells.columns and f"Cell__{m}" in red.columns]
    df = cells.set_index("object_id")[["cell_region"] + have].rename(
        columns={m: f"raw__{m}" for m in have})
    df = df.join(red.set_index("object_id").rename(
        columns={f"Cell__{m}": f"red__{m}" for m in have})[[f"red__{m}" for m in have]], how="left")
    df = df.join(qc.set_index("object_id")[["analysis_eligible"]], how="left")
    df["analysis_eligible"] = df["analysis_eligible"].fillna(False).astype(bool)
    return df, have, root.name


def reference_panel() -> dict[str, str]:
    if not REF_GENE_PANEL.exists():
        return {}
    gp = json.loads(REF_GENE_PANEL.read_text())
    return {t["type"]["data"]["name"]: t["source"]["category"]
            for t in gp["payload"]["targets"] if t["type"]["descriptor"] == "gene"}


def resolve_genes(marker: str, available: set[str]) -> tuple[str, list[str], list[str]]:
    """(kind, genes used, direct genes missing). kind 'none' when nothing is on the panel."""
    opts = MARKER_GENES.get(marker, [])
    missing = [g for k, gs in opts if k == "direct" for g in gs if g not in available]
    for kind, genes in opts:
        use = [g for g in genes if g in available]
        if use:
            return kind, use, missing
    return "none", [], missing


def cmd_protein(args) -> int:
    ref = reference_panel()
    base = {g for g, c in ref.items() if c == "base"}
    markers = sorted(MARKER_GENES)
    for d in args.donor:
        donor = canon(d)
        df, have, rootname = load_protein(donor, markers)
        rows = []
        elig = df["analysis_eligible"].to_numpy()
        for m in have:
            raw = df[f"raw__{m}"].to_numpy(float)
            red = df[f"red__{m}"].to_numpy(float)
            ok = elig & np.isfinite(raw) & np.isfinite(red)
            pos = ok & (raw > 0)
            z = pos & (red == 0)
            kind, genes, missing = resolve_genes(m, base)
            rows.append({"marker": m, "redsea_applied": m not in PASS_THROUGH,
                         "n_eligible": int(ok.sum()),
                         "frac_zeroed_of_raw_gt0": round(float(z.sum() / max(pos.sum(), 1)), 4),
                         "n_zeroed": int(z.sum()),
                         "median_raw": float(np.median(raw[ok])),
                         "median_redsea": float(np.median(red[ok])),
                         "rna_kind_stock5k": kind, "rna_genes_stock5k": ",".join(genes),
                         "direct_genes_absent_stock5k": ",".join(missing)})
        tab = pd.DataFrame(rows).sort_values("frac_zeroed_of_raw_gt0", ascending=False)
        out = OUT / donor
        out.mkdir(parents=True, exist_ok=True)
        tab.to_csv(out / "protein_zeroing.csv", index=False)
        log(f"{donor} ({rootname}): {len(df):,} cells, {int(elig.sum()):,} analysis-eligible; "
            f"wrote {out/'protein_zeroing.csv'}")
        with pd.option_context("display.width", 200, "display.max_rows", 100):
            print(tab.head(args.show).to_string(index=False))
    return 0


# --------------------------------------------------------------------------------------------
# Xenium side
# --------------------------------------------------------------------------------------------

def xenium_nuclei(bundle: Path) -> pd.DataFrame:
    """One nucleus centroid per Xenium cell (largest nucleus when a cell has several)."""
    import pyarrow.parquet as pq
    names = pq.ParquetFile(bundle / "nucleus_boundaries.parquet").schema.names
    cols = ["cell_id", "vertex_x", "vertex_y"] + (["label_id"] if "label_id" in names else [])
    nb = pq.read_table(bundle / "nucleus_boundaries.parquet", columns=cols).to_pandas()
    cid = nb["cell_id"].astype(str).to_numpy()
    key = nb["label_id"].to_numpy() if "label_id" in nb.columns else pd.factorize(cid)[0]
    k, gx, gy, area = grouped_ring_centroids(np.asarray(key), nb["vertex_x"].to_numpy(float),
                                             nb["vertex_y"].to_numpy(float))
    starts = np.flatnonzero(np.r_[True, np.asarray(key)[1:] != np.asarray(key)[:-1]])
    df = pd.DataFrame({"cell_id": cid[starts], "nx": gx, "ny": gy, "narea": area})
    df = df.sort_values("narea", ascending=False).drop_duplicates("cell_id")
    return df.reset_index(drop=True)


def xenium_counts(bundle: Path, genes: list[str]) -> tuple[pd.DataFrame, int, int]:
    """Counts for ``genes`` per cell from cell_feature_matrix.h5 (10x CSC: features x cells)."""
    import h5py
    with h5py.File(bundle / "cell_feature_matrix.h5", "r") as f:
        g = f["matrix"]
        names = g["features/name"][:].astype(str)
        ftype = g["features/feature_type"][:].astype(str)
        barcodes = g["barcodes"][:].astype(str)
        indptr = g["indptr"][:]
        indices = g["indices"][:]
        data = g["data"][:]
    is_gene = ftype == "Gene Expression"
    row = {n: i for i, n in enumerate(names) if is_gene[i]}
    want = [x for x in genes if x in row]
    out = np.zeros((len(barcodes), len(want)), dtype=np.int32)
    if want:
        rows = np.array([row[x] for x in want])
        lut = np.full(len(names), -1, dtype=np.int64)
        lut[rows] = np.arange(len(want))
        sel = np.flatnonzero(lut[indices] >= 0)
        col = np.searchsorted(indptr, sel, side="right") - 1
        np.add.at(out, (col, lut[indices[sel]]), data[sel].astype(np.int32))
    df = pd.DataFrame(out, columns=want)
    df.insert(0, "cell_id", barcodes)
    return df, int(is_gene.sum()), int((ftype == "Genomic Control").sum())


def gdna_lambda(cells: pd.DataFrame, n_probes: int) -> np.ndarray:
    """Per-cell expected gDNA counts per gene, shrunk toward the cell-area-decile mean."""
    gc = cells["genomic_control_counts"].to_numpy(float)
    dec = pd.qcut(cells["cell_area"].to_numpy(float), 10, labels=False, duplicates="drop")
    mu = pd.Series(gc).groupby(dec).transform("mean").to_numpy() / n_probes
    return (gc + n_probes * mu) / (2 * n_probes)


def rna_status(counts: np.ndarray, lam: np.ndarray, k: int) -> np.ndarray:
    """+1 RNA-positive, -1 RNA-negative, 0 ambiguous."""
    from scipy.stats import poisson
    p = poisson.sf(counts - 1, k * lam)
    s = np.zeros(len(counts), dtype=np.int8)
    s[(counts >= 2) & (p < RNA_POS_P)] = 1
    s[counts == 0] = -1
    return s


# --------------------------------------------------------------------------------------------
# registration
# --------------------------------------------------------------------------------------------

#: Coarse-to-fine block matching for the GLOBAL affine: (window, step, search, pixel) in um.
#: register_points.fit's rigid seed is on a 5-degree grid, and nearest-neighbour ICP cannot
#: recover a residual rotation: 1 degree is already ~100 um at the ends of a 22 mm section,
#: while ICP only sees offsets below half the nucleus spacing (~2.5 um).
GLOBAL_PASSES = ((1000.0, 500.0, 400.0, 8.0), (600.0, 300.0, 120.0, 5.0), (500.0, 500.0, 40.0, 2.0))
#: Local field pass: 6533 was physically sheared (kimwipe drag, up to ~90 um) between the runs.
FIELD_PASS = (500.0, 250.0, 80.0, 4.0)
BLOCK_MIN_Z = 6.0


def block_shifts(tf, moving, fixed, *, win_um, step_um, search_um, px_um, min_pts=150):
    """Per-window translation that lands the warped moving nuclei on the fixed nuclei.

    Phase correlation of nucleus-density rasters per window, so it uses tissue structure rather
    than nearest neighbours. Returns one row per scored window: its fixed-frame centre, the
    shift to ADD to the warped moving points, and the peak's z-score within the search box.
    """
    from scipy.ndimage import gaussian_filter
    from scipy.spatial import cKDTree
    cur = tf.apply(moving)
    lo = np.floor(np.minimum(cur.min(0), fixed.min(0)) / step_um) * step_um
    hi = np.maximum(cur.max(0), fixed.max(0))
    xs = np.arange(lo[0], hi[0] + step_um, step_um)
    ys = np.arange(lo[1], hi[1] + step_um, step_um)
    half, reach = win_um / 2, win_um / 2 + search_um
    n = int(np.ceil(2 * reach / px_um))
    lim = int(round(search_um / px_um))
    ctr = n // 2
    tm, tfx = cKDTree(cur), cKDTree(fixed)

    def rast(pts, x0, y0):
        ix = np.floor((pts[:, 0] - x0) / px_um).astype(int)
        iy = np.floor((pts[:, 1] - y0) / px_um).astype(int)
        ok = (ix >= 0) & (ix < n) & (iy >= 0) & (iy < n)
        a = np.zeros((n, n), np.float32)
        np.add.at(a, (iy[ok], ix[ok]), 1.0)
        a = gaussian_filter(a, 1.0)
        return a - a.mean()

    rows = []
    for cy in ys:
        for cx in xs:
            im = tm.query_ball_point([cx, cy], r=half * np.sqrt(2))
            if len(im) < min_pts:
                continue
            m = cur[im]
            m = m[(np.abs(m[:, 0] - cx) <= half) & (np.abs(m[:, 1] - cy) <= half)]
            jf = tfx.query_ball_point([cx, cy], r=reach * np.sqrt(2))
            if len(m) < min_pts or len(jf) < min_pts:
                continue
            f = fixed[jf]
            f = f[(np.abs(f[:, 0] - cx) <= reach) & (np.abs(f[:, 1] - cy) <= reach)]
            x0, y0 = cx - reach, cy - reach
            # Phase correlation: whitening the cross-power spectrum removes the broad density
            # texture that otherwise flattens the peak, leaving a sharp one at the shift.
            cp = np.fft.rfft2(rast(f, x0, y0)) * np.conj(np.fft.rfft2(rast(m, x0, y0)))
            cp /= np.abs(cp) + 1e-3 * np.abs(cp).max()
            cc = np.fft.fftshift(np.fft.irfft2(cp, s=(n, n)))
            sub = cc[ctr - lim:ctr + lim + 1, ctr - lim:ctr + lim + 1]
            sd = sub.std()
            if not sd > 0:
                continue
            k = np.unravel_index(np.argmax(sub), sub.shape)
            rows.append((cx, cy, (k[1] - lim) * px_um, (k[0] - lim) * px_um,
                         float((sub[k] - sub.mean()) / sd), len(m), len(f)))
    return pd.DataFrame(rows, columns=["cx", "cy", "dx", "dy", "z", "n_moving", "n_fixed"])


def fit_affine_robust(src, dst, *, k=3.0, floor_um=10.0, n_iter=10):
    """Least-squares 2x3 affine src->dst, iteratively dropping residuals > max(k*MAD, floor)."""
    X = np.column_stack([src, np.ones(len(src))])
    keep = np.ones(len(src), bool)
    for _ in range(n_iter):
        sol, *_ = np.linalg.lstsq(X[keep], dst[keep], rcond=None)
        r = np.linalg.norm(dst - X @ sol, axis=1)
        new = r <= max(k * 1.4826 * np.median(r[keep]), floor_um)
        if new.sum() < 6 or (new == keep).all():
            break
        keep = new
    return sol.T, keep, r


def fit_global(moving, fixed):
    """Rigid seed (register_points.fit) -> block-matched affine passes -> similarity ICP."""
    rp, _rf, Transform, _G = senior_imports()
    tf, prov = rp.fit(moving, fixed, verbose=False)
    passes = []
    for win, step, search, px in GLOBAL_PASSES:
        sh = block_shifts(tf, moving, fixed, win_um=win, step_um=step, search_um=search, px_um=px)
        good = sh[sh.z >= BLOCK_MIN_Z]
        rec = {"win_um": win, "search_um": search, "px_um": px, "windows": int(len(sh)),
               "windows_z_ok": int(len(good))}
        if len(good) >= 6:
            src = good[["cx", "cy"]].to_numpy()
            U, keep, r = fit_affine_robust(src, src + good[["dx", "dy"]].to_numpy())
            L = U[:, :2] @ tf.affine[:, :2]
            t = U[:, :2] @ tf.affine[:, 2] + U[:, 2]
            tf = Transform(affine=np.column_stack([L, t]))
            rec.update({"inliers": int(keep.sum()), "inlier_resid_p50_um": float(np.median(r[keep])),
                        "rotation_deg": tf.rotation_deg, "scale": tf.scale})
        passes.append(rec)
    for gate in (5.0, 3.0):
        tf = rp.icp(moving, fixed, initial=tf, max_dist_um=gate)
    return tf, {"register_points_seed": prov, "block_affine_passes": passes,
                "icp_gates_um": [5.0, 3.0], "rotation_deg": tf.rotation_deg, "scale": tf.scale,
                "translation_um": list(tf.translation_um), "mirrored": tf.is_mirrored}


def pair_field(tf, moving, fixed, *, gates=(6.0, 4.5, 3.5), pixel_um=250.0, min_pairs=20):
    """Field refinement from mutual nucleus pairs, re-estimated relative to the affine."""
    from scipy.ndimage import distance_transform_edt, gaussian_filter
    from scipy.stats import binned_statistic_2d
    _rp, rf, Transform, GridSpec = senior_imports()
    aff = Transform(affine=tf.affine)
    base = aff.apply(moving)
    x0, y0 = base.min(0) - pixel_um
    x1, y1 = base.max(0) + pixel_um
    w = int(np.ceil((x1 - x0) / pixel_um)) + 1
    h = int(np.ceil((y1 - y0) / pixel_um)) + 1
    xe = x0 - 0.5 * pixel_um + np.arange(w + 1) * pixel_um
    ye = y0 - 0.5 * pixel_um + np.arange(h + 1) * pixel_um
    grid = GridSpec(float(x0), float(y0), float(pixel_um), h, w)
    hist = []
    for gate in gates:
        i, j, _d = rf.mutual_pairs(tf.apply(moving), fixed, max_dist_um=gate)
        pts, res = base[i], fixed[j] - base[i]
        cnt = binned_statistic_2d(pts[:, 1], pts[:, 0], None, "count", bins=[ye, xe])[0]
        pop = cnt >= min_pairs
        if not pop.any():
            break
        _dd, idx = distance_transform_edt(~pop, return_indices=True)
        field = np.zeros((2, h, w), np.float32)
        for q in (0, 1):
            st = binned_statistic_2d(pts[:, 1], pts[:, 0], res[:, q], "median", bins=[ye, xe])[0]
            st = np.where(np.isfinite(st), st, 0.0)
            field[q] = gaussian_filter(np.where(pop, st, st[tuple(idx)]), 1.0)
        tf = Transform(affine=aff.affine, displacement=field, grid=grid,
                       meta={"stage": "pair_field", "pixel_um": pixel_um})
        hist.append({"gate_um": gate, "n_pairs": int(i.size), "bins_populated": int(pop.sum())})
    return tf, hist


def fit_local(tf_global, moving, fixed):
    """Block-matched displacement field on top of the global affine, then pair refinement."""
    import warnings
    from scipy.ndimage import distance_transform_edt, gaussian_filter, generic_filter
    _rp, _rf, Transform, GridSpec = senior_imports()
    win, step, search, px = FIELD_PASS
    aff = Transform(affine=tf_global.affine)
    sh = block_shifts(aff, moving, fixed, win_um=win, step_um=step, search_um=search, px_um=px)
    meta = {"windows": int(len(sh)), "windows_z_ok": int((sh.z >= BLOCK_MIN_Z).sum())}
    if meta["windows_z_ok"] < 4:
        tf, hist = pair_field(aff, moving, fixed)
        return tf, {"block_field": meta, "pair_field": hist}
    x0, y0 = sh.cx.min(), sh.cy.min()
    ix = np.round((sh.cx - x0) / step).astype(int).to_numpy()
    iy = np.round((sh.cy - y0) / step).astype(int).to_numpy()
    H, W = iy.max() + 1, ix.max() + 1
    shift = np.full((2, H, W), np.nan)
    good = np.zeros((H, W), bool)
    ok = (sh.z >= BLOCK_MIN_Z).to_numpy()
    shift[0, iy[ok], ix[ok]] = sh.dx.to_numpy()[ok]
    shift[1, iy[ok], ix[ok]] = sh.dy.to_numpy()[ok]
    good[iy[ok], ix[ok]] = True
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        for q in (0, 1):
            med = generic_filter(shift[q], np.nanmedian, size=3, mode="constant", cval=np.nan)
            good &= ~(np.isfinite(med) & (np.abs(shift[q] - med) > 3 * px))
    meta["windows_kept"] = int(good.sum())
    meta["max_abs_shift_um"] = float(np.nanmax(np.abs(shift[:, good]))) if good.any() else 0.0
    field = np.zeros((2, H, W), np.float32)
    _d, idx = distance_transform_edt(~good, return_indices=True)
    for q in (0, 1):
        v = np.where(good, shift[q], 0.0)
        field[q] = gaussian_filter(np.where(good, v, v[tuple(idx)]), 1.0)
    tf = Transform(affine=aff.affine, displacement=field,
                   grid=GridSpec(float(x0), float(y0), float(step), int(H), int(W)),
                   meta={"stage": "block_field"})
    tf, hist = pair_field(tf, moving, fixed)
    return tf, {"block_field": meta, "pair_field": hist}


def checkerboard_heldout(fit_fn, moving, fixed, *, square_um=500.0) -> dict:
    """Fit on one checkerboard colour, score on the other (both ways)."""
    _rp, rf, _T, _G = senior_imports()
    colour = ((np.floor(moving[:, 0] / square_um) + np.floor(moving[:, 1] / square_um)) % 2
              ).astype(bool)
    out = {}
    for name, fit_mask in (("fit_black_score_white", ~colour), ("fit_white_score_black", colour)):
        tf = fit_fn(moving[fit_mask])
        s = rf.score_transform(tf, moving[~fit_mask], fixed)
        out[name] = {"pair_p50_um": s["pair_distance_um"].get("p50"),
                     "pair_p90_um": s["pair_distance_um"].get("p90"),
                     "matched_within": s["matched_within"],
                     "matched_within_null": s["matched_within_null"]}
    p50s = [v["pair_p50_um"] for v in out.values() if isinstance(v, dict)]
    out["heldout_pair_p50_um"] = (None if any(p is None for p in p50s)
                                  else float(max(p50s)))
    return out


def register(moving, fixed, *, transform_dir=None, transform_stem="transform", out_dir=None):
    rp, rf, Transform, _G = senior_imports()
    rec: dict = {}
    if transform_dir:
        tf = Transform.load(Path(transform_dir), stem=transform_stem)
        rec["source"] = f"loaded {transform_dir}/{transform_stem}"
        rec["score"] = rf.score_transform(tf, moving, fixed)
        # A loaded transform was fitted on other data (DAPI images / v2 nuclei), so the v3
        # nuclei are an independent test set; its pair p50 stands in for the held-out median.
        rec["heldout_pair_p50_um"] = rec["score"]["pair_distance_um"].get("p50")
        return tf, rec
    t0 = time.time()
    tf_g, gmeta = fit_global(moving, fixed)
    sg = rf.score_transform(tf_g, moving, fixed)
    rec.update({"source": "refit: register_points seed -> block-matched affine -> ICP",
                "global": {"meta": gmeta, "score": sg,
                           "held_out_split_half": rp.held_out(moving, fixed, tf_g)}})
    log(f"global fit: rot {gmeta['rotation_deg']:.3f} deg, pair p50 "
        f"{sg['pair_distance_um'].get('p50')}  matched_within {sg['matched_within']}  "
        f"({time.time()-t0:.0f}s)")
    tf_l, lmeta = fit_local(tf_g, moving, fixed)
    sl = rf.score_transform(tf_l, moving, fixed)
    rec["local"] = {"meta": lmeta, "score": sl}
    log(f"local field: pair p50 {sl['pair_distance_um'].get('p50')}  "
        f"matched_within {sl['matched_within']}  ({time.time()-t0:.0f}s)")
    # Same parsimony rule as register_points: the field only when it clearly helps.
    p_g = sg["pair_distance_um"].get("p50", 9e9)
    p_l = sl["pair_distance_um"].get("p50", 9e9)
    use_local = (p_l < p_g - 0.05 or sl["matched_within"]["2um"] > sg["matched_within"]["2um"] + 0.02)
    if use_local:
        def refit(mv):
            return fit_local(rp.icp(mv, fixed, initial=tf_g, max_dist_um=3.0), mv, fixed)[0]
        rec["held_out_checkerboard"] = checkerboard_heldout(refit, moving, fixed)
        rec["heldout_pair_p50_um"] = rec["held_out_checkerboard"]["heldout_pair_p50_um"]
        tf, rec["chosen"] = tf_l, "global+local_field"
    else:
        ho = rec["global"]["held_out_split_half"]
        vals = [v for k, v in ho.items() if k.startswith("fit_") and v is not None]
        rec["heldout_pair_p50_um"] = float(max(vals)) if vals else None
        tf, rec["chosen"] = tf_g, "global"
    rec["score"] = rf.score_transform(tf, moving, fixed)
    if out_dir is not None:
        tf.save(Path(out_dir), stem="xenium_to_v3")
    return tf, rec

# --------------------------------------------------------------------------------------------
# agreement metrics
# --------------------------------------------------------------------------------------------

def auc(vals: np.ndarray, pos: np.ndarray) -> float:
    from scipy.stats import rankdata
    n1 = int(pos.sum())
    n0 = len(vals) - n1
    if n1 == 0 or n0 == 0:
        return np.nan
    r = rankdata(vals)
    return float((r[pos].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def marker_metrics(raw, red, status, fold) -> dict:
    """Agreement of raw vs REDSEA protein with RNA on one marker's scorable cells."""
    known = status != 0
    raw, red, st, fd = raw[known], red[known], status[known], fold[known]
    pos = st == 1
    res = {"n_scored": int(known.sum()), "n_rna_pos": int(pos.sum()),
           "rna_pos_rate": float(pos.mean()) if len(pos) else np.nan}
    if pos.sum() < 20 or (~pos).sum() < 20:
        return res
    a_raw, a_red = auc(raw, pos), auc(red, pos)
    jk = []
    for k in np.unique(fd):
        keep = fd != k
        jk.append(auc(red[keep], pos[keep]) - auc(raw[keep], pos[keep]))
    jk = np.array(jk)
    K = len(jk)
    se = float(np.sqrt((K - 1) / K * np.sum((jk - jk.mean()) ** 2))) if K > 1 else np.nan
    res.update({"auc_raw": a_raw, "auc_redsea": a_red, "delta_auc": a_red - a_raw,
                "delta_auc_jackknife_se": se,
                "delta_auc_z": (a_red - a_raw) / se if se and se > 0 else np.nan})
    rawpos = raw > 0
    Z = rawpos & (red == 0)
    Kp = rawpos & (red > 0)
    res.update({"n_zeroed": int(Z.sum()), "frac_zeroed": float(Z.mean()),
                "rna_pos_rate_zeroed": float(pos[Z].mean()) if Z.any() else np.nan,
                "zero_rate_rna_pos": float(Z[pos & rawpos].mean()) if (pos & rawpos).any() else np.nan,
                "zero_rate_rna_neg": float(Z[~pos & rawpos].mean()) if (~pos & rawpos).any() else np.nan,
                "n_rna_pos_zeroed": int((Z & pos).sum())})
    lo_q, hi_q = np.quantile(raw, [0.1, 0.9])
    p_bottom = float(pos[raw <= lo_q].mean())
    res["rna_pos_rate_raw_bottom10"] = p_bottom
    res["rna_pos_rate_raw_top10"] = float(pos[raw >= hi_q].mean())
    if Z.sum() >= 50:
        # 20 quantile bins of the zeroed cells' raw intensity; kept cells are confined to the
        # same range, so no open-ended edge bin can import the bright (or dim) tail.
        nb = 20
        qs = np.quantile(raw[Z], np.linspace(0, 1, nb + 1))
        inK = Kp & (raw >= qs[0]) & (raw <= qs[-1])
        bz = np.clip(np.searchsorted(qs, raw[Z], side="right") - 1, 0, nb - 1)
        bk = np.clip(np.searchsorted(qs, raw[inK], side="right") - 1, 0, nb - 1)
        pz, pk, wz = [], [], []
        for b in range(nb):
            nk = (bk == b).sum()
            nz = (bz == b).sum()
            if nk >= 20 and nz > 0:
                pz.append(pos[Z][bz == b].mean())
                pk.append(pos[inK][bk == b].mean())
                wz.append(nz)
        if wz:
            wz = np.array(wz, float)
            pzm = float(np.dot(wz, pz) / wz.sum())
            pkm = float(np.dot(wz, pk) / wz.sum())
            res["rna_pos_rate_zeroed_matched"] = pzm
            res["rna_pos_rate_kept_matched"] = pkm
            # 0 = zeroed cells look like the raw-dim background (REDSEA removed spillover);
            # 1 = as RNA-positive as kept cells of the same raw intensity (real signal lost).
            res["zeroed_rna_index"] = ((pzm - p_bottom) / (pkm - p_bottom)
                                       if pkm - p_bottom > 1e-3 else np.nan)
    return res


def cmd_run(args) -> int:
    donor = canon(args.donor)
    bundle = Path(args.bundle)
    need = ["cells.parquet", "nucleus_boundaries.parquet", "cell_feature_matrix.h5"]
    miss = [x for x in need if not (bundle / x).is_file()]
    if miss:
        raise SystemExit(f"bundle {bundle} is missing {miss}. The same-slide bundle for {donor} is "
                         f"on HiPerGator: /blue/maigan/smith6jt/Xenium_Analysis/data/"
                         f"{DONORS[donor]['bundle']}")
    out = Path(args.out) / donor if args.out else OUT / donor
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    nuc = load_v3_nuclei(donor)
    nuc = nuc[np.isfinite(nuc.nuc_x_um)].reset_index(drop=True)
    fixed = nuc[["nuc_x_um", "nuc_y_um"]].to_numpy()
    xn = xenium_nuclei(bundle)
    moving = xn[["nx", "ny"]].to_numpy()
    log(f"{donor}: {len(moving):,} Xenium nuclei -> {len(fixed):,} v3 nuclei")

    tf, reg = register(moving, fixed, transform_dir=args.transform_dir,
                       transform_stem=args.transform_stem, out_dir=out / "transform")
    ho = reg.get("heldout_pair_p50_um")
    reg["cell_gate_um"] = CELL_GATE_UM
    reg["cell_gate_pass"] = bool(ho is not None and ho <= CELL_GATE_UM)
    (out / "registration.json").write_text(json.dumps(reg, indent=2, default=float) + "\n")
    log(f"registration: {reg['source']} | held-out pair p50 {ho} um | gate "
        f"{'PASS' if reg['cell_gate_pass'] else 'FAIL'} (<= {CELL_GATE_UM} um)")
    if not reg["cell_gate_pass"] and not args.force:
        raise SystemExit("registration does not support a per-cell statement; stopping "
                         "(use --force to score anyway, and do not quote the result)")

    _rp, rf, _T, _G = senior_imports()
    warped = tf.apply(moving)
    i, j, dist = rf.mutual_pairs(warped, fixed, max_dist_um=args.cap_um)
    ni, _nj, _nd = rf.mutual_pairs(warped + np.array([25.0, 0.0]), fixed, max_dist_um=args.cap_um)
    match = {"cap_um": args.cap_um, "n_pairs": int(i.size),
             "frac_of_xenium_nuclei": float(i.size / len(moving)),
             "frac_of_v3_nuclei": float(i.size / len(fixed)),
             "frac_of_xenium_nuclei_null_25um": float(ni.size / len(moving)),
             "pair_distance_um": {q: float(np.percentile(dist, p)) for q, p in
                                  (("p50", 50), ("p90", 90), ("p99", 99))},
             "n_pairs_within_2um": int((dist <= 2.0).sum())}
    log(f"match: {match}")
    pairs = pd.DataFrame({"cell_id": xn["cell_id"].to_numpy()[i],
                          "object_id": nuc["object_id"].to_numpy()[j], "dist_um": dist,
                          "x_um": fixed[j, 0], "y_um": fixed[j, 1],
                          "parent": nuc["parent"].to_numpy()[j]})

    # RNA
    markers = sorted(MARKER_GENES)
    all_genes = sorted({g for m in markers for _k, gs in MARKER_GENES[m] for g in gs})
    counts, n_genes, n_probes = xenium_counts(bundle, all_genes)
    available = set(counts.columns) - {"cell_id"}
    xc = pd.read_parquet(bundle / "cells.parquet",
                         columns=["cell_id", "transcript_counts", "genomic_control_counts",
                                  "cell_area"])
    xc["cell_id"] = xc["cell_id"].astype(str)
    n_probes = n_probes or 21
    xc["gdna_lambda"] = gdna_lambda(xc, n_probes)
    xc["corrected_counts"] = xc["transcript_counts"] - xc["gdna_lambda"] * n_genes
    rna = xc.merge(counts, on="cell_id", how="left")
    df = pairs.merge(rna, on="cell_id", how="left")

    prot, have, rootname = load_protein(donor, markers)
    df = df.merge(prot, left_on="object_id", right_index=True, how="left")
    base_ok = (df["analysis_eligible"].fillna(False).to_numpy(bool)
               & (df["corrected_counts"].to_numpy(float) >= MIN_CORRECTED_COUNTS))
    bx = np.floor(df["x_um"].to_numpy() / FOLD_BLOCK_UM).astype(int)
    by = np.floor(df["y_um"].to_numpy() / FOLD_BLOCK_UM).astype(int)
    fold = (bx * 7 + by * 13) % N_FOLDS
    lam = df["gdna_lambda"].to_numpy(float)

    dist_all = df["dist_um"].to_numpy(float)
    for label, sub in (("", np.ones(len(df), bool)),
                       (f"_le{args.strict_um:g}um", dist_all <= args.strict_um)):
        rows = []
        for m in have:
            kind, genes, missing = resolve_genes(m, available)
            row = {"marker": m, "redsea_applied": m not in PASS_THROUGH, "rna_kind": kind,
                   "rna_genes": ",".join(genes), "direct_genes_absent": ",".join(missing)}
            if kind == "none":
                rows.append(row)
                continue
            cnt = df[genes].fillna(0).to_numpy().sum(axis=1)
            st = rna_status(cnt, lam, len(genes))
            raw = df[f"raw__{m}"].to_numpy(float)
            red = df[f"red__{m}"].to_numpy(float)
            ok = sub & base_ok & np.isfinite(raw) & np.isfinite(red)
            row.update(marker_metrics(raw[ok], red[ok], st[ok], fold[ok]))
            rows.append(row)
            if not label:
                log(f"  {m:16s} {kind:6s} {row.get('rna_genes',''):22s} "
                    f"dAUC {row.get('delta_auc', np.nan):+.4f} (z {row.get('delta_auc_z', np.nan):+.1f})"
                    f"  zeroed RNA+ {row.get('rna_pos_rate_zeroed', np.nan):.3f} vs kept-matched "
                    f"{row.get('rna_pos_rate_kept_matched', np.nan):.3f}")
        pd.DataFrame(rows).to_csv(out / f"marker_agreement{label}.csv", index=False)
    pairs.to_parquet(out / "matches.parquet", index=False)
    summary = {"donor": donor, "bundle": str(bundle), "senior_root": rootname,
               "registration_gate_pass": reg["cell_gate_pass"],
               "heldout_pair_p50_um": ho, "match": match,
               "n_scorable_cells": int(base_ok.sum()),
               "n_gene_features": n_genes, "n_genomic_control_probes": n_probes,
               "rna_positive_rule": f"counts>=2 and Poisson(k*lambda) p<{RNA_POS_P}; negative = 0",
               "min_corrected_counts": MIN_CORRECTED_COUNTS, "seconds": round(time.time() - t0, 1)}
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=float) + "\n")
    log(f"wrote {out}/marker_agreement.csv, matches.parquet, registration.json, summary.json")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("nuclei", help="v3 nucleus centroids + placement null")
    a.add_argument("--donor", nargs="+", default=["6533", "0112"])
    a.add_argument("--force", action="store_true")
    a.set_defaults(fn=cmd_nuclei)
    b = sub.add_parser("protein", help="REDSEA zeroing per marker + stock-5K gene availability")
    b.add_argument("--donor", nargs="+", default=["6533", "0112"])
    b.add_argument("--show", type=int, default=25)
    b.set_defaults(fn=cmd_protein)
    c = sub.add_parser("run", help="register, match and score one donor")
    c.add_argument("--donor", required=True)
    c.add_argument("--bundle", required=True, help="same-slide Xenium bundle directory")
    c.add_argument("--transform-dir", default=None)
    c.add_argument("--transform-stem", default="transform")
    c.add_argument("--cap-um", type=float, default=CAP_UM)
    c.add_argument("--strict-um", type=float, default=2.0,
                   help="also score the pairs within this distance (sensitivity)")
    c.add_argument("--out", default=None)
    c.add_argument("--force", action="store_true")
    c.set_defaults(fn=cmd_run)
    args = ap.parse_args(argv)
    return args.fn(args)


if __name__ == "__main__":
    raise SystemExit(main())
