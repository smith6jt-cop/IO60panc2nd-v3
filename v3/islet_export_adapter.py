#!/usr/bin/env python3
"""islet_export_adapter.py -- hand the v3 QuPath exports to the Senior typing lane in the shape it reads.

The typing lane (Islet-Explorer-Senior, `core/phenocycler/scripts/typing`) reads two things per image
that pass J does not write in its shape. This script derives both, per image, from the exports; it
decides nothing about cells and changes nothing in QuPath.

1. Islet outlines  ->  <export>/typing_islets/Islets__<tag>.geojson
     the lane knows ONE kind of islet object: a feature named `Islet_<N>`, classification `Islet`
     (`islet_membership.outlines`, `build_islet_geometry.islet_outlines`). The reviewed QuPath objects
     are Islet / EndocrineCluster / EndocrineSingle, named Islet_N / Cluster_N / Single_N with one
     number space per class. Which of them count as islets is the maintainer's ruling (--islets):
         all     every endocrine object, numbered  Islet_N -> N, Cluster_N -> 100000 + N,
                 Single_N -> 200000 + N  (so a number still says which QuPath object it is)
         islet   only the >= 35 um Islet objects, numbered N
     The size class (step 1b's designation: Islet / EndocrineCluster / EndocrineSingle), the QuPath name and id,
     the count of whole cells and the diameter travel as properties.
     With --expansion outlines each islet also gets `Islet_<N>_exp20um` (classification IsletExpanded):
     the 20 um zone outside it, minus every islet, cut where another islet is nearer -- the nearest-
     islet split the maintainer's 4_islet_expander.groovy made, computed here because no band is
     drawn in QuPath any more (2026-09-29).

2. The revised-export table  ->  <senior>/data/RevisedQuPathExports/<donor>.tsv
     the eleven columns `rebuild_donor_from_revised_export.py` reads, from cells ⋈ context on Object ID:
     Parent is rewritten to the islet names of (1) (a cell in an object that does not count as an
     islet is a Tissue cell); the signed distance to "Islet" is the distance to the nearest object
     that counts as one; centroids and distances are rescaled from QuPath's pixel size to the Senior
     manifest's, so a centroid lands on the pixel it came from (QuPath wrote um at the project
     entry's 0.5078 um/px, the manifest carries the qptiff's own 0.5077663810243286).

    islet_export_adapter.py --export <results/v3/export> --senior <worktree> --islets all|islet \\
        --expansion outlines|none [--images <tag> ...]
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import time
from pathlib import Path

import shapely
from shapely import STRtree
from shapely.geometry import MultiPoint, mapping, shape
from shapely.ops import unary_union

PROJECT = Path("/home/smith6jt/IO60panc2nd/project.qpproj")
PREFIX = {"Islet": ("Islet", 0), "EndocrineCluster": ("Cluster", 100_000), "EndocrineSingle": ("Single", 200_000)}
ISLET_CLASSES = {"all": ("Islet", "EndocrineCluster", "EndocrineSingle"), "islet": ("Islet",)}
EXPANSION_UM = 20.0
EXPANSION_SUFFIX = "_exp20um"
DISTANCE = "Signed distance to annotation with {} µm"
TSV_COLUMNS = ["Object ID", "Image", "Parent", "Classification", "Centroid X µm", "Centroid Y µm",
               "Nucleus: Area µm^2", "Cell: Area µm^2", "Cell: Solidity", "Nucleus/Cell area ratio",
               DISTANCE.format("Islet")]


def tag_of(image: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", image)


def class_of(feature: dict) -> str:
    return ((feature.get("properties") or {}).get("classification") or {}).get("name") or "Unclassified"


def typing_number(qupath_class: str, name: str) -> int:
    """Islet_12 -> 12, Cluster_7 -> 100007, Single_3 -> 200003; anything else is refused."""
    prefix, offset = PREFIX[qupath_class]
    match = re.fullmatch(rf"{prefix}_(\d+)", name or "")
    if not match:
        raise ValueError(f"{qupath_class} object named {name!r}: expected {prefix}_<N> (run step 4 after the import)")
    number = int(match.group(1))
    if not 0 < number < 100_000:
        raise ValueError(f"{name}: number outside 1..99999")
    return offset + number


def select_islets(features: list[dict], rule: str) -> dict[int, dict]:
    """typing number -> {geom, qupath_name, qupath_class, qupath_id, measurements} for the objects that count."""
    chosen: dict[int, dict] = {}
    for feature in features:
        qclass = class_of(feature)
        if qclass not in ISLET_CLASSES[rule]:
            continue
        props = feature.get("properties") or {}
        number = typing_number(qclass, props.get("name"))
        if number in chosen:
            raise ValueError(f"two objects resolve to Islet_{number}")
        geom = shape(feature["geometry"])
        if not geom.is_valid:
            geom = shapely.make_valid(geom)
        chosen[number] = {"geom": geom, "qupath_name": props.get("name"), "qupath_class": qclass,
                          "qupath_id": feature.get("id"), "measurements": props.get("measurements") or {}}
    return chosen


def nearest_islet_bands(cores: dict[int, object], width_px: float, spacing_px: float = 2.0) -> dict[int, object]:
    """number -> its expansion: within `width_px` of the islet, outside every islet, and nearer to this
    islet than to any other (a Voronoi partition of points sampled every `spacing_px` along the outlines).

    An islet whose 20 um zone meets no other islet's keeps the plain zone, so the partition is only
    computed where it changes something.
    """
    numbers = list(cores)
    geoms = [cores[n] for n in numbers]
    zones = [g.buffer(width_px) for g in geoms]
    tree = STRtree(zones)
    left, right = tree.query(zones, predicate="intersects")
    contested = {int(i) for i, j in zip(left, right) if i != j}
    all_cores = unary_union(geoms)
    bands: dict[int, object] = {}
    cells_of: dict[int, list] = {}
    if contested:
        points, owner = [], []
        for i in sorted(contested):
            boundary = geoms[i].boundary
            n = max(8, int(boundary.length / spacing_px))
            for k in range(n):
                p = boundary.interpolate(k / n, normalized=True)
                points.append((p.x, p.y))
                owner.append(i)
        # duplicate sample points (two islets sharing a vertex) would make the cells ambiguous
        seen, keep_pts, keep_owner = set(), [], []
        for pt, o in zip(points, owner):
            if pt not in seen:
                seen.add(pt)
                keep_pts.append(pt)
                keep_owner.append(o)
        hull = unary_union([zones[i] for i in contested]).envelope.buffer(width_px)
        cells = shapely.voronoi_polygons(MultiPoint(keep_pts), extend_to=hull, ordered=True)
        for cell, o in zip(cells.geoms, keep_owner):
            cells_of.setdefault(o, []).append(cell)
    for i, number in enumerate(numbers):
        band = zones[i].difference(all_cores)
        if i in contested:
            band = band.intersection(unary_union(cells_of[i]))
        if not band.is_valid:
            band = shapely.make_valid(band)
        # keep area only (a cut can leave a stray line)
        polys = [p for p in (band.geoms if hasattr(band, "geoms") else [band])
                 if p.geom_type in ("Polygon", "MultiPolygon") and not p.is_empty]
        bands[number] = unary_union(polys) if polys else None
    return bands


def islet_features(chosen: dict[int, dict], bands: dict[int, object] | None) -> list[dict]:
    out = []
    for number in sorted(chosen):
        o = chosen[number]
        m = o["measurements"]
        out.append({"type": "Feature", "id": o["qupath_id"], "geometry": mapping(o["geom"]),
                    "properties": {"objectType": "annotation", "name": f"Islet_{number}",
                                   "classification": {"name": "Islet"},
                                   "qupath_name": o["qupath_name"], "size_class": o["qupath_class"],
                                   # step 1b's own measurements (2026-10-02): the whole cells in the object and
                                   # how far across it is; the retired detector's keys are no longer written
                                   "n_cells": m.get("Whole cells"), "mean_diameter_um": m.get("Mean caliper diameter µm")}})
        band = (bands or {}).get(number)
        if band is not None and not band.is_empty:
            out.append({"type": "Feature", "geometry": mapping(band),
                        "properties": {"objectType": "annotation", "name": f"Islet_{number}{EXPANSION_SUFFIX}",
                                       "classification": {"name": "IsletExpanded"}}})
    return out


def write_tsv(cells: Path, context: Path, out: Path, chosen: dict[int, dict], rule: str, scale: float) -> dict:
    """cells ⋈ context -> the eleven columns, Parent and distance expressed in the islets that count."""
    import duckdb
    con = duckdb.connect()
    con.execute("CREATE TABLE names (qupath_name VARCHAR PRIMARY KEY, typing_name VARCHAR)")
    con.executemany("INSERT INTO names VALUES (?, ?)",
                    [(o["qupath_name"], f"Islet_{n}") for n, o in chosen.items()])
    have = {r[0] for r in con.execute("DESCRIBE SELECT * FROM read_parquet(?)", [str(context)]).fetchall()}
    wanted = [DISTANCE.format(c) for c in ISLET_CLASSES[rule]]
    missing = [c for c in wanted if c not in have]
    if missing:
        raise KeyError(f"{context.name}: no {missing} -- step 3 (signed distance) has not run for these classes")
    quoted = [f'x."{c}"' for c in wanted]
    # nearest object that counts as an islet: the smallest signed distance (inside one object it is negative)
    distance = quoted[0] if len(quoted) == 1 else "least(" + ", ".join(f"nullif({q}, 'nan'::DOUBLE)" for q in quoted) + ")"
    endocrine_name = r"^(Islet|Cluster|Single)_[0-9]+$"
    sql = f"""
        SELECT c."Object ID", c."Image",
               coalesce(n.typing_name,
                        CASE WHEN regexp_matches(x."Parent", '{endocrine_name}') THEN 'Annotation (Tissue)'
                             ELSE x."Parent" END) AS "Parent",
               coalesce(c."Classification", '') AS "Classification",
               x."Centroid X µm" * {scale!r} AS "Centroid X µm", x."Centroid Y µm" * {scale!r} AS "Centroid Y µm",
               c."Nucleus: Area µm^2", c."Cell: Area µm^2", c."Cell: Solidity", c."Nucleus/Cell area ratio",
               ({distance}) * {scale!r} AS "{DISTANCE.format('Islet')}"
        FROM read_parquet('{cells}') c
        JOIN read_parquet('{context}') x USING ("Object ID")
        LEFT JOIN names n ON n.qupath_name = x."Parent"
    """
    n_cells = con.execute(f"SELECT count(*) FROM read_parquet('{cells}')").fetchone()[0]
    n_context = con.execute(f"SELECT count(*), count(DISTINCT \"Object ID\") FROM read_parquet('{context}')").fetchone()
    part = out.with_suffix(".tsv.part")
    con.execute(f"COPY ({sql}) TO '{part}' (HEADER, DELIMITER '\t')")
    stats = con.execute(f"""
        SELECT count(*), count(*) FILTER (WHERE "Parent" LIKE 'Islet\\_%' ESCAPE '\\'),
               count(*) FILTER (WHERE "Parent" = 'Annotation (pLN)'),
               count(*) FILTER (WHERE "{DISTANCE.format('Islet')}" <= 0),
               count(*) FILTER (WHERE "{DISTANCE.format('Islet')}" IS NULL)
        FROM read_csv('{part}', delim='\t', header=true)""").fetchone()
    if not (n_cells == n_context[0] == n_context[1] == stats[0]):
        part.unlink()
        raise ValueError(f"{out.name}: cells {n_cells:,}, context {n_context[0]:,} ({n_context[1]:,} distinct ids), "
                         f"joined {stats[0]:,} -- the two exports are not of the same cells")
    part.replace(out)
    return {"rows": stats[0], "parent_is_islet": stats[1], "parent_is_pln": stats[2],
            "distance_inside_islet": stats[3], "distance_missing": stats[4]}


def one_image(entry: dict, export: Path, senior: Path, manifest: dict, rule: str, expansion: str) -> dict:
    name = entry["imageName"]
    tag = tag_of(name)
    row = manifest.get(name)
    if row is None:
        raise KeyError(f"{name}: not in the Senior manifest")
    donor = str(row["donor_id"])
    px_manifest = float(row["pixel_size_um_x"])
    builder = entry["serverBuilder"]
    meta = builder.get("metadata") or builder.get("builder", {}).get("metadata") or {}
    px_qupath = float(meta["pixelCalibration"]["pixelWidth"]["value"])
    started = time.time()

    raw = json.loads((export / "islets" / f"{tag}.geojson").read_text())
    features = raw["features"] if isinstance(raw, dict) else raw
    chosen = select_islets(features, rule)
    if not chosen:
        raise ValueError(f"{name}: no object counts as an islet under --islets {rule}")
    bands = nearest_islet_bands({n: o["geom"] for n, o in chosen.items()}, EXPANSION_UM / px_qupath) \
        if expansion == "outlines" else None
    target = export / "typing_islets" / f"Islets__{tag}.geojson"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps({"type": "FeatureCollection", "features": islet_features(chosen, bands)}))

    tsv = senior / "data/RevisedQuPathExports" / f"{donor}.tsv"
    tsv.parent.mkdir(parents=True, exist_ok=True)
    stats = write_tsv(export / "cells" / f"{tag}.parquet", export / "context" / f"{tag}.parquet", tsv,
                      chosen, rule, px_manifest / px_qupath)
    by_class: dict[str, int] = {}
    for o in chosen.values():
        by_class[o["qupath_class"]] = by_class.get(o["qupath_class"], 0) + 1
    return {"image": name, "donor_id": donor, "islets_rule": rule, "expansion": expansion,
            "islets": len(chosen), "islets_by_qupath_class": by_class,
            "expansions": sum(1 for b in (bands or {}).values() if b is not None and not b.is_empty),
            "pixel_um_qupath": px_qupath, "pixel_um_manifest": px_manifest, "typing_islets": str(target),
            "tsv": str(tsv), **stats, "seconds": round(time.time() - started, 1)}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--export", type=Path, required=True)
    ap.add_argument("--senior", type=Path, required=True, help="the Senior checkout whose data/ receives the TSVs")
    ap.add_argument("--islets", choices=sorted(ISLET_CLASSES), required=True,
                    help="which QuPath objects count as islets (the maintainer's ruling; no default)")
    ap.add_argument("--expansion", choices=("outlines", "none"), required=True,
                    help="write the 20 um nearest-islet outlines, or none (the maintainer's ruling; no default)")
    ap.add_argument("--images", nargs="*", default=None, help="image tags; default every image in the Senior manifest")
    ap.add_argument("--summary", type=Path, default=None,
                    help="where to write this run's report (default <export>/typing_islets/adapter_summary.json); several "
                         "runs over disjoint --images can work side by side when each has its own")
    args = ap.parse_args()

    manifest = {r["image_id"]: r for r in json.loads(
        (args.senior / "core/phenocycler/data/qupath_manifest.json").read_text())["images"]}
    project = json.loads(PROJECT.read_text())
    entries = [e for e in project["images"] if e["imageName"] in manifest
               and (not args.images or tag_of(e["imageName"]) in args.images)]
    if args.images and len(entries) != len(args.images):
        raise SystemExit(f"asked for {len(args.images)} images, found {len(entries)} in the manifest")
    reports, failed = [], []
    for entry in entries:
        try:
            report = one_image(entry, args.export, args.senior, manifest, args.islets, args.expansion)
        except Exception as error:                                    # noqa: BLE001 -- one image must not hide the rest
            failed.append((entry["imageName"], f"{type(error).__name__}: {error}"))
            print(f"FAILED {entry['imageName']}: {type(error).__name__}: {error}", flush=True)
            continue
        reports.append(report)
        print(f"{report['donor_id']:>5}: {report['islets']:>5,} islets {report['islets_by_qupath_class']}, "
              f"{report['expansions']:>5,} expansions; {report['rows']:>9,} cells, {report['parent_is_islet']:>7,} in an islet, "
              f"{report['parent_is_pln']:>7,} in a pLN  [{report['seconds']} s]", flush=True)
    summary = args.summary or (args.export / "typing_islets" / "adapter_summary.json")
    summary.parent.mkdir(parents=True, exist_ok=True)
    summary.write_text(json.dumps({"islets_rule": args.islets, "expansion": args.expansion,
                                   "images": reports, "failed": failed}, indent=2))
    print(f"{len(reports)} image(s) written, {len(failed)} failed -> {summary}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
