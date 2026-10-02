#!/usr/bin/env python3
"""pixel_levels.py -- what each stored image holds, per channel, before any of our processing.

For the pre-processing audit of 2026-10-02 ("be sure all 25 were pre-processed in the same way prior to UniFORM
processing"). Everything downstream -- QuPath's cell means, REDSEA, the expression stage, UniFORM -- starts from the
pixels of the qptiff as the scanner software wrote them. This reads those pixels directly, at full resolution, on a
fixed sample of tiles, and writes one row per image, channel and region:

  tissue   48 tiles lying wholly inside the image's Tissue annotation, evenly spaced through it
  glass    up to 24 tiles of scanned, empty slide: farther than 256 px from every annotation, scanned (under 2 % of
           DAPI pixels exactly 0 -- an unscanned area is stored as 0 everywhere), and the lowest DAPI 99th percentile
           among those (no nuclei). The DAPI 99th percentile of the chosen tiles is written beside them so a reader
           can see they are empty.

Per row: the fraction of pixels stored as exactly 0 (the scanner subtracts and stores unsigned, so 0 means "at or
below zero after subtraction"), percentiles, mean and maximum. No cell, no annotation class other than Tissue, and
no measurement of ours enters -- so a difference between images in this table was in the file before we touched it.

It gives levels and no verdict.

    python pixel_levels.py --manifest <worktree>/core/phenocycler/data/qupath_manifest.json \\
        --annotations <project>/results/v3/export/annotations_2026-10-01/manual_annotations \\
        --out <project>/results/v3/typing/pixel_levels.tsv [--jobs 8]
"""
from __future__ import annotations

import argparse
import json
import re
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

TILE = 512                      # the qptiff's own tile size at full resolution: one read decodes one stored tile
N_TISSUE, N_GLASS = 48, 24
GLASS_MARGIN_PX = 256           # a glass tile is at least this far from every annotation
SCANNED_MAX_ZERO = 0.02         # ... and has under this fraction of DAPI pixels exactly 0 (an unscanned area is all 0)
PERCENTILES = (1, 5, 25, 50, 75, 95, 99, 99.9)


def image_tag(image_id: str) -> str:
    """The file stem the v3 exports use for a QuPath image name."""
    return re.sub(r"[^A-Za-z0-9._-]+", "_", image_id).strip("_")


def evenly(items: list, n: int) -> list:
    if len(items) <= n:
        return list(items)
    return [items[i] for i in np.linspace(0, len(items) - 1, n).round().astype(int)]


def from_histogram(hist: np.ndarray) -> dict:
    """Zero fraction, percentiles, mean and maximum of the uint16 values a 65,536-bin histogram counts. Exact."""
    n = int(hist.sum())
    if n == 0:
        return {"n_pixels": 0}
    cum = np.cumsum(hist)
    out = {"n_pixels": n, "frac_zero": float(hist[0]) / n}
    for p in PERCENTILES:
        out[f"p{p:g}"] = int(np.searchsorted(cum, np.ceil(p / 100.0 * n), side="left"))
    values = np.arange(hist.size, dtype=np.float64)
    out["mean"] = float((hist * values).sum() / n)
    out["max"] = int(np.flatnonzero(hist)[-1])
    return out


def one_image(job: tuple) -> list:
    donor, image_id, qptiff, annotation_file, channel_map = job
    import tifffile  # noqa: PLC0415
    import zarr  # noqa: PLC0415
    from shapely.geometry import box, shape  # noqa: PLC0415
    from shapely.ops import unary_union  # noqa: PLC0415
    from shapely.prepared import prep  # noqa: PLC0415

    features = json.loads(Path(annotation_file).read_text())["features"]
    klass = lambda f: ((f.get("properties") or {}).get("classification") or {}).get("name")  # noqa: E731
    tissue = [shape(f["geometry"]) for f in features if klass(f) == "Tissue"]
    if not tissue:
        raise SystemExit(f"{donor}: no Tissue annotation in {annotation_file}")
    tissue = unary_union(tissue)
    everything = unary_union([shape(f["geometry"]) for f in features])
    in_tissue, near_any = prep(tissue), prep(everything.buffer(GLASS_MARGIN_PX))

    marker_of = {int(c["channel_index"]): c["marker"] for c in channel_map}
    rows = []
    with tifffile.TiffFile(qptiff) as tf:
        series = tf.series[0]
        names = []
        for page in series.levels[0].pages:
            m = re.search(r"<Biomarker>(.*?)</Biomarker>", page.description or "")
            names.append(m.group(1) if m else "?")
        z = zarr.open(series.aszarr(level=0), mode="r")
        n_channels, height, width = z.shape
        if names[0] != "DAPI":
            raise SystemExit(f"{donor}: channel 0 is {names[0]!r}, not DAPI")

        tissue_tiles, off_tiles = [], []
        for y0 in range(0, height - TILE + 1, TILE):
            for x0 in range(0, width - TILE + 1, TILE):
                b = box(x0, y0, x0 + TILE, y0 + TILE)
                if in_tissue.contains(b):
                    tissue_tiles.append((y0, x0))
                elif not near_any.intersects(b):
                    off_tiles.append((y0, x0))
        tissue_tiles = evenly(tissue_tiles, N_TISSUE)

        # glass: scanned and empty, judged on DAPI alone
        glass = []
        for (y0, x0) in off_tiles:
            a = np.asarray(z[0, y0:y0 + TILE, x0:x0 + TILE])
            if (a == 0).mean() < SCANNED_MAX_ZERO:
                glass.append((float(np.percentile(a, 99)), y0, x0))
        glass.sort()
        glass = glass[:N_GLASS]
        glass_dapi_p99_max = max((g[0] for g in glass), default=np.nan)
        glass_tiles = [(y0, x0) for _, y0, x0 in glass]

        for c in range(n_channels):
            for region, tiles in (("tissue", tissue_tiles), ("glass", glass_tiles)):
                hist = np.zeros(65536, dtype=np.int64)
                for (y0, x0) in tiles:
                    a = np.asarray(z[c, y0:y0 + TILE, x0:x0 + TILE])
                    hist += np.bincount(a.ravel(), minlength=65536)
                row = {"donor_id": donor, "image": Path(qptiff).name, "channel_index": c, "biomarker": names[c],
                       "marker": marker_of.get(c, ""), "region": region, "n_tiles": len(tiles),
                       "glass_dapi_p99_max": glass_dapi_p99_max if region == "glass" else np.nan}
                row.update(from_histogram(hist))
                rows.append(row)
    return rows


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", type=Path, required=True)
    ap.add_argument("--annotations", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--jobs", type=int, default=8)
    args = ap.parse_args(argv)

    manifest = json.loads(args.manifest.read_text())
    jobs = []
    for image in manifest["images"]:
        annotation = args.annotations / f"{image_tag(image['image_id'])}.geojson"
        if not annotation.is_file():
            raise SystemExit(f"{image['donor_id']}: no annotation export at {annotation}")
        jobs.append((str(image["donor_id"]), image["image_id"], image["qptiff"]["path"], str(annotation), image["channel_map"]))

    rows = []
    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        for part in pool.map(one_image, jobs):
            rows.extend(part)
            print(f"{part[0]['donor_id']}: {len(part) // 2} channels, tissue tiles {part[0]['n_tiles']}, glass tiles {part[1]['n_tiles']}", flush=True)
    table = pd.DataFrame(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.out, sep="\t", index=False, float_format="%.6g")
    print(f"-> {args.out}: {len(table)} rows, {table['donor_id'].nunique()} images")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
