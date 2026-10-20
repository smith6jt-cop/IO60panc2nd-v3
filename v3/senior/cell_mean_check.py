#!/usr/bin/env python3
"""cell_mean_check.py -- is each image's measurement of a marker that image's channel of that marker?

For the pre-processing audit of 2026-10-02. The four 70-channel images store their channels in another order than
the 58/59-channel images, so the one fault that would make a marker read differently on them through OUR processing
is a measurement taken from the wrong channel. This checks it directly, on every image, with nothing in between:

  for a sample of cells, the mean of the stored full-resolution pixels of the page whose own <Biomarker> is M,
  inside the cell's exported outline,
  against QuPath's "Cell: M: Mean" for the same cell (the value every later stage starts from).

Per image and marker it writes the number of cells, the median ratio (QuPath / pixels), the share of cells within
2 %, the correlation -- and the same correlation against the NEXT page of the file, which is what a channel shifted
by one would give. It reads; it changes nothing.

By default six markers are checked; --all-channels checks every page of every image (the next-page column is then
left out: it doubles the reading and the six-marker run already carries it).

    python cell_mean_check.py --manifest <worktree>/core/phenocycler/data/qupath_manifest.json \\
        --out <project>/results/v3/typing/preprocessing_cell_mean_check.tsv [--cells 240] [--jobs 8] [--all-channels]
"""
from __future__ import annotations

import argparse
import json
import re
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

MARKERS = ("DAPI", "CD8", "CD163", "Podoplanin", "CD3e", "INS")
N_OFFSETS = 8                    # places in the outline file the sample is drawn from
CHUNK = 6 * 1024 * 1024


def features_at(path: Path, offset: int, want: int) -> list:
    """Whole cell features starting after `offset` in a QuPath GeoJSON export (the file is too large to parse)."""
    out = []
    with open(path, "rb") as fh:
        fh.seek(offset)
        text = fh.read(CHUNK).decode("utf-8", errors="ignore")
    for m in re.finditer(r'\{"type":\s*"Feature"', text):
        depth, i, in_string, escape = 0, m.start(), False, False
        while i < len(text):
            ch = text[i]
            if in_string:
                if escape:
                    escape = False
                elif ch == "\\":
                    escape = True
                elif ch == '"':
                    in_string = False
            elif ch == '"':
                in_string = True
            elif ch == "{":
                depth += 1
            elif ch == "}":
                depth -= 1
                if depth == 0:
                    break
            i += 1
        if depth != 0:
            break                                       # the chunk ended inside this feature
        try:
            out.append(json.loads(text[m.start():i + 1]))
        except json.JSONDecodeError:
            continue
        if len(out) >= want:
            break
    return out


def one_image(job: tuple) -> list:
    donor, qptiff, geojson, cells_table, n_cells, all_channels = job
    import pyarrow.dataset as ds  # noqa: PLC0415
    import shapely  # noqa: PLC0415
    import tifffile  # noqa: PLC0415
    import zarr  # noqa: PLC0415
    from shapely.geometry import shape  # noqa: PLC0415

    size = Path(geojson).stat().st_size
    per = max(1, n_cells // N_OFFSETS)
    cells = {}
    for k in range(N_OFFSETS):
        for f in features_at(Path(geojson), int(size * (k + 0.5) / N_OFFSETS), per):
            if f.get("id") and f.get("geometry"):
                cells[f["id"]] = shape(f["geometry"])
    ids = list(cells)

    with tifffile.TiffFile(qptiff) as tf:
        series = tf.series[0]
        names = []
        for page in series.levels[0].pages:
            m = re.search(r"<Biomarker>(.*?)</Biomarker>", page.description or "")
            names.append(m.group(1) if m else "?")
        z = zarr.open(series.aszarr(level=0), mode="r")
        n_channels, height, width = z.shape
        present = list(names) if all_channels else [m for m in MARKERS if m in names]
        if len(set(names)) != len(names):
            raise SystemExit(f"{donor}: two pages carry the same <Biomarker> ({sorted(n for n in set(names) if names.count(n) > 1)})")
        index = {m: names.index(m) for m in present}
        ours = {m: np.full(len(ids), np.nan) for m in present}
        nxt = {m: np.full(len(ids), np.nan) for m in present}
        for j, oid in enumerate(ids):
            geom = cells[oid]
            x0, y0, x1, y1 = geom.bounds
            c0, r0, c1, r1 = max(int(np.floor(x0)), 0), max(int(np.floor(y0)), 0), min(int(np.ceil(x1)), width), min(int(np.ceil(y1)), height)
            if c1 <= c0 or r1 <= r0:
                continue
            xs, ys = np.meshgrid(np.arange(c0, c1) + 0.5, np.arange(r0, r1) + 0.5)
            inside = shapely.contains_xy(geom, xs, ys)
            if not inside.any():
                continue
            for m in present:
                c = index[m]
                ours[m][j] = np.asarray(z[c, r0:r1, c0:c1])[inside].mean()
                if not all_channels:
                    nxt[m][j] = np.asarray(z[(c + 1) % n_channels, r0:r1, c0:c1])[inside].mean()

    dataset = ds.dataset(cells_table)
    measured = set(dataset.schema.names)
    columns = ["Object ID"] + [f"Cell: {m}: Mean" for m in present if f"Cell: {m}: Mean" in measured]
    table = dataset.to_table(columns=columns, filter=ds.field("Object ID").isin(ids)).to_pandas().set_index("Object ID")
    rows = []
    for m in present:
        if f"Cell: {m}: Mean" not in measured:                 # a stored channel with no measurement of that name
            rows.append({"donor_id": donor, "marker": m, "page": index[m], "n_cells": 0})
            continue
        q = table[f"Cell: {m}: Mean"].reindex(ids).to_numpy(dtype=float)
        o, n = ours[m], nxt[m]
        ok = np.isfinite(q) & np.isfinite(o)
        ratio = q[ok] / np.where(o[ok] > 0, o[ok], np.nan)
        rows.append({"donor_id": donor, "marker": m, "page": index[m], "next_page_biomarker": names[(index[m] + 1) % len(names)],
                     "n_cells": int(ok.sum()), "median_ratio_qupath_over_pixels": float(np.nanmedian(ratio)),
                     "share_within_2pct": float(np.nanmean(np.abs(ratio[np.isfinite(ratio)] - 1) <= 0.02)),
                     "median_abs_difference": float(np.median(np.abs(q[ok] - o[ok]))),
                     "r_named_page": float(np.corrcoef(q[ok], o[ok])[0, 1]) if ok.sum() > 2 else np.nan,
                     "r_next_page": float(np.corrcoef(q[ok], n[ok])[0, 1]) if ok.sum() > 2 and not all_channels else np.nan})
    return rows


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--cells", type=int, default=240)
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--all-channels", action="store_true")
    args = ap.parse_args(argv)
    manifest = json.loads(args.manifest.read_text())
    jobs = [(str(i["donor_id"]), i["qptiff"]["path"], i["cell_geojson"]["path"], i["measurement_csv"]["path"], args.cells, args.all_channels)
            for i in manifest["images"]]
    rows = []
    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        for part in pool.map(one_image, jobs):
            rows.extend(part)
            print(f"{part[0]['donor_id']}: {part[0]['n_cells']} cells", flush=True)
    table = pd.DataFrame(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.out, sep="\t", index=False, float_format="%.5g")
    print(f"-> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
