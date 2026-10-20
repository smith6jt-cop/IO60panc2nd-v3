#!/usr/bin/env python3
"""make_cohort_spec.py -- the Senior cohort specification for the v3 exports (25 donors).

Starts from the TRACKED spec of the Senior checkout it is pointed at (`core/phenocycler/cohort-spec.json`:
donor ids and panel ids for 23 images) and changes only what v3 changes:

  - each image's measurement source is its own cells Parquet from the v3 export root, and its cell
    outlines the v3 GeoJSON (one combined file per image: cell + nucleus geometry);
  - the two same-slide donors the 23 do not list are added, with their donor ids DECLARED here
    (never parsed from a file name): 6347 and 6580, both on the 70-channel panel of 6533 and 0112;
  - with --context, each image also names its context Parquet (pass J), making it a (cells, context) pair.

6457 is never listed: it is centrally excluded (`phenocycler.cohort.DONOR_EXCLUSIONS`) and is not processed.
An image whose export files are not on disk yet is reported and the spec is not written unless --allow-missing.

    make_cohort_spec.py --senior <checkout> --export <results/v3/export> --out <spec.json> [--context]
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

PROJECT = Path("/home/smith6jt/IO60panc2nd/project.qpproj")
#: declared, hand-recorded: image -> (donor, panel). The other 23 come from the tracked spec.
ADDED = {
    "0080579_6347_pLNpanc_Scan1.er.qptiff - resolution #1": ("6347", "panel-metabolic-70ch"),
    "0080502_6580_pLNpanc_Scan1.er.qptiff - resolution #1": ("6580", "panel-metabolic-70ch"),
}
SEGMENTATION = "qupath-0.7.0__instanseg-fluorescence_nuclei_and_cells-0.1.1__v3-13ch"


def tag_of(image: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", image)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--senior", type=Path, required=True, help="the Senior checkout (worktree) to read the tracked spec from")
    ap.add_argument("--export", type=Path, required=True, help="the v3 export root (cells/, cell_geojson/, context/)")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--context", action="store_true", help="pair every image with context/<tag>.parquet")
    ap.add_argument("--allow-missing", action="store_true")
    args = ap.parse_args()

    sys.path.insert(0, str(args.senior / "core/phenocycler"))
    from phenocycler.cohort import DONOR_EXCLUSIONS  # noqa: PLC0415

    tracked = json.loads((args.senior / "core/phenocycler/cohort-spec.json").read_text())
    known = {row["image_id"]: (str(row["donor_id"]), row["panel_id"]) for row in tracked["images"]}
    overlap = set(known) & set(ADDED)
    if overlap:
        raise SystemExit(f"already in the tracked spec, remove from ADDED: {sorted(overlap)}")
    known.update(ADDED)

    project = json.loads(PROJECT.read_text())
    images, missing, skipped = [], [], []
    for entry in project["images"]:
        name = entry["imageName"]
        if name not in known:
            skipped.append(name)
            continue
        donor, panel = known[name]
        if donor in DONOR_EXCLUSIONS:
            raise SystemExit(f"{name}: donor {donor} is centrally excluded ({DONOR_EXCLUSIONS[donor]})")
        builder = entry["serverBuilder"]
        qptiff = (builder.get("uri") or builder["builder"]["uri"]).replace("file://", "").replace("file:", "")
        tag = tag_of(name)
        row = {"donor_id": donor, "image_id": name, "qptiff": qptiff, "panel_id": panel,
               "measurement_csv": str(args.export / "cells" / f"{tag}.parquet"),
               "cell_geojson": str(args.export / "cell_geojson" / f"{tag}.geojson")}
        need = [row["measurement_csv"], row["cell_geojson"], qptiff]
        if args.context:
            row["measurement_context"] = str(args.export / "context" / f"{tag}.parquet")
            need.append(row["measurement_context"])
        gone = [p for p in need if not Path(p).is_file()]
        if gone:
            missing.append((donor, gone))
        images.append(row)
    not_in_project = sorted(set(known) - {e["imageName"] for e in project["images"]})
    if not_in_project:
        raise SystemExit(f"in the spec but not in the QuPath project: {not_in_project}")
    donors = [row["donor_id"] for row in images]
    if len(set(donors)) != len(donors):
        raise SystemExit("a donor has two images")

    spec = {"spec_version": tracked["spec_version"], "cohort_name": "islet-phenocycler-v3",
            "expected_donors": sorted(donors),
            "defaults": {"segmentation_version": SEGMENTATION, "nucleus_geojson": None,
                         "pixel_size_um_x": tracked["defaults"]["pixel_size_um_x"],
                         "pixel_size_um_y": tracked["defaults"]["pixel_size_um_y"], "channel_map": None},
            "images": sorted(images, key=lambda r: r["donor_id"])}
    print(f"{len(images)} images, {len(set(donors))} donors; not listed (by design): {skipped}")
    for donor, gone in missing:
        print(f"  not on disk yet for {donor}: {[Path(p).name for p in gone]}")
    if missing and not args.allow_missing:
        print("spec NOT written (files missing; --allow-missing to write it anyway)")
        return 1
    args.out.write_text(json.dumps(spec, indent=2) + "\n")
    print("wrote", args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
