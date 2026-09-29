#!/usr/bin/env python3
"""islet_benchmark_report.py -- coverage and leak check for the v3 islet benchmark (islet_benchmark_tiles.groovy).

Coverage: tiles per CATEGORY x STRATUM (acquisition group x disease) against the quota, from
OUT_DIR/private/key.tsv and the per-image pool records in OUT_DIR/private/runs/*.json. Every cell that is not
filled is flagged:
  NOT_RUN     no image of that stratum has been tiled yet
  PARTIAL     below the quota while some of the stratum's images are still untiled
  POOL_EMPTY  every image of the stratum is tiled and every pool was empty (e.g. no insulin-deficient islets
              in ND): a real absence
  SHORT       every image is tiled, candidates existed, but fewer tiles than the quota could be placed
Written to OUT_DIR/private/coverage.tsv (long form) and printed as a matrix.

Leak check (every file under OUT_DIR/annotate/, which is what annotators see):
  - only tiles/<code>.ome.tif files; codes unique; nothing else (no key, runs, preds);
  - every TIFF tag string of every page (OME-XML ImageDescription, Software, ...) is scanned for project image
    names, tags, stems, donor IDs (whole tokens), 'qptiff', 'Scan<N>', 'panc', 'pLN', 'HDL' and pixel sizes other
    than the fixed 0.5078; the raw bytes are scanned for the long identifiers (names, tags, stems);
  - OME-XML: Image Name == file code, PhysicalSizeX/Y == 0.5078, channels == DAPI/CD99/INS/GCG/SST with the fixed
    palette, one tile size for all tiles, no acquisition/instrument/description elements;
  - XResolution/YResolution consistent with 0.5078 um; every file mtime equal (the fixed date).
Written to OUT_DIR/private/leak_check.tsv; exit status 1 on any failure.

Usage:
  islet_benchmark_report.py --out-dir /scratch/io60/benchmark_v2 [--project /home/smith6jt/IO60panc2nd/project.qpproj]
"""
from __future__ import annotations

import argparse
import json
import math
import re
import sys
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path

# user ruling 2026-09-28 (same table as islet_benchmark_tiles.groovy): donor -> (group, disease[, excluded])
STRATA = {
    "6476": ("E1", "ND"), "HDL115": ("E1", "ND"), "6521": ("E1", "Aab+"), "6579": ("E1", "T1D"), "6380": ("E1", "T1D"),
    "6551": ("E1", "T1D"), "6482": ("E2", "ND"), "6374": ("E2", "ND"), "6493": ("E2", "ND"), "6539": ("E2", "ND"),
    "6442": ("E2", "ND"), "6461": ("E2", "ND"), "6593": ("E2", "T1D"), "6436": ("E2", "T1D"), "6457": ("E2", "T1D", "excluded"),
    "6505": ("E3", "Aab+"), "6450": ("E3", "Aab+"), "6623": ("E3", "T1D"), "6534": ("E3", "T1D"), "6566": ("E3", "T1D"),
    "6591": ("E3", "T1D"), "6523": ("E3", "T1D"), "HDL112": ("E4", "ND"), "6347": ("E4", "Aab+"), "6533": ("E4", "T1D"),
    "6580": ("E4", "T1D"),
}
CATEGORIES = [
    "pos_islet_small", "pos_islet_medium", "pos_islet_large", "pos_ins_rich", "pos_ins_deficient", "pos_infiltrated",
    "pos_cluster", "pos_single", "pos_cd99_only", "pos_hormone_only", "pos_legacy_only", "hard_near_threshold",
    "neg_exocrine", "neg_exo_cd99", "neg_speckle", "neg_duct", "neg_vessel", "neg_immune", "neg_edge", "neg_fold",
]
OUT_CHANNELS = ["DAPI", "CD99", "INS", "GCG", "SST"]
PIXEL_SIZE = "0.5078"
OME_NS = "{http://www.openmicroscopy.org/Schemas/OME/2016-06}"
# OME elements a tile may contain; anything else (AcquisitionDate, Instrument, Description, ...) is flagged
OME_ALLOWED = {"OME", "Image", "Pixels", "Channel", "LightPath", "TiffData", "UUID", "StructuredAnnotations",
               "MapAnnotation", "Value", "M"}


def strata_table(path: Path | None) -> dict:
    if path is None:
        return dict(STRATA)
    out = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith(("#", "donor\t")):
            continue
        f = line.split("\t")
        out[f[0].strip()] = (f[1].strip(), f[2].strip()) + (("excluded",) if len(f) > 3 and f[3].strip() in ("1", "true", "excluded") else ())
    return out


def stratum_sizes(table: dict) -> Counter:
    return Counter(f"{v[0]}_{v[1]}" for v in table.values() if len(v) < 3)


def read_key(path: Path) -> list[dict]:
    if not path.exists():
        return []
    lines = [l for l in path.read_text(encoding="utf-8").splitlines() if l.strip()]
    head = lines[0].split("\t")
    return [dict(zip(head, l.split("\t"))) for l in lines[1:]]


def coverage(out_dir: Path, table: dict) -> list[dict]:
    priv = out_dir / "private"
    key = read_key(priv / "key.tsv")
    runs = [json.loads(p.read_text(encoding="utf-8")) for p in sorted((priv / "runs").glob("*.json"))]
    quota = max([int(r["params"]["QUOTA_PER_STRATUM"]) for r in runs] or [3])
    got = Counter((r["tile_kind"], r["stratum"]) for r in key)
    pooled = defaultdict(list)     # (cat, stratum) -> pool sizes of the tiled images
    tiled = defaultdict(set)
    for r in runs:
        tiled[r["stratum"]].add(r["image"])
        for cat, st in (r.get("pools") or {}).items():
            pooled[(cat, r["stratum"])].append(int(st["pool_size"]))
    sizes = stratum_sizes(table)
    cats = CATEGORIES + sorted({c for c, _ in pooled} - set(CATEGORIES))
    rows = []
    for s in sorted(sizes):
        for c in cats:
            n = got.get((c, s), 0)
            pools = pooled.get((c, s), [])
            if not tiled.get(s):
                status = "NOT_RUN"
            elif n >= quota:
                status = "OK"
            elif len(tiled[s]) < sizes[s]:
                status = "PARTIAL"          # more images of this stratum still to tile
            elif pools and max(pools) == 0:
                status = "POOL_EMPTY"
            else:
                status = "SHORT"
            rows.append(dict(category=c, stratum=s, tiles=n, quota=quota, images_tiled=len(tiled.get(s, ())),
                             images_in_stratum=sizes[s], pool_sizes=",".join(map(str, pools)), status=status))
    return rows


def print_matrix(rows: list[dict]) -> None:
    strata = sorted({r["stratum"] for r in rows})
    cats = list(dict.fromkeys(r["category"] for r in rows))
    cell = {(r["category"], r["stratum"]): r for r in rows}
    code = {"OK": "", "NOT_RUN": "-", "POOL_EMPTY": "!0", "SHORT": "!", "PARTIAL": "+"}
    w = max(len(c) for c in cats)
    print(" " * w + " | " + " | ".join(f"{s:>9s}" for s in strata))
    for c in cats:
        vals = []
        for s in strata:
            r = cell[(c, s)]
            vals.append(f"{'-' if r['status'] == 'NOT_RUN' else str(r['tiles']) + '/' + str(r['quota'])}{code[r['status']] if r['status'] != 'NOT_RUN' else ''}".rjust(9))
        print(f"{c:<{w}s} | " + " | ".join(vals))
    n = Counter(r["status"] for r in rows)
    print(f"cells: {dict(n)}  ('+' PARTIAL: stratum images still untiled, '!' SHORT, '!0' POOL_EMPTY, '-' NOT_RUN)")


# ------------------------------------------------------------------------------------------ leak check

def identifiers(project: Path | None, key_rows: list[dict], table: dict) -> tuple[list[str], list[str]]:
    """(long identifiers scanned in the raw bytes, token patterns scanned in the metadata strings)."""
    names = {r["image"] for r in key_rows}
    if project is not None and project.exists():
        proj = json.loads(project.read_text(encoding="utf-8"))
        names |= {e["imageName"] for e in proj["images"]}
    longs = set()
    for n in names:
        tag = re.sub(r"[^A-Za-z0-9._-]+", "_", n)
        stem = re.sub(r"(?i)\.(er\.)?qptiff.*$", "", n)
        stem = re.sub(r"(?i)\.ome\.tif+.*$", "", stem)
        longs |= {x for x in (n, tag, stem) if len(x) >= 8}
    donors = set(table) | {r.get("donor", "") for r in key_rows if r.get("donor")}
    tokens = [rf"(?<![0-9A-Za-z]){re.escape(d)}(?![0-9])" for d in sorted(donors) if d]
    tokens += [r"(?i)qptiff", r"(?i)scan[0-9]", r"(?i)panc", r"pLN", r"(?i)HDL", r"0\.50776", r"0\.5077[0-9]"]
    return sorted(longs), tokens


def tag_strings(tf) -> list[tuple[str, str]]:
    out = []
    for i, page in enumerate(tf.pages):
        for t in page.tags.values():
            v = t.value
            if isinstance(v, bytes):
                v = v.decode("latin-1", "replace")
            if isinstance(v, str):
                out.append((f"page{i}:{t.name}", v))
    return out


def check_tile(path: Path, longs: list[str], tokens: list[str], palette: dict | None) -> tuple[list[str], dict]:
    import tifffile
    probs = []
    info = {}
    raw = path.read_bytes()
    for s in longs:
        if s.encode("utf-8") in raw or s.encode("utf-16-le") in raw:
            probs.append(f"raw bytes contain '{s}'")
    with tifffile.TiffFile(path) as tf:
        strings = tag_strings(tf)
        p0 = tf.pages[0]
        info["shape"] = tuple(tf.series[0].shape)
        xr = p0.tags.get("XResolution")
        ru = p0.tags.get("ResolutionUnit")
        if xr is not None and ru is not None:
            num, den = xr.value
            per_unit = num / den
            unit_um = {2: 25400.0, 3: 10000.0}.get(int(ru.value), None)
            if unit_um:
                px = unit_um / per_unit
                info["xres_um"] = px
                if abs(px - float(PIXEL_SIZE)) > 5e-5:
                    probs.append(f"XResolution gives {px:.7f} um/px (expected {PIXEL_SIZE})")
    for where, s in strings:
        s = re.sub(r"(?i)urn:uuid:[0-9a-f-]+", "urn:uuid:<random>", s)   # random UUIDs may contain digit runs
        for pat in tokens:
            m = re.search(pat, s)
            if m:
                probs.append(f"{where} contains '{m.group(0)}'")
    desc = next((s for w, s in strings if w == "page0:ImageDescription"), None)
    if desc is None or "<OME" not in desc:
        probs.append("no OME-XML ImageDescription")
        return probs, info
    root = ET.fromstring(desc)
    for el in root.iter():
        tagname = el.tag.replace(OME_NS, "").split("}")[-1]
        if tagname not in OME_ALLOWED:
            probs.append(f"OME element <{tagname}> present")
    img = root.find(f"{OME_NS}Image")
    pix = img.find(f"{OME_NS}Pixels") if img is not None else None
    code = path.name.split(".")[0]
    if img is None or img.get("Name") != code:
        probs.append(f"OME Image Name {None if img is None else img.get('Name')!r} != code {code}")
    if pix is not None:
        for a in ("PhysicalSizeX", "PhysicalSizeY"):
            if pix.get(a) != PIXEL_SIZE:
                probs.append(f"{a}={pix.get(a)} (expected {PIXEL_SIZE})")
        chans = pix.findall(f"{OME_NS}Channel")
        names = [c.get("Name") for c in chans]
        colors = [c.get("Color") for c in chans]
        info["channels"] = names
        info["colors"] = colors
        if names != OUT_CHANNELS:
            probs.append(f"channels {names} != {OUT_CHANNELS}")
        if palette is not None and colors != palette.get("colors"):
            probs.append(f"channel colors {colors} differ from the first tile's {palette.get('colors')}")
        info["size"] = (pix.get("SizeX"), pix.get("SizeY"))
    return probs, info


def leak_check(out_dir: Path, project: Path | None, table: dict) -> list[dict]:
    ann = out_dir / "annotate"
    key_rows = read_key(out_dir / "private" / "key.tsv")
    longs, tokens = identifiers(project, key_rows, table)
    rows = []
    files = sorted(p for p in ann.rglob("*") if p.is_file())
    tiles = [p for p in files if p.parent == ann / "tiles" and re.fullmatch(r"T[A-Z0-9]{8}\.ome\.tif", p.name)]
    for p in sorted(set(files) - set(tiles)):
        rows.append(dict(file=str(p.relative_to(out_dir)), ok=False, problems="not a tile: annotate/ must hold tiles only"))
    palette, sizes, mtimes = None, Counter(), Counter()
    codes = Counter(p.name.split(".")[0] for p in tiles)
    for p in tiles:
        probs, info = check_tile(p, longs, tokens, palette)
        if palette is None and info.get("colors"):
            palette = {"colors": info["colors"]}
        sizes[info.get("size")] += 1
        mtimes[int(p.stat().st_mtime)] += 1
        if codes[p.name.split(".")[0]] > 1:
            probs.append("duplicate code")
        rows.append(dict(file=str(p.relative_to(out_dir)), ok=not probs, problems="; ".join(probs)))
    if len(sizes) > 1:
        rows.append(dict(file="(all tiles)", ok=False, problems=f"tile sizes differ: {dict(sizes)}"))
    if len(mtimes) > 1:
        rows.append(dict(file="(all tiles)", ok=False, problems=f"file mtimes differ ({len(mtimes)} values): write order leaks"))
    keyed = {r["code"] for r in key_rows}
    missing = sorted(keyed - set(codes))
    orphans = sorted(set(codes) - keyed)
    if missing or orphans:
        rows.append(dict(file="(key)", ok=False, problems=f"key rows without tile: {missing[:5]}... ; tiles without key row: {orphans[:5]}..."))
    return rows


def write_tsv(rows: list[dict], path: Path) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    cols = list(rows[0])
    path.write_text("\t".join(cols) + "\n" + "".join("\t".join(str(r[c]) for c in cols) + "\n" for r in rows), encoding="utf-8")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--project", type=Path, default=Path("/home/smith6jt/IO60panc2nd/project.qpproj"),
                    help="read-only: the image names to scan for")
    ap.add_argument("--strata-tsv", type=Path)
    ap.add_argument("--deliver", type=Path,
                    help="after a clean leak check, write the tiles to this .zip in random order with one fixed "
                         "timestamp: hand annotators the zip, never the annotate/ folder, whose creation times and "
                         "directory order group the tiles by image (review 2026-09-29)")
    a = ap.parse_args(argv)
    table = strata_table(a.strata_tsv)
    cov = coverage(a.out_dir, table)
    write_tsv(cov, a.out_dir / "private" / "coverage.tsv")
    print_matrix(cov)
    leaks = leak_check(a.out_dir, a.project, table)
    write_tsv(leaks, a.out_dir / "private" / "leak_check.tsv")
    bad = [r for r in leaks if not r["ok"]]
    n_tiles = sum(1 for r in leaks if r["file"].startswith("annotate/tiles/"))
    print(f"leak check: {n_tiles} tile(s), {len(bad)} problem row(s)" + ("" if not bad else ":"))
    for r in bad[:20]:
        print(f"  {r['file']}: {r['problems']}")
    if a.deliver and not bad:
        n = deliver(a.out_dir, a.deliver)
        print(f"delivered {n} tile(s) in random order -> {a.deliver}")
    return 1 if bad else 0


def deliver(out_dir: Path, zip_path: Path) -> int:
    """Pack annotate/tiles/*.ome.tif into a zip in a random order, every entry dated 1980-01-01, so neither the
    member order, nor timestamps, nor the ctimes/inode order of an unpacked copy say which tiles share an image."""
    import random
    import zipfile
    tiles = sorted((out_dir / "annotate" / "tiles").glob("*.ome.tif"))
    random.SystemRandom().shuffle(tiles)
    zip_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = zip_path.with_name(zip_path.name + ".part")
    with zipfile.ZipFile(tmp, "w", compression=zipfile.ZIP_STORED) as z:
        for t in tiles:
            info = zipfile.ZipInfo(f"tiles/{t.name}", date_time=(1980, 1, 1, 0, 0, 0))
            info.external_attr = 0o644 << 16
            with open(t, "rb") as f:
                z.writestr(info, f.read())
    tmp.replace(zip_path)
    return len(tiles)


if __name__ == "__main__":
    sys.exit(main())
