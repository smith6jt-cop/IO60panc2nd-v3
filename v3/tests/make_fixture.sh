#!/usr/bin/env bash
# make_fixture.sh -- build the run_v3.sh test fixture on /scratch (never touches the real project).
#
#   1. copies project.qpproj to a private file and crops a 2048x2048 full-resolution, all-channel
#      window of 6450panc (40 v2 islets inside the Tissue) to an OME-TIFF (once; reused afterwards);
#   2. checks the OME-TIFF kept the entry's channel names and the pixel size;
#   3. writes the fixture annotations for the chosen variant (below), plus one stale
#      'Other: Gal3Norm' for preflight to remove and the v2 islets that lie inside the Tissue;
#   4. creates /scratch/io60/fixture/<project dir>/project.qpproj with that image and those annotations.
# The fixture project gets a COPY of classifiers/ and resources/ (a symlink would let a step write
# into the real classifiers/) and a symlink to the real scripts/ (steps only read lib/common.groovy).
#
# Variants (--variant NAME -> /scratch/io60/fixture/project_NAME; default -> /scratch/io60/fixture/project):
#   (default)    Tissue = the whole crop minus 16 px; one pLN INSIDE it (contained: segmented with Tissue)
#   pln          Tissue = the left 60%; 'pLN_disjoint' on the right, 100 px (51 um) away from Tissue;
#                'pLN_contained' inside Tissue. Tests v3 script 0 segmenting a disjoint lymph node, and
#                step 8 parenting its cells (and the contained one's) to pLN.
#   pln-partial  like pln, but 'pLN_disjoint' reaches 150 px into Tissue (~19% of its area): preflight fails
#   pln-sliver   like pln, but it reaches 5 px into Tissue (<1% of its area, ~1700 um^2): preflight fails
#                with a warning, v3 script 0 refuses (more than 10 um^2 would be segmented twice)
#
# Usage: make_fixture.sh [--force] [--variant pln|pln-partial|pln-sliver]
#   (--force deletes an existing fixture project directory of that variant first)
set -euo pipefail
PROJ=/home/smith6jt/IO60panc2nd
QUPATH=/home/smith6jt/QuPath/bin/QuPath
PY=/scratch/io60/venv-v3/bin/python
FIX=/scratch/io60/fixture
HERE=$(cd "$(dirname "$0")" && pwd)

FORCE=0 VARIANT=""
while (($#)); do
    case $1 in
        --force) FORCE=1 ;;
        --variant) VARIANT=$2; shift ;;
        *) echo "unknown option $1" >&2; exit 2 ;;
    esac
    shift
done
case $VARIANT in
    ''|pln|pln-partial|pln-sliver) ;;
    *) echo "unknown variant '$VARIANT'" >&2; exit 2 ;;
esac
PDIR=$FIX/project${VARIANT:+_${VARIANT//-/_}}

IMAGE="6450panc_Scan1.er.qptiff - resolution #1"
X=13312 Y=4096 W=2048 H=2048          # 40 v2 islets, window fully inside the Tissue
PIXEL=0.5078                           # the 6450 entry says 0.50776638; mimic the other 23 images
SRC_GEOJSON=$PROJ/results/v3/annotations_all/6450panc_Scan1.er.qptiff_-_resolution_1.geojson

if [[ -e $PDIR ]]; then
    ((FORCE)) || { echo "$PDIR exists; rerun with --force to rebuild" >&2; exit 1; }
    rm -rf -- "$PDIR"
fi
BUILD=$FIX/build${VARIANT:+/$VARIANT}
mkdir -p "$BUILD"
cp -f "$PROJ/project.qpproj" "$FIX/build/master_copy.qpproj"   # private copy: QuPath never opens the master

tif=$FIX/fixture_6450.ome.tif
if [[ ! -s $tif ]]; then
    JAVA_TOOL_OPTIONS=-Xmx16g "$QUPATH" script -a "$FIX/build/master_copy.qpproj" -a "$IMAGE" \
        -a $X -a $Y -a $W -a $H -a "${tif%.ome.tif}.tmp.ome.tif" -a $PIXEL "$HERE/fixture_crop.groovy" \
        > "$FIX/build/fixture_crop.log" 2>&1 || true
    # Bio-Formats picks the writer from the extension, so the temporary name keeps .ome.tif
    [[ -s ${tif%.ome.tif}.tmp.ome.tif ]] || { echo "crop failed; see $FIX/build/fixture_crop.log" >&2; exit 1; }
    mv "${tif%.ome.tif}.tmp.ome.tif" "$tif"
fi

# channel names + pixel size must survive the OME-TIFF round trip
"$PY" - "$tif" "$FIX/build/master_copy.qpproj" "$IMAGE" "$PIXEL" <<'EOF'
import json, sys, xml.etree.ElementTree as ET, tifffile
tif, proj, image, pixel = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4])
entry = next(e for e in json.load(open(proj))['images'] if e['imageName'] == image)
want = [c['name'] for c in entry['serverBuilder']['metadata']['channels']]
with tifffile.TiffFile(tif) as t:
    root = ET.fromstring(t.ome_metadata)
    ns = {'ome': root.tag.split('}')[0].strip('{')}
    px = root.find('.//ome:Pixels', ns)
    got = [c.get('Name') for c in px.findall('ome:Channel', ns)]
    levels = len(t.series[0].levels)
    shape, dtype = t.series[0].shape, t.series[0].dtype
psx, psy = float(px.get('PhysicalSizeX')), float(px.get('PhysicalSizeY'))
assert got == want, f'channel names differ: {got[:5]}... vs {want[:5]}...'
assert abs(psx - pixel) < 1e-9 and abs(psy - pixel) < 1e-9, (psx, psy)
print(f'[make_fixture] OME-TIFF ok: {len(got)} channels match the entry, pixel {psx} x {psy} '
      f'{px.get("PhysicalSizeXUnit")}, shape {shape} {dtype}, {levels} pyramid levels')
EOF

# fixture annotations in crop coordinates
"$PY" - "$SRC_GEOJSON" "$BUILD/fixture_annotations.geojson" $X $Y $W $H "$VARIANT" <<'EOF'
import json, sys
from shapely.geometry import shape, box, mapping
from shapely.affinity import translate
src, out = sys.argv[1], sys.argv[2]
x, y, w, h = map(int, sys.argv[3:7])
variant = sys.argv[7] if len(sys.argv) > 7 else ''
feats = json.load(open(src))['features']
cls = lambda f: (f['properties'].get('classification') or {}).get('name')
def feat(geom, classification, name=None):
    props = {'objectType': 'annotation', 'classification': classification}
    if name: props['name'] = name
    return {'type': 'Feature', 'geometry': mapping(geom), 'properties': props}
PLN = {'name': 'pLN', 'color': [0, 128, 255]}          # QuPath 0.7's GeoJSON reader rejects a class without a color
if not variant:
    tissue = box(16, 16, w - 16, h - 16)
    plns = [(box(w - 420, h - 420, w - 120, h - 120), 'fixture_pLN')]
else:
    tissue = box(16, 16, 1250, h - 16)                 # left ~60% of the crop
    left = {'pln': 1350, 'pln-partial': 1100, 'pln-sliver': 1245}[variant]
    plns = [(box(left, 200, w - 16, 1500), 'pLN_disjoint'),     # 'pln': 100 px (51 um) right of Tissue
            (box(900, 1700, 1150, 1950), 'pLN_contained')]      # inside Tissue
out_feats = [feat(tissue, {'name': 'Tissue', 'color': [102, 153, 153]})]
out_feats += [feat(g, PLN, n) for g, n in plns]
out_feats.append(feat(box(700, 1500, 900, 1640), {'names': ['Other', 'Gal3Norm'], 'color': [255, 0, 0]}))
# default: exactly the original fixture (islets 24 px inside the crop); variants: inside Tissue, off every pLN
inner = box(24, 24, w - 24, h - 24) if not variant else tissue.buffer(-24)
avoid = [] if not variant else [g for g, _ in plns]
n = 0
for f in feats:
    if cls(f) != 'Islet': continue
    g = translate(shape(f['geometry']), -x, -y)
    if not inner.contains(g) or any(g.intersects(p) for p in avoid): continue
    out_feats.append(feat(g, {'name': 'Islet', 'color': [255, 0, 0]}, f['properties'].get('name')))
    n += 1
json.dump({'type': 'FeatureCollection', 'features': out_feats}, open(out, 'w'))
for g, name in plns:
    inside = g.intersection(tissue).area
    print(f'[make_fixture] {name}: {g.area * 0.5078 ** 2:.0f} um^2, {inside * 0.5078 ** 2:.0f} um^2 = '
          f'{inside / g.area:.4f} of it inside Tissue')
print(f'[make_fixture] annotations: Tissue, {len(plns)} pLN, Other: Gal3Norm, {n} v2 islets -> {out}')
EOF

JAVA_TOOL_OPTIONS=-Xmx16g "$QUPATH" script -a "$PDIR" -a "$tif" \
    -a "$BUILD/fixture_annotations.geojson" "$HERE/fixture_project.groovy" \
    > "$BUILD/fixture_project.log" 2>&1 || true
[[ -s $PDIR/project.qpproj ]] || { echo "project creation failed; see $BUILD/fixture_project.log" >&2; exit 1; }
grep -E '^\[fixture_project\]' "$BUILD/fixture_project.log"

for d in classifiers resources; do mkdir -p "$PDIR/$d" && cp -r "$PROJ/$d/." "$PDIR/$d/"; done
ln -sfn "$PROJ/scripts" "$PDIR/scripts"
echo "[make_fixture] fixture project ready: $PDIR/project.qpproj"
