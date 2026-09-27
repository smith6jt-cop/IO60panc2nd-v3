#!/usr/bin/env bash
# make_fixture.sh -- build the run_v3.sh test fixture on /scratch (never touches the real project).
#
#   1. copies project.qpproj to a private file and crops a 2048x2048 full-resolution, all-channel
#      window of 6450panc (40 v2 islets inside the Tissue) to an OME-TIFF;
#   2. checks the OME-TIFF kept the entry's channel names and the pixel size;
#   3. writes the fixture annotations (Tissue rectangle, one pLN inside it, one stale
#      'Other: Gal3Norm' for preflight to remove, and the v2 islets inside the crop);
#   4. creates /scratch/io60/fixture/project/project.qpproj with that image and those annotations.
# The fixture project gets a COPY of classifiers/ and resources/ (a symlink would let a step write
# into the real classifiers/) and a symlink to the real scripts/ (steps only read lib/common.groovy).
#
# Usage: make_fixture.sh [--force]      (--force deletes an existing /scratch/io60/fixture/project)
set -euo pipefail
PROJ=/home/smith6jt/IO60panc2nd
QUPATH=/home/smith6jt/QuPath/bin/QuPath
PY=/scratch/io60/venv-v3/bin/python
FIX=/scratch/io60/fixture
HERE=$(cd "$(dirname "$0")" && pwd)

IMAGE="6450panc_Scan1.er.qptiff - resolution #1"
X=13312 Y=4096 W=2048 H=2048          # 40 v2 islets, window fully inside the Tissue
PIXEL=0.5078                           # the 6450 entry says 0.50776638; mimic the other 23 images
SRC_GEOJSON=$PROJ/results/v3/annotations_all/6450panc_Scan1.er.qptiff_-_resolution_1.geojson

if [[ -e $FIX/project ]]; then
    [[ ${1:-} == --force ]] || { echo "$FIX/project exists; rerun with --force to rebuild" >&2; exit 1; }
    rm -rf -- "$FIX/project"
fi
mkdir -p "$FIX/build"
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
"$PY" - "$SRC_GEOJSON" "$FIX/build/fixture_annotations.geojson" $X $Y $W $H <<'EOF'
import json, sys
from shapely.geometry import shape, box, mapping
from shapely.affinity import translate
src, out = sys.argv[1], sys.argv[2]
x, y, w, h = map(int, sys.argv[3:7])
feats = json.load(open(src))['features']
cls = lambda f: (f['properties'].get('classification') or {}).get('name')
inner = box(x + 24, y + 24, x + w - 24, y + h - 24)
def feat(geom, classification, name=None):
    props = {'objectType': 'annotation', 'classification': classification}
    if name: props['name'] = name
    return {'type': 'Feature', 'geometry': mapping(geom), 'properties': props}
# QuPath 0.7's GeoJSON reader rejects a classification without a color
out_feats = [feat(box(16, 16, w - 16, h - 16), {'name': 'Tissue', 'color': [102, 153, 153]}),
             feat(box(w - 420, h - 420, w - 120, h - 120), {'name': 'pLN', 'color': [0, 128, 255]}, 'fixture_pLN'),
             feat(box(700, 1500, 900, 1640), {'names': ['Other', 'Gal3Norm'], 'color': [255, 0, 0]})]
n = 0
for f in feats:
    if cls(f) != 'Islet': continue
    g = shape(f['geometry'])
    if not inner.contains(g): continue
    out_feats.append(feat(translate(g, -x, -y), {'name': 'Islet', 'color': [255, 0, 0]}, f['properties'].get('name')))
    n += 1
json.dump({'type': 'FeatureCollection', 'features': out_feats}, open(out, 'w'))
print(f'[make_fixture] annotations: Tissue, pLN, Other: Gal3Norm, {n} v2 islets -> {out}')
EOF

JAVA_TOOL_OPTIONS=-Xmx16g "$QUPATH" script -a "$FIX/project" -a "$tif" \
    -a "$FIX/build/fixture_annotations.geojson" "$HERE/fixture_project.groovy" \
    > "$FIX/build/fixture_project.log" 2>&1 || true
[[ -s $FIX/project/project.qpproj ]] || { echo "project creation failed; see $FIX/build/fixture_project.log" >&2; exit 1; }
grep -E '^\[fixture_project\]' "$FIX/build/fixture_project.log"

for d in classifiers resources; do mkdir -p "$FIX/project/$d" && cp -r "$PROJ/$d/." "$FIX/project/$d/"; done
ln -sfn "$PROJ/scripts" "$FIX/project/scripts"
echo "[make_fixture] fixture project ready: $FIX/project/project.qpproj"
