/**
 * export_benchmark_truth.groovy -- write the annotator's outlines on one blinded benchmark tile to
 * <OUT_DIR>/truth/<code>.geojson, the ground truth islet_benchmark_score.py reads (truth_dir/<code>.geojson).
 *
 * Run it on the ANNOTATION project (the 534 tiles of benchmark_tiles_for_annotation.zip), never on the
 * master project: Automate > Script editor > Run > Run for project. Every tile gets a file; a tile with
 * no outline is written as an empty FeatureCollection, which the scorer reads as a true negative.
 *
 * Exported: every area annotation, in the tile's own pixel coordinates -- the frame of the detector's
 * private/pred/<label>/<code>.geojson (0..788 px at 0.5078 um). The scorer recomputes the size class
 * (EndocrineSingle / EndocrineCluster / Islet) from geometry, so the class set while drawing does not
 * matter, except "Ambiguous": it is kept in the file so those objects can be excluded at scoring.
 */

// -a OUT_DIR=<dir> overrides the destination (tests); the GUI's Run for project uses the default.
def ARGS = binding.hasVariable('args') ? (args as List).collect { it.toString() } : []
String OUT_DIR = ARGS.find { it.startsWith('OUT_DIR=') }?.substring(8) ?: '/home/smith6jt/IO60panc2nd/results/v3/benchmark_annotation'

def entry = getProjectEntry()
String name = entry != null ? entry.getImageName() : getCurrentServer().getMetadata().getName()
String code = name.tokenize('.')[0].trim()
if (!(code ==~ /T[A-Z0-9]{8}/))
    throw new IllegalArgumentException("not a benchmark tile name: '${name}' -- run this on the annotation project only")

def all = getAnnotationObjects()
def objects = all.findAll { it.getROI() != null && it.getROI().isArea() }
def dir = new File(OUT_DIR, 'truth')
dir.mkdirs()
def out = new File(dir, code + '.geojson')
if (objects.isEmpty())
    out.text = '{"type": "FeatureCollection", "features": []}\n'
else
    exportObjectsToGeoJson(objects, out.getAbsolutePath(), "FEATURE_COLLECTION", "EXCLUDE_MEASUREMENTS")
int ambiguous = objects.count { it.getPathClass()?.getName() == 'Ambiguous' } as int
println "${code}\t${objects.size()} objects (${ambiguous} Ambiguous)\t${all.size() - objects.size()} non-area skipped\t-> ${out}"
