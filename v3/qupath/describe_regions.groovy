/*
 * describe_regions.groovy -- what an image's data file holds about its regions and its cells. Reads only.
 *
 * Written for the moment a data file has been edited outside the driver (the maintainer redraws a Tissue or pLN
 * annotation in the QuPath window): before anything is re-run, say what the file holds now.
 *
 *   annotations     count per class
 *   regions         every Tissue and pLN annotation: object ID, name, area, pieces, holes
 *   overlaps        for every pair of regions the area they share (the maintainer, 2026-10-02: "no Tissue annotation
 *                   should" overlap a pLN)
 *   detections      how many; how many per class of parent in the hierarchy; how many have their centroid inside
 *                   each region, inside more than one, and inside none
 *
 * and writes the regions themselves as GeoJSON beside the report, so they can be compared with an earlier export.
 *
 *   QuPath script -p <PRIVATE project copy> -i "<image>" -a <report.json> -a <regions.geojson> describe_regions.groovy
 */

import org.locationtech.jts.geom.Geometry
import org.locationtech.jts.geom.Polygon
import org.locationtech.jts.geom.prep.PreparedGeometryFactory
import qupath.lib.io.GsonTools   // groovy-json is not bundled with QuPath 0.7

def a = (binding.hasVariable('args') && args != null) ? (args as List).collect { it.toString() } : []
if (a.size() < 2)
    throw new IllegalArgumentException('usage: describe_regions.groovy <report.json> <regions.geojson>')
def reportFile = new File(a[0])
def regionsFile = new File(a[1])
Set<String> REGION_CLASSES = ['Tissue', 'pLN'] as Set

def server = getCurrentServer()
double px = server.getPixelCalibration().getAveragedPixelSizeMicrons()
def clsOf = { o -> o.getPathClass()?.toString() ?: 'Unclassified' }
def anns = getAnnotationObjects()
def regions = anns.findAll { clsOf(it) in REGION_CLASSES }.sort { clsOf(it) + (it.getName() ?: '') + it.getID() }
def dets = getDetectionObjects()

def holesOf = { Geometry g ->
    int n = 0
    for (int i = 0; i < g.getNumGeometries(); i++) {
        def p = g.getGeometryN(i)
        if (p instanceof Polygon) n += ((Polygon) p).getNumInteriorRing()
    }
    n
}
def geoms = regions.collect { it.getROI().getGeometry() }
def regionRows = []
regions.eachWithIndex { r, i ->
    def g = geoms[i]
    regionRows << [index: i, 'class': clsOf(r), name: r.getName(), id: r.getID().toString(),
                   area_um2: g.getArea() * px * px, pieces: g.getNumGeometries(), holes: holesOf(g),
                   vertices: g.getNumPoints(), valid: g.isValid(), locked: r.isLocked()]
}
def overlaps = []
for (int i = 0; i < regions.size(); i++)
    for (int j = i + 1; j < regions.size(); j++) {
        if (!geoms[i].getEnvelopeInternal().intersects(geoms[j].getEnvelopeInternal())) continue
        double shared = 0
        try { shared = geoms[i].intersection(geoms[j]).getArea() * px * px }
        catch (Throwable t) { shared = geoms[i].buffer(0).intersection(geoms[j].buffer(0)).getArea() * px * px }
        overlaps << [a: i, b: j, a_class: clsOf(regions[i]), b_class: clsOf(regions[j]), shared_um2: shared,
                     gap_um: shared > 0 ? 0.0d : geoms[i].distance(geoms[j]) * px]
    }

// detections: parent class in the hierarchy, and which region the centroid is in
def byParent = new TreeMap(dets.countBy { d ->
    def p = d.getParent()
    p == null ? 'NONE' : (p.isRootObject() ? 'ROOT' : (p.isAnnotation() ? "Annotation (${clsOf(p)})".toString() : p.getClass().getSimpleName()))
})
def prepared = geoms.collect { PreparedGeometryFactory.prepare(it) }
def gf = new org.locationtech.jts.geom.GeometryFactory()
int[] inside = new int[regions.size()]
java.util.concurrent.atomic.AtomicInteger none = new java.util.concurrent.atomic.AtomicInteger()
java.util.concurrent.atomic.AtomicInteger several = new java.util.concurrent.atomic.AtomicInteger()
def counts = (0..<regions.size()).collect { new java.util.concurrent.atomic.AtomicInteger() }
def parentDisagrees = new java.util.concurrent.atomic.AtomicInteger()
dets.parallelStream().forEach { d ->
    def roi = d.getROI()
    def pt = gf.createPoint(new org.locationtech.jts.geom.Coordinate(roi.getCentroidX(), roi.getCentroidY()))
    int n = 0
    String in = null
    for (int i = 0; i < prepared.size(); i++)
        if (prepared[i].contains(pt)) { counts[i].incrementAndGet(); n++; in = clsOf(regions[i]) }
    if (n == 0) none.incrementAndGet()
    if (n > 1) several.incrementAndGet()
    def p = d.getParent()
    String parentClass = (p != null && p.isAnnotation()) ? clsOf(p) : null
    if (n == 1 && parentClass in REGION_CLASSES && parentClass != in) parentDisagrees.incrementAndGet()
}
regionRows.eachWithIndex { row, i -> row.detections_with_centroid_inside = counts[i].get() }

exportObjectsToGeoJson(regions, regionsFile.getAbsolutePath(), 'FEATURE_COLLECTION')
def report = [
        image: getProjectEntry()?.getImageName(), entry_id: getProjectEntry()?.getID(), pixel_um: px,
        annotations_by_class: new TreeMap(anns.countBy { clsOf(it) }),
        regions: regionRows, overlaps: overlaps,
        detections: dets.size(), detections_by_parent: byParent,
        detections_with_centroid_in_no_region: none.get(), detections_with_centroid_in_several_regions: several.get(),
        detections_whose_region_parent_is_not_the_region_they_are_in: parentDisagrees.get(),
        regions_geojson: regionsFile.getAbsolutePath(),
]
reportFile.getParentFile()?.mkdirs()
reportFile.text = GsonTools.getInstance(true).toJson(report)
println "[describe_regions] ${report.annotations_by_class}; ${dets.size()} detections; overlaps ${overlaps.findAll { it.shared_um2 > 0 }.size()}; -> ${reportFile}"
