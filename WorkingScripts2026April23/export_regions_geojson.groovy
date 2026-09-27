/**
 * Export Follicle and PALS annotations as a GeoJSON FeatureCollection.
 * Runs per-image when invoked with QuPath script -p (no --save needed).
 *
 *   QuPath script -p project.qpproj scripts/export_regions_geojson.groovy
 */

import qupath.lib.io.GsonTools
import com.google.gson.stream.JsonWriter

def imageData = getCurrentImageData()
def hierarchy = imageData.getHierarchy()
def annotations = hierarchy.getAnnotationObjects()

def filtered = annotations.findAll { anno ->
    def cls = anno.getPathClass()?.toString()
    cls == 'Follicle' || cls == 'PALS'
}

if (filtered.isEmpty()) {
    println "No Follicle/PALS annotations — skipping"
    return
}

def nFollicle = filtered.count { it.getPathClass()?.toString() == 'Follicle' }
def nPALS = filtered.count { it.getPathClass()?.toString() == 'PALS' }

def project = getProject()
def outputDir = new File(project.getPath().getParent().toFile(), "analysis/geojson")
outputDir.mkdirs()

def imageName = getProjectEntry().getImageName()
def baseName = imageName
    .replaceAll(/\.ome\.tiff?$/, '')
    .replaceAll(/\.tiff?$/, '')
def outputFile = new File(outputDir, "${baseName}.geojson")

def gson = GsonTools.getInstance(true)
def writer = new JsonWriter(new java.io.FileWriter(outputFile))
writer.setIndent("  ")
writer.beginObject()
writer.name("type").value("FeatureCollection")
writer.name("features")
writer.beginArray()
for (obj in filtered) {
    gson.toJson(obj, obj.getClass(), writer)
}
writer.endArray()
writer.endObject()
writer.close()

println "${nFollicle} Follicle + ${nPALS} PALS → ${outputFile.getName()}"
