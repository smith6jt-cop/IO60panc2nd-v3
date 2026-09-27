/*
 * 1_panc_objects.groovy
 *
 * Workflow step 2 after creating Tissue annotation and, if necessary, "Islets" pixel classifier:
 *   - sets image type and pixel calibration
 *   - normalizes channel names (strips "(Cxx)" suffixes, applies alias map)
 *   - creates Islet annotations via the Islets pixel classifier
 *   - applies size-adaptive erode/dilate to separate touching islets
 *   - refines small fragments
 *
 * Tissue annotations are assumed to already exist (provided by the user).
 */

evaluate(new File(buildFilePath(PROJECT_BASE_DIR, 'scripts', 'lib', 'common.groovy')))

import qupath.lib.objects.PathObjects
import qupath.lib.roi.GeometryTools
import qupath.lib.roi.RoiTools
import qupath.lib.roi.ROIs
import qupath.lib.images.ImageData
import qupath.lib.objects.hierarchy.PathObjectHierarchy
import qupath.opencv.ml.pixel.PixelClassifierTools

setImageType('FLUORESCENCE')
setPixelSizeMicrons(0.5078, 0.5078)

// ---- Edge/inner analysis area (um) ----
double edgeMarginUm = 5.0

// ---- Size-adaptive erosion parameters ----
double SMALL_ISLET_AREA_UM2  = 2000.0
double MEDIUM_ISLET_AREA_UM2 = 10000.0
double SMALL_EROSION_UM      = 1
double MEDIUM_EROSION_UM     = 2
double LARGE_EROSION_UM      = 3
double DILATION_FACTOR       = 0.85

def server = getCurrentServer()
def cal = server?.getPixelCalibration()
double sx = (cal?.getPixelWidthMicrons())  ?: (cal?.getAveragedPixelSizeMicrons()) ?: 1.0
double sy = (cal?.getPixelHeightMicrons()) ?: (cal?.getAveragedPixelSizeMicrons()) ?: 1.0
double marginPxX = edgeMarginUm / sx
double pxUm = Math.sqrt(sx * sy)   // um per pixel; ROI geometries are in pixel units
double marginPxY = edgeMarginUm / sy

def erosionMap = [:]

def innerGeomForPlane = { plane ->
    double w = server.getWidth()
    double h = server.getHeight()
    double iw = Math.max(0.0, w - 2*marginPxX)
    double ih = Math.max(0.0, h - 2*marginPxY)
    def innerROI = ROIs.createRectangleROI(marginPxX, marginPxY, iw, ih, plane)
    innerROI.getGeometry()
}

def clipObjectsToInner = { Collection objs ->
    def toRemove = []
    objs.each { o ->
        def roi = o.getROI(); if (roi == null) return
        def plane = roi.getImagePlane()
        def inner = innerGeomForPlane(plane)
        def inter = roi.getGeometry().intersection(inner)
        if (inter.isEmpty()) {
            toRemove << o
        } else if (!inter.equalsTopo(roi.getGeometry())) {
            o.setROI(GeometryTools.geometryToROI(inter, plane))
        }
    }
    if (!toRemove.isEmpty()) removeObjects(toRemove, true)
}

def sizeAdaptiveErode = { Collection objs ->
    println "\n=== Size-Adaptive Erosion ==="
    println "Analyzing ${objs.size()} annotations..."

    def toRemove = []
    def toAdd = []
    int smallCount = 0, mediumCount = 0, largeCount = 0
    int splitCount = 0, erodedAwayCount = 0

    objs.each { obj ->
        def roi = obj.getROI()
        if (roi == null) return

        def plane = roi.getImagePlane()
        def geom = roi.getGeometry()
        double areaUm2 = geom.getArea() * sx * sy

        double erosionUm
        String sizeCategory
        if (areaUm2 < SMALL_ISLET_AREA_UM2) {
            erosionUm = SMALL_EROSION_UM; sizeCategory = "small"; smallCount++
        } else if (areaUm2 < MEDIUM_ISLET_AREA_UM2) {
            erosionUm = MEDIUM_EROSION_UM; sizeCategory = "medium"; mediumCount++
        } else {
            erosionUm = LARGE_EROSION_UM; sizeCategory = "large"; largeCount++
        }

        try {
            def eroded = geom.buffer(-erosionUm / pxUm)
            if (eroded.isEmpty() || eroded.getArea() < 1.0) {
                println "  ${sizeCategory} islet (${String.format('%.0f', areaUm2)} um^2) eroded away with ${erosionUm} um"
                toRemove << obj
                erodedAwayCount++
                return
            }
            if (eroded.getNumGeometries() > 1) {
                println "  ${sizeCategory} cluster (${String.format('%.0f', areaUm2)} um^2) split into ${eroded.getNumGeometries()} parts"
                toRemove << obj
                splitCount++
                for (int i = 0; i < eroded.getNumGeometries(); i++) {
                    def part = eroded.getGeometryN(i)
                    if (!part.isEmpty() && part.getArea() >= 1.0) {
                        def newROI = GeometryTools.geometryToROI(part, plane)
                        def newObj = PathObjects.createAnnotationObject(newROI, obj.getPathClass())
                        if (obj.getName()) newObj.setName(obj.getName())
                        erosionMap[newObj] = erosionUm
                        toAdd << newObj
                    }
                }
            } else {
                def newROI = GeometryTools.geometryToROI(eroded, plane)
                obj.setROI(newROI)
                erosionMap[obj] = erosionUm
            }
        } catch (Exception e) {
            println "  ERROR eroding ${sizeCategory} islet: ${e.message}"
        }
    }

    println "\nErosion summary:"
    println "  Small islets (<${SMALL_ISLET_AREA_UM2} um^2): ${smallCount} (eroded ${SMALL_EROSION_UM} um)"
    println "  Medium islets (<${MEDIUM_ISLET_AREA_UM2} um^2): ${mediumCount} (eroded ${MEDIUM_EROSION_UM} um)"
    println "  Large islets/clusters (>=${MEDIUM_ISLET_AREA_UM2} um^2): ${largeCount} (eroded ${LARGE_EROSION_UM} um)"
    println "  Split into parts: ${splitCount}"
    println "  Eroded away: ${erodedAwayCount}"

    if (!toRemove.isEmpty()) removeObjects(toRemove, true)
    if (!toAdd.isEmpty()) {
        addObjects(toAdd)
        println "  Added ${toAdd.size()} new objects from splits"
    }
}

def sizeAdaptiveDilate = { Collection objs ->
    println "\n=== Size-Adaptive Dilation ==="
    println "Dilating ${objs.size()} annotations..."

    int dilatedCount = 0
    objs.each { obj ->
        def roi = obj.getROI()
        if (roi == null) return

        double erosionAmount = erosionMap.get(obj) ?: MEDIUM_EROSION_UM
        double dilationUm = erosionAmount * DILATION_FACTOR

        def plane = roi.getImagePlane()
        def geom = roi.getGeometry()

        try {
            def dilated = geom.buffer(dilationUm / pxUm)
            if (!dilated.isEmpty()) {
                obj.setROI(GeometryTools.geometryToROI(dilated, plane))
                dilatedCount++
            }
        } catch (Exception e) {
            println "  ERROR dilating object: ${e.message}"
        }
    }

    println "Dilated ${dilatedCount} annotations (restored ${DILATION_FACTOR*100}% of erosion)"
}

// --- Channel normalization ---
def meta   = server.getMetadata()
def names0 = meta.getChannels().collect { it.getName() }
println "Channels before normalize:"
names0.eachWithIndex { n,i -> println String.format("  C%02d: %s", i+1, n) }

def names1 = names0.collect { n ->
    if (n == null) return null
    def s = n
    s = s.replaceAll('\\s*\\(C\\d+\\)\\s*$', '')
    s = s.replaceAll('(?i)\\.(tif|tiff)\\s*$', '')
    s = s.trim()
    return s
}
setChannelNames(*names1 as String[])

def aliasMap = [
    'SST': ['SST','Somatostatin','Stt'],
    'GCG': ['GCG','Glucagon','Gcg', 'GCG', 'Opal 520'],
    'INS': ['INS','Insulin','Ins', 'INS', 'Opal 690'],
    'CD99': ['CD99']
]
def names2 = getCurrentServer().getMetadata().getChannels().collect { it.getName() } as List<String>
aliasMap.each { canonical, aliases ->
    if (names2.contains(canonical)) return
    def idx = names2.findIndexOf { it != null && aliases.any { a -> it.equalsIgnoreCase(a) } }
    if (idx >= 0) {
        names2[idx] = canonical
        setChannelNames(*names2 as String[])
    }
}

println "Channels after normalize:"
getCurrentServer().getMetadata().getChannels().collect{ it.getName() }
    .eachWithIndex { n,i -> println String.format("  C%02d: %s", i+1, n) }

// --- Islet classifier on existing Tissue ---

// Run the classifier through a private Bio-Formats server instance. In QuPath <= 0.7.0 an
// interrupted tile read (live classifier/threshold preview, viewer tile cancellation) closes the
// FileChannel of a pooled Bio-Formats reader, and that reader is reused, so every later read on the
// viewer's server fails. A separate instance has its own reader pool and pending-tile map.
// Re-run safety: createAnnotationsFromPixelClassifier (without DELETE_EXISTING) adds new islets next to
// any existing ones, and the erode/dilate step below would then erode the old islets a second time.
// Remove previous Islet annotations only; keepChildren=true moves any parented detections up to Tissue.
// (DELETE_EXISTING is not an option: it clears ALL children of Tissue, detections included.)
def oldIslets = getAnnotationObjects().findAll { it.getPathClass()?.toString() == 'Islet' }
if (!oldIslets.isEmpty()) {
    println "Removing ${oldIslets.size()} existing Islet annotations before re-classifying"
    removeObjects(oldIslets, true)
}
selectObjectsByClassification("Tissue")
def imgData = getCurrentImageData()
def isletClassifier = loadPixelClassifier("Islets2")
def privateServer = imgData.getServerBuilder().build()
try {
    privateServer.setMetadata(imgData.getServer().getMetadata())   // carry normalized channel names + calibration
    def privateData = new ImageData(privateServer, new PathObjectHierarchy(), imgData.getImageType())
    def classServer = PixelClassifierTools.createPixelClassificationServer(privateData, isletClassifier)
    // Objects go into the real hierarchy, constrained to the selected Tissue annotations
    PixelClassifierTools.createAnnotationsFromPixelClassifier(
        imgData.getHierarchy(), classServer, 50.0, -1.0,
        PixelClassifierTools.CreateObjectOptions.SPLIT)
} finally {
    privateServer.close()
}

clipObjectsToInner(getAnnotationObjects().findAll { it.getPathClass()?.toString() == 'Islet' })

// ---- SIZE-ADAPTIVE ERODE-DILATE ----
def isletAnnotations = getAnnotationObjects().findAll { it.getPathClass()?.toString() == 'Islet' }
println "\n=== Applying size-adaptive erosion-dilation to separate islets ==="
println "Found ${isletAnnotations.size()} Islet annotations before processing"

if (!isletAnnotations.isEmpty()) {
    try {
        sizeAdaptiveErode(isletAnnotations)
        isletAnnotations = getAnnotationObjects().findAll { it.getPathClass()?.toString() == 'Islet' }
        println "\nAfter erosion: ${isletAnnotations.size()} Islet annotations"
        if (!isletAnnotations.isEmpty()) {
            sizeAdaptiveDilate(isletAnnotations)
        }
    } catch (Exception e) {
        println "ERROR during size-adaptive processing: ${e.message}"
        e.printStackTrace()
    }
}

println "\n=== Size-adaptive processing complete ===\n"

selectObjectsByClassification("Islet")
runPlugin('qupath.lib.plugins.objects.RefineAnnotationsPlugin',
    '{"minFragmentSizeMicrons":50.0,"maxHoleSizeMicrons":-1.0}')

print "Done."
