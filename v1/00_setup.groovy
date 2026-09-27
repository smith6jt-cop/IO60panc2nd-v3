/*
 * 00_setup.groovy
 *
 * Stage 00: per-image setup.
 *   - sets image type and pixel calibration
 *   - normalizes channel names (strips "(Cxx)" suffixes, applies alias map)
 *   - creates Tissue annotation (optional, uncomment if needed)
 *   - creates Islet annotations via pixel classifier
 *   - applies size-adaptive erode/dilate to separate touching islets
 *   - refines small fragments
 *
 * Run via: Run -> Run for project (select all 15 images).
 */

evaluate(new File(buildFilePath(PROJECT_BASE_DIR, 'scripts', 'lib', 'common.groovy')))

import qupath.lib.objects.PathObjects
import qupath.lib.roi.GeometryTools
import qupath.lib.roi.RoiTools
import qupath.lib.roi.ROIs

// ==================== PARAMETERS ====================
double EDGE_MARGIN_UM        = 5.0
double SMALL_ISLET_AREA_UM2  = 4000.0
double MEDIUM_ISLET_AREA_UM2 = 10000.0
double SMALL_EROSION_UM      = 1.0
double MEDIUM_EROSION_UM     = 2.0
double LARGE_EROSION_UM      = 3.0
double DILATION_FACTOR       = 0.85
String TISSUE_CLASSIFIER     = 'Tissue'     // set to null to skip Tissue creation
String ISLET_CLASSIFIER      = 'Islets'
double ISLET_MIN_AREA_UM2    = 100.0
double TISSUE_MIN_AREA_UM2   = 1000.0
double TISSUE_MIN_FRAGMENT   = 10000.0
double ISLET_MIN_FRAGMENT    = 50.0
// ====================================================

setImageType('FLUORESCENCE')
setPixelSizeMicrons(0.5078, 0.5078)

def server = getCurrentServer()
def cal = server.getPixelCalibration()
double sx = (cal?.getPixelWidthMicrons())  ?: (cal?.getAveragedPixelSizeMicrons()) ?: 1.0
double sy = (cal?.getPixelHeightMicrons()) ?: (cal?.getAveragedPixelSizeMicrons()) ?: 1.0
double marginPxX = EDGE_MARGIN_UM / sx
double marginPxY = EDGE_MARGIN_UM / sy

// ---------- Channel normalization ----------
def meta   = server.getMetadata()
def names0 = meta.getChannels().collect { it.getName() }
println "[00] Channels before normalize:"
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

// Alias map: canonical -> [alternative spellings].
// Extended for the object-classifier channels (M2Gal3, CD34, Iba1).
def aliasMap = [
    'SST'   : ['SST','Somatostatin','Stt'],
    'GCG'   : ['GCG','Glucagon','Gcg'],
    'INS'   : ['INS','Insulin','Ins'],
    'M2Gal3': ['M2Gal3','MacGal3','Gal3','Galectin3','Galectin-3','LGALS3'],
    'CD34'  : ['CD34'],
    'CD99'  : ['CD99'],
    'Iba1'  : ['Iba1','IBA1','AIF1','Aif1']
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

println "[00] Channels after normalize:"
getCurrentServer().getMetadata().getChannels().collect{ it.getName() }
    .eachWithIndex { n,i -> println String.format("  C%02d: %s", i+1, n) }

// Verify classifier channels are present
def finalNames = getCurrentServer().getMetadata().getChannels().collect { it.getName() }
def missingClf = OBJECT_CLF_CHANNELS.findAll { !finalNames.contains(it) }
if (!missingClf.isEmpty()) {
    println "[00] WARN: object-classifier channels missing after alias: ${missingClf}"
    println "[00] Check the alias map above and re-run."
}

// ---------- Inner-area clipping ----------
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

// ---------- Size-adaptive erode/dilate ----------
def erosionMap = [:]

def sizeAdaptiveErode = { Collection objs ->
    println "\n[00] Size-adaptive erosion on ${objs.size()} annotations..."
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
        if (areaUm2 < SMALL_ISLET_AREA_UM2)        { erosionUm = SMALL_EROSION_UM;  sizeCategory = 'small';  smallCount++ }
        else if (areaUm2 < MEDIUM_ISLET_AREA_UM2)  { erosionUm = MEDIUM_EROSION_UM; sizeCategory = 'medium'; mediumCount++ }
        else                                       { erosionUm = LARGE_EROSION_UM;  sizeCategory = 'large';  largeCount++ }

        try {
            def eroded = geom.buffer(-erosionUm)
            if (eroded.isEmpty() || eroded.getArea() < 1.0) {
                toRemove << obj; erodedAwayCount++; return
            }
            if (eroded.getNumGeometries() > 1) {
                toRemove << obj; splitCount++
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
                obj.setROI(GeometryTools.geometryToROI(eroded, plane))
                erosionMap[obj] = erosionUm
            }
        } catch (Exception e) {
            println "  ERROR eroding ${sizeCategory}: ${e.message}"
        }
    }

    println "  small<${SMALL_ISLET_AREA_UM2}: ${smallCount}  medium<${MEDIUM_ISLET_AREA_UM2}: ${mediumCount}  large: ${largeCount}"
    println "  split: ${splitCount}  eroded-away: ${erodedAwayCount}"
    if (!toRemove.isEmpty()) removeObjects(toRemove, true)
    if (!toAdd.isEmpty())    addObjects(toAdd)
}

def sizeAdaptiveDilate = { Collection objs ->
    println "\n[00] Size-adaptive dilation on ${objs.size()} annotations..."
    int dilatedCount = 0
    objs.each { obj ->
        def roi = obj.getROI(); if (roi == null) return
        double erosionAmount = erosionMap.get(obj) ?: MEDIUM_EROSION_UM
        double dilationUm = erosionAmount * DILATION_FACTOR
        def plane = roi.getImagePlane()
        def geom = roi.getGeometry()
        try {
            def dilated = geom.buffer(dilationUm)
            if (!dilated.isEmpty()) {
                obj.setROI(GeometryTools.geometryToROI(dilated, plane))
                dilatedCount++
            }
        } catch (Exception e) {
            println "  ERROR dilating: ${e.message}"
        }
    }
    println "  dilated ${dilatedCount} (restored ${DILATION_FACTOR*100}% of erosion)"
}

// ---------- Tissue (optional) ----------
if (TISSUE_CLASSIFIER != null) {
    def existingTissue = getAnnotationObjects().findAll { it.getPathClass()?.toString() == CLASS_TISSUE }
    if (existingTissue.isEmpty()) {
        println "[00] Creating ${CLASS_TISSUE} via classifier '${TISSUE_CLASSIFIER}'"
        createAnnotationsFromPixelClassifier(TISSUE_CLASSIFIER, TISSUE_MIN_AREA_UM2, -1.0)
        clipObjectsToInner(getAnnotationObjects().findAll { it.getPathClass()?.toString() == CLASS_TISSUE })
        selectObjectsByClassification(CLASS_TISSUE)
        runPlugin('qupath.lib.plugins.objects.RefineAnnotationsPlugin',
            "{\"minFragmentSizeMicrons\":${TISSUE_MIN_FRAGMENT},\"maxHoleSizeMicrons\":-1}")
    } else {
        println "[00] ${existingTissue.size()} existing ${CLASS_TISSUE} annotations - skipping creation."
    }
}

// ---------- Islets ----------
def tissueAnns = getAnnotationObjects().findAll { it.getPathClass()?.toString() == CLASS_TISSUE }
if (tissueAnns.isEmpty()) {
    println "[00] ERROR: no ${CLASS_TISSUE} annotations - cannot create islets."
    return
}

// Remove existing Islet annotations to keep this stage idempotent
def existingIslets = getAnnotationObjects().findAll { it.getPathClass()?.toString() == CLASS_ISLET }
if (!existingIslets.isEmpty()) {
    println "[00] Removing ${existingIslets.size()} existing ${CLASS_ISLET} annotations."
    removeObjects(existingIslets, true)
}

selectObjectsByClassification(CLASS_TISSUE)
createAnnotationsFromPixelClassifier(ISLET_CLASSIFIER, ISLET_MIN_AREA_UM2, -1.0, "SPLIT")
clipObjectsToInner(getAnnotationObjects().findAll { it.getPathClass()?.toString() == CLASS_ISLET })

def isletAnnotations = getAnnotationObjects().findAll { it.getPathClass()?.toString() == CLASS_ISLET }
println "[00] ${isletAnnotations.size()} ${CLASS_ISLET} annotations before erode/dilate"

if (!isletAnnotations.isEmpty()) {
    try {
        sizeAdaptiveErode(isletAnnotations)
        isletAnnotations = getAnnotationObjects().findAll { it.getPathClass()?.toString() == CLASS_ISLET }
        if (!isletAnnotations.isEmpty()) {
            sizeAdaptiveDilate(isletAnnotations)
        }
    } catch (Exception e) {
        println "[00] ERROR during erode/dilate: ${e.message}"
    }
}

selectObjectsByClassification(CLASS_ISLET)
runPlugin('qupath.lib.plugins.objects.RefineAnnotationsPlugin',
    "{\"minFragmentSizeMicrons\":${ISLET_MIN_FRAGMENT},\"maxHoleSizeMicrons\":-1.0}")
addShapeMeasurements("AREA", "LENGTH", "CIRCULARITY", "SOLIDITY", "MAX_DIAMETER", "MIN_DIAMETER", "NUCLEUS_CELL_RATIO")
fireHierarchyUpdate()
println "[00] Stage 00 complete. Islet count: ${getAnnotationObjects().findAll{ it.getPathClass()?.toString() == CLASS_ISLET }.size()}"
