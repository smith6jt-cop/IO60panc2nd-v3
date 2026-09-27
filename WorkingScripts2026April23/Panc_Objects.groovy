import qupath.lib.objects.PathObjects
import qupath.lib.roi.GeometryTools
import qupath.lib.roi.RoiTools

setImageType('FLUORESCENCE')
setPixelSizeMicrons(0.5078, 0.5078)

// ---- Edge/inner analysis area (µm) ----
double edgeMarginUm = 5.0

// ---- Size-adaptive erosion parameters ----
// ADJUSTED: More conservative to avoid losing large islets
double SMALL_ISLET_AREA_UM2 = 5000.0      // Islets < 1000 µm² considered small
double MEDIUM_ISLET_AREA_UM2 = 20000.0    // Islets < 30000 µm² considered medium
// Large islets >= 30000 µm²

double SMALL_EROSION_UM = 1           // Very gentle for small islets (REDUCED from 1.0)
double MEDIUM_EROSION_UM = 2          // Moderate for medium islets (REDUCED from 3.0)
double LARGE_EROSION_UM = 3            // Conservative for large islets (REDUCED from 6.0)

double DILATION_FACTOR = 0.85            // Restore 85% of erosion (INCREASED from 0.8)

// Calibration
def server = getCurrentServer()
def cal = server?.getPixelCalibration()
double sx = (cal?.getPixelWidthMicrons())  ?: (cal?.getAveragedPixelSizeMicrons()) ?: 1.0
double sy = (cal?.getPixelHeightMicrons()) ?: (cal?.getAveragedPixelSizeMicrons()) ?: 1.0
double marginPxX = edgeMarginUm / sx
double marginPxY = edgeMarginUm / sy

// Map to store erosion amounts (FIXED: use external map instead of setMetadata)
def erosionMap = [:]

// Per-plane helpers
def innerGeomForPlane = { plane ->
    double w = server.getWidth()
    double h = server.getHeight()
    double iw = Math.max(0.0, w - 2*marginPxX)
    double ih = Math.max(0.0, h - 2*marginPxY)
    def innerROI = ROIs.createRectangleROI(marginPxX, marginPxY, iw, ih, plane)
    innerROI.getGeometry()
}

// Clip a collection of objects to the inner area
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

// ---- Size-adaptive erosion/dilation ----
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
        double areaUm2 = geom.getArea() * sx * sy  // Convert to µm²
        
        // Determine erosion amount based on size
        double erosionUm
        String sizeCategory
        if (areaUm2 < SMALL_ISLET_AREA_UM2) {
            erosionUm = SMALL_EROSION_UM
            sizeCategory = "small"
            smallCount++
        } else if (areaUm2 < MEDIUM_ISLET_AREA_UM2) {
            erosionUm = MEDIUM_EROSION_UM
            sizeCategory = "medium"
            mediumCount++
        } else {
            erosionUm = LARGE_EROSION_UM
            sizeCategory = "large"
            largeCount++
        }
        
        try {
            // Apply size-specific erosion
            def eroded = geom.buffer(-erosionUm)
            
            // Check if completely eroded away
            if (eroded.isEmpty() || eroded.getArea() < 1.0) {
                println "  ${sizeCategory} islet (${String.format('%.0f', areaUm2)} µm²) eroded away with ${erosionUm} µm"
                toRemove << obj
                erodedAwayCount++
                return
            }
            
            // Check if split into multiple parts
            if (eroded.getNumGeometries() > 1) {
                println "  ${sizeCategory} cluster (${String.format('%.0f', areaUm2)} µm²) split into ${eroded.getNumGeometries()} parts"
                toRemove << obj
                splitCount++
                
                // Create separate objects for each part
                for (int i = 0; i < eroded.getNumGeometries(); i++) {
                    def part = eroded.getGeometryN(i)
                    if (!part.isEmpty() && part.getArea() >= 1.0) {
                        def newROI = GeometryTools.geometryToROI(part, plane)
                        def newObj = PathObjects.createAnnotationObject(newROI, obj.getPathClass())
                        if (obj.getName()) newObj.setName(obj.getName())
                        // FIXED: Store erosion amount in external map
                        erosionMap[newObj] = erosionUm
                        toAdd << newObj
                    }
                }
            } else {
                // Single part - update in place
                def newROI = GeometryTools.geometryToROI(eroded, plane)
                obj.setROI(newROI)
                // FIXED: Store erosion amount in external map
                erosionMap[obj] = erosionUm
            }
        } catch (Exception e) {
            println "  ERROR eroding ${sizeCategory} islet: ${e.message}"
        }
    }
    
    println "\nErosion summary:"
    println "  Small islets (<${SMALL_ISLET_AREA_UM2} µm²): ${smallCount} (eroded ${SMALL_EROSION_UM} µm)"
    println "  Medium islets (<${MEDIUM_ISLET_AREA_UM2} µm²): ${mediumCount} (eroded ${MEDIUM_EROSION_UM} µm)"
    println "  Large islets/clusters (≥${MEDIUM_ISLET_AREA_UM2} µm²): ${largeCount} (eroded ${LARGE_EROSION_UM} µm)"
    println "  Split into parts: ${splitCount}"
    println "  Eroded away: ${erodedAwayCount}"
    
    if (!toRemove.isEmpty()) {
        removeObjects(toRemove, true)
    }
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
        
        // FIXED: Get erosion amount from external map
        double erosionAmount = erosionMap.get(obj) ?: MEDIUM_EROSION_UM
        
        // Dilate by a fraction of the erosion to maintain separation
        double dilationUm = erosionAmount * DILATION_FACTOR
        
        def plane = roi.getImagePlane()
        def geom = roi.getGeometry()
        
        try {
            def dilated = geom.buffer(dilationUm)
            
            if (!dilated.isEmpty()) {
                def newROI = GeometryTools.geometryToROI(dilated, plane)
                obj.setROI(newROI)
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
    'GCG': ['GCG','Glucagon','Gcg'],
    'INS': ['INS','Insulin','Ins']
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

// --- Run classifiers & refinement ---
//createAnnotationsFromPixelClassifier("Tissue", 1000.0, -1.0)
//clipObjectsToInner(getAnnotationObjects().findAll { it.getPathClass()?.toString() == 'Tissue' })
//
//selectObjectsByClassification("Tissue")
//runPlugin('qupath.lib.plugins.objects.RefineAnnotationsPlugin',
//    '{"minFragmentSizeMicrons":10000.0,"maxHoleSizeMicrons":-1}')

selectObjectsByClassification("Tissue")
createAnnotationsFromPixelClassifier("Islets", 300.0, -1.0, "SPLIT")

clipObjectsToInner(getAnnotationObjects().findAll { it.getPathClass()?.toString() == 'Islets' })

// ---- SIZE-ADAPTIVE ERODE-DILATE ----
def isletAnnotations = getAnnotationObjects().findAll { it.getPathClass()?.toString() == 'Islet' }
println "\n=== Applying size-adaptive erosion-dilation to separate islets ==="
println "Found ${isletAnnotations.size()} Islet annotations before processing"

if (!isletAnnotations.isEmpty()) {
    try {
        // Size-adaptive erosion
        sizeAdaptiveErode(isletAnnotations)
        
        // Refresh collection
        isletAnnotations = getAnnotationObjects().findAll { it.getPathClass()?.toString() == 'Islet' }
        println "\nAfter erosion: ${isletAnnotations.size()} Islet annotations"
        
        // Size-adaptive dilation
        if (!isletAnnotations.isEmpty()) {
            sizeAdaptiveDilate(isletAnnotations)
        }
    } catch (Exception e) {
        println "ERROR during size-adaptive processing: ${e.message}"
        e.printStackTrace()
    }
}

println "\n=== Size-adaptive processing complete ===\n"

// Continue with refinement
selectObjectsByClassification("Islet")
runPlugin('qupath.lib.plugins.objects.RefineAnnotationsPlugin',
    '{"minFragmentSizeMicrons":250.0,"maxHoleSizeMicrons":-1.0}')

print "Done."