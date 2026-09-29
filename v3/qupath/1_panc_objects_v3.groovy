/*
 * 1_panc_objects_v3.groovy   (v3 copy of scripts/1_panc_objects.groovy; the original is unchanged)
 *
 * Workflow step 2 after the Tissue annotation exists:
 *   - sets image type and pixel calibration
 *   - normalizes channel names (strips "(Cxx)" suffixes, applies alias map)
 *   - creates endocrine objects inside Tissue with one of two methods (ISLET_METHOD):
 *
 *   'Islets2'      legacy path, unchanged: the trained Islets2 pixel classifier at 1.0156 um, min area
 *                  50 um^2, SPLIT, clip 5 um inside the image edge, size-adaptive erode/dilate, refine.
 *                  Kept only as a benchmark baseline.
 *
 *   'drotar_cd99'  (default) Drotar 2025 / ADIIDAS endocrine objects, CD99 primary (scripts/v3/ISLETS.md):
 *     1. per-image anchors and clamp models from results/v3/islet_anchors.csv (scripts/v3/islet_anchors.py,
 *        method 2.x), estimated in THIS script's pixel space (RESOLUTION_UM, Gaussian SIGMA_UM) on the same
 *        analysis area. The CSV row carries a hash of the Tissue/pLN/Artifact annotations it used; the hash is
 *        recomputed here and a mismatch (or a different SIGMA_UM / RESOLUTION_UM) stops the run
 *        (ANCHOR_CHECK='error'; 'warn' only reports it).
 *     2. pixels (ImageOps, at RESOLUTION_UM): x -> Gaussian(SIGMA_UM); channel c is positive when
 *          x_c >= max(tau * anchor_c,  clamp_mu_c + FLOOR_K * clamp_sigma_c)
 *        tau = TAU_CD99 for CD99 and TAU_H for INS, GCG, SST. The second term is the per-image FLOOR: the
 *        scanner clamps over-subtracted background to 0 (INS: up to ~80% of level-0 tissue pixels), which
 *        shrinks the half-max anchor, so tau*anchor alone falls into the noise on the most clamped images.
 *        clamp_mu/sigma is the latent background normal behind the clamped counts (zeros left-censored).
 *          C = CD99 positive      H = INS or GCG or SST positive (hormone union)
 *        IAPP (failed QC) and HLA-A (disease-associated) are never used. The same ops are saved as per-image
 *        pixel classifiers for review in the GUI.
 *     3. geometry (JTS), inside Tissue minus EXCLUDE_CLASSES (pLN, Artifact), clipped 5 um inside the image:
 *          U = H + (C - H); connected parts (sharing an edge = one object); holes filled; < MIN_AREA_UM2
 *          dropped. CD99-only pixels are NOT opened: CD99-only cell-sized blobs are endocrine candidates.
 *        The legacy size-adaptive erode/dilate is NOT applied: it deletes single endocrine cells.
 *        Optional, off by default (ISLETS.md; for calibration sweeps):
 *          NECK_SPLIT    cut thin CD99-only connections (narrower than 2 x NECK_R_UM) that join >= 2
 *                        separate hormone cores (hormone parts >= CORE_MIN_AREA_UM2) inside one object;
 *          SINGLE_COPOS  drop EndocrineSingles that are positive ONLY in one heavily clamped hormone channel
 *                        (level-0 clamp fraction >= CLAMP_HEAVY_FRAC), i.e. without CD99 or a second hormone;
 *          OPEN_H_R_UM   open the hormone pieces by this radius.
 *     4. size class by MEAN CALIPER DIAMETER = (Max diameter + Min diameter) / 2, QuPath's shape features
 *        (Drotar): 'Max diameter' = diameter of the minimum bounding circle, 'Min diameter' = minimum caliper
 *        width. Islet >= 35 um;  EndocrineCluster > 10 and < 35 um;  EndocrineSingle <= 10 um.
 *     5. per-object measurements: Mean caliper diameter, Equivalent diameter, Circularity, Solidity, Area,
 *        Endocrine route (1 CD99 only, 2 hormone only, 3 both), CD99+/INS+/GCG+/SST+/Hormone+ area, CD99-only
 *        fraction, Hormone cores + Multi-core flag (>= 2 cores), Neck split, Positive channels, Clamped-only
 *        flag, Bridge flag (area >= BRIDGE_AREA_UM2 and CD99-only fraction >= BRIDGE_CD99ONLY_FRAC), Tissue
 *        border flag (within BORDER_UM of the analysis-area edge; kept, as Heide/Bodenmiller), Distance to
 *        tissue border.
 *     6. image sanity: endocrine area / tissue area, flagged above MAX_ENDOCRINE_FRACTION (ADIIDAS).
 *
 * Previous Islet / EndocrineCluster / EndocrineSingle annotations (and their now-orphaned IsletExpanded
 * bands) are removed only after the new objects are built. Their child detections move up; any detection
 * that ends up directly under the root (e.g. from a band, which lives at the root) is re-parented to the
 * Tissue (else pLN) annotation containing its nucleus centroid. Objects are named and numbered by
 * 4_islet_expander_v3.groovy; script 8 parents cells. Nothing is saved here: the caller (run_v3.sh footer,
 * or the user in the GUI) decides.
 *
 * Parameters: edit the block below, pass script args KEY=VALUE (headless: -a TAU_H=4 -a FLOOR_K=4), or let
 * run_v3 rewrite a declaration. Calibration sweeps: the *_GRID parameters take comma lists; every combination
 * is built from ONE pixel pass, written next to WRITE_GEOJSON, and the hierarchy is left untouched.
 * SIGMA_UM or RESOLUTION_UM changes need their own anchors CSV (islet_anchors.py --sigma-um/--resolution-um)
 * and one run each.
 */

import qupath.lib.analysis.features.ObjectMeasurements
import qupath.lib.analysis.images.ContourTracing
import qupath.lib.images.ImageData
import qupath.lib.images.servers.ColorTransforms
import qupath.lib.images.servers.PixelCalibration
import qupath.lib.images.servers.PixelType
import qupath.lib.io.GsonTools
import qupath.lib.io.PathIO
import qupath.lib.objects.PathObjects
import qupath.lib.objects.PathObjectTools
import qupath.lib.objects.classes.PathClass
import qupath.lib.objects.hierarchy.PathObjectHierarchy
import qupath.lib.regions.RegionRequest
import qupath.lib.roi.GeometryTools
import qupath.lib.roi.RoiTools
import qupath.lib.roi.ROIs
import qupath.opencv.ml.pixel.PixelClassifierTools
import qupath.opencv.ml.pixel.PixelClassifiers
import qupath.opencv.ops.ImageOps
import com.google.gson.JsonParser
import org.locationtech.jts.algorithm.MinimumBoundingCircle
import org.locationtech.jts.algorithm.MinimumDiameter
import org.locationtech.jts.geom.Coordinate
import org.locationtech.jts.geom.Geometry
import org.locationtech.jts.geom.Location
import org.locationtech.jts.geom.Polygon
import org.locationtech.jts.geom.prep.PreparedGeometryFactory
import org.locationtech.jts.geom.util.PolygonExtracter
import org.locationtech.jts.index.strtree.STRtree
import org.locationtech.jts.operation.distance.IndexedFacetDistance
import org.locationtech.jts.operation.union.CascadedPolygonUnion
import org.locationtech.jts.operation.union.UnaryUnionOp
import java.nio.file.Paths
import java.security.MessageDigest
import java.util.concurrent.ConcurrentHashMap
import java.util.concurrent.ConcurrentLinkedQueue

evaluate(new File(buildFilePath(PROJECT_BASE_DIR, 'scripts', 'lib', 'common.groovy')))

// ==================== PARAMETERS ====================
String ISLET_METHOD            = 'drotar_cd99'   // 'drotar_cd99' | 'Islets2' (legacy baseline)
double TAU_CD99                = 3.0     // CD99 fold threshold -- provisional until the benchmark (ISLETS.md)
double TAU_H                   = 5.0     // hormone fold threshold (INS, GCG, SST each) -- provisional
double FLOOR_K                 = 4.0     // clamp floor: positive needs >= clamp_mu + FLOOR_K * clamp_sigma too (0 = off)
double SIGMA_UM                = 1.5     // Gaussian sigma before thresholding, um (the anchors must match)
double MIN_AREA_UM2            = 20.0    // smaller connected parts are dropped (a 5 um disc)
double RESOLUTION_UM           = 1.0156  // classifier resolution (pyramid level 1), as Islets2 (anchors must match)
boolean NECK_SPLIT             = false   // optional: cut thin CD99-only necks joining >= 2 hormone cores
double NECK_R_UM               = 4.0     // ... a neck is narrower than 2 x NECK_R_UM
double CORE_MIN_AREA_UM2       = 200.0   // a hormone core: one connected hormone+ part at least this large
boolean SINGLE_COPOS           = false   // optional: singles positive only in one heavily clamped hormone are dropped
double CLAMP_HEAVY_FRAC        = 0.6     // ... "heavily clamped": level-0 clamp fraction >= this
double OPEN_H_R_UM             = 0.0     // optional opening of hormone pixels (0 = off)
String ANCHORS_CSV             = '/home/smith6jt/IO60panc2nd/results/v3/islet_anchors.csv'
String ANCHOR_CHECK            = 'error' // annotation hash / pixel-space mismatch with the anchors: 'error' | 'warn'
String OUT_CLASSIFIER_DIR      = '/home/smith6jt/IO60panc2nd/classifiers/generated'
boolean SAVE_CLASSIFIERS       = true
double BRIDGE_AREA_UM2         = 5000.0  // bridge flag: large ...
double BRIDGE_CD99ONLY_FRAC    = 0.5     // ... and mostly CD99-only
double BORDER_UM               = 10.0
double EDGE_MARGIN_UM          = 5.0     // objects are clipped 5 um inside the image edge (as legacy)
String EXCLUDE_CLASSES         = 'pLN,Artifact'   // area classes cut out of the analysis area (lymph node,
                                                  // folds/debris), whether or not they overlap Tissue
double MAX_ENDOCRINE_FRACTION  = 0.10    // ADIIDAS-style sanity flag
double ISLET_MIN_DIAM_UM       = 35.0
double SINGLE_MAX_DIAM_UM      = 10.0
int TILE_PX                    = 1024    // classifier tile at RESOLUTION_UM (fewer seams than 512)
int THREADS                    = 0       // tile workers; 0 = the ForkJoin parallelism run_v3 sets, else QuPath's
String WRITE_GEOJSON           = ''      // tests: write objects + <name>.summary.json here
String TAU_CD99_GRID           = ''      // sweeps, e.g. '2,3,5' -> hierarchy untouched
String TAU_H_GRID              = ''
String FLOOR_K_GRID            = ''      // e.g. '0,3,4'
String NECK_SPLIT_GRID         = ''      // e.g. 'false,true'
String SINGLE_COPOS_GRID       = ''      // e.g. 'false,true'
String OPEN_H_GRID             = ''
// ====================================================

def P = new LinkedHashMap([
        ISLET_METHOD: ISLET_METHOD, TAU_CD99: TAU_CD99, TAU_H: TAU_H, FLOOR_K: FLOOR_K, SIGMA_UM: SIGMA_UM,
        MIN_AREA_UM2: MIN_AREA_UM2, RESOLUTION_UM: RESOLUTION_UM, NECK_SPLIT: NECK_SPLIT, NECK_R_UM: NECK_R_UM,
        CORE_MIN_AREA_UM2: CORE_MIN_AREA_UM2, SINGLE_COPOS: SINGLE_COPOS, CLAMP_HEAVY_FRAC: CLAMP_HEAVY_FRAC,
        OPEN_H_R_UM: OPEN_H_R_UM, ANCHORS_CSV: ANCHORS_CSV, ANCHOR_CHECK: ANCHOR_CHECK,
        OUT_CLASSIFIER_DIR: OUT_CLASSIFIER_DIR, SAVE_CLASSIFIERS: SAVE_CLASSIFIERS, BRIDGE_AREA_UM2: BRIDGE_AREA_UM2,
        BRIDGE_CD99ONLY_FRAC: BRIDGE_CD99ONLY_FRAC, BORDER_UM: BORDER_UM, EDGE_MARGIN_UM: EDGE_MARGIN_UM,
        EXCLUDE_CLASSES: EXCLUDE_CLASSES,
        MAX_ENDOCRINE_FRACTION: MAX_ENDOCRINE_FRACTION, ISLET_MIN_DIAM_UM: ISLET_MIN_DIAM_UM,
        SINGLE_MAX_DIAM_UM: SINGLE_MAX_DIAM_UM, TILE_PX: TILE_PX, THREADS: THREADS, WRITE_GEOJSON: WRITE_GEOJSON,
        TAU_CD99_GRID: TAU_CD99_GRID, TAU_H_GRID: TAU_H_GRID, FLOOR_K_GRID: FLOOR_K_GRID,
        NECK_SPLIT_GRID: NECK_SPLIT_GRID, SINGLE_COPOS_GRID: SINGLE_COPOS_GRID, OPEN_H_GRID: OPEN_H_GRID])
// Script args KEY=VALUE override the block above; an unknown key is an error (catches typos).
def scriptArgs = (binding.hasVariable('args') && args != null) ? (args as List).collect { it.toString() } : []
scriptArgs.each { String a ->
    int eq = a.indexOf('=')
    if (eq <= 0) throw new IllegalArgumentException("script argument '${a}' is not KEY=VALUE")
    String k = a.substring(0, eq).trim(), v = a.substring(eq + 1).trim()
    if (!P.containsKey(k)) throw new IllegalArgumentException("unknown parameter '${k}'; known: ${P.keySet()}")
    def old = P[k]
    if (old instanceof Boolean && !(v.toLowerCase() in ['true', 'false']))
        throw new IllegalArgumentException("parameter ${k} takes true or false, not '${v}'")
    P[k] = (old instanceof Double) ? Double.parseDouble(v) : (old instanceof Integer) ? Integer.parseInt(v) :
            (old instanceof Boolean) ? Boolean.parseBoolean(v) : v
}
println "[01v3] Parameters: ${P}"
if (!(P.ISLET_METHOD in ['drotar_cd99', 'Islets2']))
    throw new IllegalArgumentException("ISLET_METHOD must be 'drotar_cd99' or 'Islets2', not '${P.ISLET_METHOD}'")
if (!(P.ANCHOR_CHECK in ['error', 'warn']))
    throw new IllegalArgumentException("ANCHOR_CHECK must be 'error' or 'warn', not '${P.ANCHOR_CHECK}'")

setImageType('FLUORESCENCE')
setPixelSizeMicrons(0.5078, 0.5078)

// ---- Edge/inner analysis area (um) ----
double edgeMarginUm = P.EDGE_MARGIN_UM as double

// ---- Size-adaptive erosion parameters (legacy Islets2 path only) ----
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

// Endocrine classes this script owns. A re-run replaces all of them (and the bands derived from them).
def ENDOCRINE_CLASSES = ['Islet', 'EndocrineCluster', 'EndocrineSingle']

if (P.ISLET_METHOD == 'Islets2') {
    // =============================== LEGACY PATH (unchanged) ===============================
    // Run the classifier through a private Bio-Formats server instance. In QuPath <= 0.7.0 an
    // interrupted tile read (live classifier/threshold preview, viewer tile cancellation) closes the
    // FileChannel of a pooled Bio-Formats reader, and that reader is reused, so every later read on the
    // viewer's server fails. A separate instance has its own reader pool and pending-tile map.
    // Re-run safety: createAnnotationsFromPixelClassifier (without DELETE_EXISTING) adds new islets next to
    // any existing ones, and the erode/dilate step below would then erode the old islets a second time.
    // Remove previous Islet annotations only; keepChildren=true moves any parented detections up to Tissue.
    // (DELETE_EXISTING is not an option: it clears ALL children of Tissue, detections included.)
    // v3: drotar_cd99 objects (clusters/singles) are removed too, so the two methods never coexist.
    def oldIslets = getAnnotationObjects().findAll { it.getPathClass()?.toString() in ENDOCRINE_CLASSES }
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
} else {
    // ============================= drotar_cd99 PATH (v3 default) =============================
    long tStart = System.currentTimeMillis()
    def imageData = getCurrentImageData()
    def hierarchy = imageData.getHierarchy()
    def entry = getProjectEntry()
    String imageName = entry?.getImageName() ?: server.getMetadata().getName()
    String tag = imageName.replaceAll(/[^A-Za-z0-9._-]+/, '_')
    String UM = cal.getPixelWidthUnit()                  // 'µm': measurement names match QuPath's own
    double pxArea = sx * sy                              // um^2 per full-resolution pixel^2
    def fmt = { double v -> (v == Math.rint(v)) ? String.valueOf((long) v) : String.valueOf(v) }
    def parseGrid = { String s, double dflt ->
        (s == null || s.trim().isEmpty()) ? [dflt] : s.split(',').collect { it.trim() }.findAll { it }.collect { it as double }
    }
    def parseBoolGrid = { String s, boolean dflt ->
        if (s == null || s.trim().isEmpty()) return [dflt]
        s.split(',').collect { it.trim().toLowerCase() }.findAll { it }.collect { v ->
            if (!(v in ['true', 'false'])) throw new IllegalArgumentException("boolean grid value '${v}' is not true/false")
            v == 'true'
        }.unique()
    }
    def tauCs = parseGrid(P.TAU_CD99_GRID as String, P.TAU_CD99 as double)
    def tauHs = parseGrid(P.TAU_H_GRID as String, P.TAU_H as double)
    def floorKs = parseGrid(P.FLOOR_K_GRID as String, P.FLOOR_K as double)
    def openHs = parseGrid(P.OPEN_H_GRID as String, P.OPEN_H_R_UM as double)
    def necks = parseBoolGrid(P.NECK_SPLIT_GRID as String, P.NECK_SPLIT as boolean)
    def coposs = parseBoolGrid(P.SINGLE_COPOS_GRID as String, P.SINGLE_COPOS as boolean)
    def gridKeys = ['TAU_CD99_GRID', 'TAU_H_GRID', 'FLOOR_K_GRID', 'NECK_SPLIT_GRID', 'SINGLE_COPOS_GRID', 'OPEN_H_GRID']
    boolean sweep = gridKeys.any { P[it] && !(P[it] as String).trim().isEmpty() }
    if (sweep && !(P.WRITE_GEOJSON as String))
        throw new IllegalArgumentException('a sweep writes its results next to WRITE_GEOJSON; set it')
    def checkIssues = []   // anchor/annotation mismatches (ANCHOR_CHECK decides: stop or warn)

    def safeDiffTop = { a, b ->
        try { return a.difference(b) } catch (Exception e) { return a.buffer(0).difference(b.buffer(0)) }
    }

    // ---- anchors + clamp models: one row per (image, channel) in islet_anchors.csv ----
    final char DQ = (char) 34
    final char COMMA = (char) 44
    def parseCsvLine = { String line ->   // RFC 4180 fields (quotes, doubled quotes); no embedded newlines
        def out = []
        def sb = new StringBuilder()
        boolean inQ = false
        for (int i = 0; i < line.length(); i++) {
            char c = line.charAt(i)
            if (inQ) {
                if (c == DQ) {
                    if (i + 1 < line.length() && line.charAt(i + 1) == DQ) { sb.append(DQ); i++ } else inQ = false
                } else sb.append(c)
            } else if (c == DQ) inQ = true
            else if (c == COMMA) { out << sb.toString(); sb.setLength(0) }
            else sb.append(c)
        }
        out << sb.toString()
        out
    }
    def anchorsFile = new File(P.ANCHORS_CSV as String)
    if (!anchorsFile.exists()) throw new IOException("anchors CSV not found: ${anchorsFile} (run scripts/v3/islet_anchors.py)")
    def lines = anchorsFile.readLines('UTF-8').findAll { it.trim() }
    def header = parseCsvLine(lines[0])
    def col = { String n -> int i = header.indexOf(n); if (i < 0) throw new IOException("anchors CSV lacks column ${n}"); i }
    def optCol = { String n -> header.indexOf(n) }
    int iImg = col('image'), iTag = col('tag'), iCh = col('channel'), iAnc = col('anchor'), iVer = col('method_version')
    int iMu = optCol('clamp_mu'), iSd = optCol('clamp_sigma'), iFr = optCol('clamp_fraction')
    int iHash = optCol('annotation_hash'), iHcls = optCol('annotation_classes'), iSig = optCol('sigma_um'), iLvl = optCol('level_um')
    def rowsCsv = lines.drop(1).collect { parseCsvLine(it) }
    def mine = rowsCsv.findAll { it[iImg] == imageName }
    if (mine.isEmpty()) mine = rowsCsv.findAll { it[iTag] == tag }
    def CH = ['CD99', 'INS', 'GCG', 'SST']
    def HORMONES = ['INS', 'GCG', 'SST']
    def num = { List r, int i -> (i < 0 || r == null || i >= r.size() || !r[i] || r[i] in ['nan', 'NaN']) ? Double.NaN : (r[i] as double) }
    def anchors = [:], clampMu = [:], clampSd = [:], clampFrac = [:]
    def anchorVersion = mine ? mine[0][iVer] : null
    CH.each { ch ->
        def r = mine.find { it[iCh] == ch }
        double a = num(r, iAnc)
        // explicit NaN test: in Groovy `!(NaN > 0)` is false, so a missing row or a NaN anchor would slip through
        if (r == null || Double.isNaN(a) || Double.isInfinite(a) || a <= 0)
            throw new IllegalStateException("no valid anchor for ${imageName} ${ch} in ${anchorsFile}")
        anchors[ch] = a
        clampMu[ch] = num(r, iMu); clampSd[ch] = num(r, iSd); clampFrac[ch] = num(r, iFr)
    }
    def chNames = getCurrentServer().getMetadata().getChannels().collect { it.getName() }
    CH.each { ch ->
        int n = chNames.count { it == ch }
        if (n != 1) throw new IllegalStateException("channel ${ch} found ${n} times; need exactly one")
    }
    println "[01v3] ${imageName}: anchors ${anchors} clamp mu ${clampMu} sigma ${clampSd} fraction ${clampFrac} (islet_anchors ${anchorVersion})"
    // the anchors must come from this pixel space (finding: estimated at 2 um without blur, applied at 1 um blurred)
    if (iMu < 0 || iSig < 0 || iLvl < 0 || iHash < 0 || iHcls < 0) {
        checkIssues << "anchors CSV ${anchorsFile} is islet_anchors method ${anchorVersion} (no clamp model / pixel space / " +
                "annotation hash columns): recompute with scripts/v3/islet_anchors.py (method >= 2.0)"
    } else {
        double aSig = num(mine[0], iSig), aLvl = num(mine[0], iLvl)
        if (!(Math.abs(aSig - (P.SIGMA_UM as double)) < 1e-6))
            checkIssues << "anchors were estimated with sigma ${aSig} um, this run uses SIGMA_UM=${P.SIGMA_UM} " +
                    "(islet_anchors.py --sigma-um ${P.SIGMA_UM} --out <csv>, then ANCHORS_CSV=<csv>)"
        if (!(Math.abs(aLvl - (P.RESOLUTION_UM as double)) < 0.002))
            checkIssues << "anchors were estimated at ${aLvl} um/px, this run classifies at RESOLUTION_UM=${P.RESOLUTION_UM}"
        CH.each { ch -> if (Double.isNaN(clampMu[ch]) || Double.isNaN(clampSd[ch]))
            println "[01v3] WARN: no clamp model for ${ch} (fit not identifiable); its floor is off, tau * anchor alone applies" }
    }

    // ---- analysis area: Tissue minus EXCLUDE_CLASSES (pLN, Artifact), clipped EDGE_MARGIN_UM inside the image ----
    def tissueAnns = getAnnotationObjects().findAll { it.getPathClass()?.toString() == CLASS_TISSUE && it.getROI()?.isArea() }
    if (tissueAnns.isEmpty()) throw new IllegalStateException("no ${CLASS_TISSUE} annotation")
    def plane = tissueAnns[0].getROI().getImagePlane()
    def tissueGeom = UnaryUnionOp.union(tissueAnns.collect { it.getROI().getGeometry() })
    def excludeNames = (P.EXCLUDE_CLASSES as String).split(',').collect { it.trim() }.findAll { it }
    def excluded = getAnnotationObjects().findAll { it.getPathClass()?.toString() in excludeNames && it.getROI()?.isArea() }
    if (!excluded.isEmpty()) {
        def exGeom = UnaryUnionOp.union(excluded.collect { it.getROI().getGeometry() })
        double before = tissueGeom.getArea()
        tissueGeom = safeDiffTop(tissueGeom, exGeom)
        println String.format("[01v3] excluded %s: %.3f mm^2 of Tissue", excluded.countBy { it.getPathClass().toString() },
                (before - tissueGeom.getArea()) * pxArea / 1e6)
    }
    def clipGeom = toPolygonal(tissueGeom.intersection(innerGeomForPlane(plane)))
    if (clipGeom == null || clipGeom.isEmpty()) throw new IllegalStateException('Tissue lies outside the analysis area')
    double tissueAreaUm2 = clipGeom.getArea() * pxArea
    def clipROI = GeometryTools.geometryToROI(clipGeom, plane)
    println String.format("[01v3] analysis area %.2f mm^2 (Tissue clipped %.0f um inside the image edge)",
            tissueAreaUm2 / 1e6, edgeMarginUm)

    // ---- annotation hash: the anchors' analysis mask must be this analysis area ----
    // Same canonical text as islet_anchors.annotation_signature: the annotations are serialised with QuPath's own
    // GeoJSON exporter (the one export_annotations.groovy used), parsed back, every vertex printed with 2 decimals.
    def annotationHash = { List<String> classes ->
        def anns = getAnnotationObjects().findAll { it.getPathClass()?.toString() in classes && it.getROI() != null }
        def sigLines = []
        if (!anns.isEmpty()) {
            def bos = new ByteArrayOutputStream()
            PathIO.exportObjectsAsGeoJSON(bos, anns, PathIO.GeoJsonExportOptions.FEATURE_COLLECTION)
            def fc = JsonParser.parseString(bos.toString('UTF-8')).getAsJsonObject()
            fc.getAsJsonArray('features').each { f ->
                def fo = f.getAsJsonObject()
                def cl = fo.getAsJsonObject('properties')?.getAsJsonObject('classification')
                String c = cl?.get('name')?.getAsString()
                if (!(c in classes)) return
                def g = fo.getAsJsonObject('geometry')
                String type = g?.get('type')?.getAsString()
                def polys = type == 'Polygon' ? [g.getAsJsonArray('coordinates')] :
                        type == 'MultiPolygon' ? g.getAsJsonArray('coordinates').collect { it.getAsJsonArray() } : []
                polys.each { poly ->
                    def rings = poly.collect { it.getAsJsonArray() }
                    if (rings.isEmpty() || rings[0].size() < 3) return
                    String body = rings.collect { ring -> ring.collect { pt -> def xy = pt.getAsJsonArray()
                        String.format(Locale.ROOT, '%.2f,%.2f', xy.get(0).getAsDouble(), xy.get(1).getAsDouble()) }.join(';') }.join('/')
                    sigLines << (c + '|' + body)
                }
            }
        }
        sigLines.sort()
        String text = 'classes=' + new ArrayList(classes).sort().join(',') + '\n' + sigLines.join('\n')
        MessageDigest.getInstance('SHA-256').digest(text.getBytes('UTF-8')).collect { String.format('%02x', it) }.join('').substring(0, 16)
    }
    String csvHash = (iHash >= 0) ? mine[0][iHash] : null
    def csvHashClasses = (iHcls >= 0 && mine[0][iHcls]) ? mine[0][iHcls].split(',').collect { it.trim() }.findAll { it } : []
    def myClasses = ([CLASS_TISSUE] + excludeNames).unique()
    String nowHash = csvHashClasses ? annotationHash(csvHashClasses) : null
    if (csvHash != null) {
        if (new HashSet(csvHashClasses) != new HashSet(myClasses))
            checkIssues << "anchors used annotation classes ${csvHashClasses}; this run analyses ${myClasses} (EXCLUDE_CLASSES)"
        if (nowHash != csvHash)
            checkIssues << "the ${csvHashClasses.join('/')} annotations changed since the anchors were computed (hash ${csvHash} in " +
                    "${anchorsFile.getName()}, ${nowHash} now): re-export the annotations (v3/export_annotations.sh) and rerun " +
                    "islet_anchors.py --images '${imageName}' (or --all --skip-existing)"
    }
    if (!checkIssues.isEmpty()) {
        String msg = "[01v3] ANCHOR CHECK FAILED for ${imageName}:\n  - " + checkIssues.join('\n  - ')
        if (P.ANCHOR_CHECK == 'error') throw new IllegalStateException(msg + "\n(ANCHOR_CHECK=warn runs anyway)")
        println "!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!\n${msg}\n!!!!!!!!!!!!!!!! (ANCHOR_CHECK=warn: continuing) !!!!!!!!!!!!!!!!"
    } else {
        println "[01v3] anchor check OK: annotation hash ${nowHash}, sigma ${P.SIGMA_UM} um, ${P.RESOLUTION_UM} um/px"
    }

    // ---- per-channel thresholds (fold units, x / anchor): max(tau, floor / anchor) ----
    def floorCounts = { String ch, double k ->
        (k > 0 && !Double.isNaN(clampMu[ch] as double) && !Double.isNaN(clampSd[ch] as double)) ?
                (clampMu[ch] as double) + k * (clampSd[ch] as double) : Double.NEGATIVE_INFINITY
    }
    def thrFold = { String ch, double tau, double k -> Math.max(tau, floorCounts(ch, k) / (anchors[ch] as double)) }

    // ---- pixel ops: fold (divide by anchor) -> Gaussian(SIGMA_UM) -> [threshold] ----
    double res = P.RESOLUTION_UM as double
    double sigmaPx = (P.SIGMA_UM as double) / res
    def calRes = new PixelCalibration.Builder().pixelSizeMicrons(res, res).build()
    def foldOp = { List<String> chans ->
        ImageOps.buildImageDataOp(chans.collect { ColorTransforms.createChannelExtractor(it) })
                .appendOps(ImageOps.Core.multiply(chans.collect { 1.0d / (anchors[it] as double) } as double[]),
                           ImageOps.Filters.gaussianBlur(sigmaPx))
    }
    // H as one channel with threshold 1: max over hormones of fold_c / t_c (each hormone has its own threshold)
    def hRatioOp = { double tI, double tG, double tS ->
        ImageOps.Core.sequential(ImageOps.Channels.extract(1, 2, 3), ImageOps.Core.multiply(1.0d / tI, 1.0d / tG, 1.0d / tS),
                ImageOps.Channels.maximum())
    }

    // Traced data server: float folds, one output channel per (mask, threshold). Thresholds are applied by the
    // contour tracer (value >= threshold, as specified; the saved GUI classifiers' threshold op is value > threshold,
    // identical for continuous values), so a sweep costs no extra pixel reads. A PROBABILITY-type server is NOT used
    // on purpose: ContourTracing collapses multi-channel probabilities to an argmax label, which would make the masks
    // mutually exclusive, whereas C and H must overlap. ensureType(FLOAT32): the op server otherwise declares the
    // input's UINT16 and would truncate fractional folds. Identical (mask, threshold) requests share one channel.
    def layout = new LinkedHashMap()   // key -> [ch: output channel, thr: threshold, op: ImageOp, mask: name]
    def request = { String mask, String key, op, double thr ->
        if (!layout.containsKey(key)) layout[key] = [ch: layout.size(), thr: thr, op: op, mask: mask]
        key
    }
    def combos = []   // one per (tC, tH, k): the channel keys of its masks
    tauCs.each { double tC -> tauHs.each { double tH -> floorKs.each { double k ->
        double fC = thrFold('CD99', tC, k), fI = thrFold('INS', tH, k), fG = thrFold('GCG', tH, k), fS = thrFold('SST', tH, k)
        def keys = [C  : request('C', "C>${fC}", ImageOps.Channels.extract(0), fC),
                    INS: request('INS', "INS>${fI}", ImageOps.Channels.extract(1), fI),
                    GCG: request('GCG', "GCG>${fG}", ImageOps.Channels.extract(2), fG),
                    SST: request('SST', "SST>${fS}", ImageOps.Channels.extract(3), fS),
                    H  : request('H', "H>${fI},${fG},${fS}", hRatioOp(fI, fG, fS), 1.0d)]
        def thrCounts = [CD99: fC * (anchors.CD99 as double), INS: fI * (anchors.INS as double),
                         GCG: fG * (anchors.GCG as double), SST: fS * (anchors.SST as double)]
        def floorBinds = CH.findAll { ch -> floorCounts(ch, k) > (ch == 'CD99' ? tC : tH) * (anchors[ch] as double) }
        combos << [tC: tC, tH: tH, k: k, keys: keys, thr_fold: [CD99: fC, INS: fI, GCG: fG, SST: fS],
                   thr_counts: thrCounts, floor_binds: floorBinds]
        println String.format(Locale.ROOT, "[01v3] tC=%s tH=%s k=%s: thresholds (counts) CD99 %.1f INS %.1f GCG %.1f SST %.1f; floor binds for %s",
                fmt(tC), fmt(tH), fmt(k), thrCounts.CD99, thrCounts.INS, thrCounts.GCG, thrCounts.SST, floorBinds ?: 'none')
    } } }
    def dataOp = foldOp(CH).appendOps(ImageOps.Core.splitMerge(layout.values().collect { it.op }), ImageOps.Core.ensureType(PixelType.FLOAT32))
    def thresholds = layout.values().collect { ContourTracing.ChannelThreshold.createAbove(it.ch as int, it.thr as double) } as ContourTracing.ChannelThreshold[]

    int nThreads = (P.THREADS as int) > 0 ? (P.THREADS as int) :
            (System.getProperty('java.util.concurrent.ForkJoinPool.common.parallelism')?.isInteger() ?
                    System.getProperty('java.util.concurrent.ForkJoinPool.common.parallelism') as int : 0)
    if (nThreads > 0) {   // each worker holds a padded all-channel tile, so this bounds native memory too
        qupath.lib.common.ThreadTools.setParallelism(nThreads)
        println "[01v3] tile workers: ${nThreads}"
    }
    long tTrace0 = System.currentTimeMillis()
    Map<Integer, Geometry> traced
    def privateServer = imageData.getServerBuilder().build()   // own reader pool (see legacy note above)
    try {
        privateServer.setMetadata(imageData.getServer().getMetadata())
        def privateData = new ImageData(privateServer, new PathObjectHierarchy(), imageData.getImageType())
        def ds = ImageOps.buildServer(privateData, dataOp, calRes, P.TILE_PX as int, P.TILE_PX as int)
        def rq = RegionRequest.createInstance(ds.getPath(), ds.getDownsampleForResolution(0), clipROI)
        traced = ContourTracing.traceGeometries(ds, rq, clipGeom, thresholds)
    } finally {
        privateServer.close()
    }
    double traceSec = (System.currentTimeMillis() - tTrace0) / 1000.0
    def maskGeom = { String key ->
        def g = traced.get(layout[key].ch as int)
        (g == null || g.isEmpty()) ? null : g
    }
    println String.format("[01v3] traced %d masks in %.1f s: %s", layout.size(), traceSec,
            layout.collect { k, L -> def g = traced.get(L.ch as int); "${k}: " + String.format('%.3f mm^2', (g == null ? 0 : g.getArea()) * pxArea / 1e6) }.join(', '))

    // ---- JTS helpers ----
    def gf = GeometryTools.getDefaultFactory()
    def polysOf = { Geometry g -> g == null ? [] : new ArrayList(PolygonExtracter.getPolygons(g)).findAll { !it.isEmpty() } }
    def treeOf = { List polys ->
        def t = new STRtree()
        polys.each { t.insert(it.getEnvelopeInternal(), it) }
        t.build()
        t
    }
    // Area of g inside the union of a tree's polygons (the traced parts of one mask are disjoint)
    def areaIn = { Geometry g, STRtree tree ->
        double a = 0
        def prep = PreparedGeometryFactory.prepare(g)
        for (def p : tree.query(g.getEnvelopeInternal())) {
            if (!prep.intersects(p)) continue
            if (prep.covers(p)) { a += p.getArea(); continue }
            try { a += g.intersection(p).getArea() }
            catch (Exception e) { a += g.buffer(0).intersection(p.buffer(0)).getArea() }
        }
        a
    }
    // Two pieces belong to one object when they share an EDGE (or overlap); a single shared corner
    // (diagonal pixels) does not join them -- the tracer's own 4-connectivity.
    def shareEdge = { Geometry a, Geometry b ->
        try {
            def im = a.relate(b)
            return im.get(Location.BOUNDARY, Location.BOUNDARY) >= 1 || im.get(Location.INTERIOR, Location.INTERIOR) >= 0
        } catch (Exception e) {   // robustness fallback: treat any contact as a join
            return a.intersects(b)
        }
    }
    def safeUnion = { List gs ->
        try { return CascadedPolygonUnion.union(gs) }
        catch (Exception e) { return UnaryUnionOp.union(gs.collect { it.buffer(0) }) }
    }
    def safeDiff = { Geometry a, Geometry b ->
        try { return a.difference(b) }
        catch (Exception e) { return a.buffer(0).difference(b.buffer(0)) }
    }
    def openPoly = { Polygon p, double rPx ->   // morphological opening, never larger than p
        def e = p.buffer(-rPx)
        if (e.isEmpty()) return null
        try { return e.buffer(rPx).intersection(p) }
        catch (Exception ex) { return e.buffer(rPx).buffer(0).intersection(p.buffer(0)) }
    }
    def filledPolys = { Geometry g -> polysOf(g).collect { Polygon poly -> gf.createPolygon(poly.getExteriorRing().getCoordinates()) } }
    // connected components of pieces by shared edges (union-find over STRtree candidate pairs); returns index lists
    def components = { List pieces ->
        int n = pieces.size()
        def idxTree = new STRtree()
        pieces.eachWithIndex { g, i -> idxTree.insert(g.getEnvelopeInternal(), i) }
        idxTree.build()
        def pairs = new ConcurrentLinkedQueue()
        (0..<n).toList().parallelStream().forEach { int i ->
            def gi = pieces[i]
            for (def j : idxTree.query(gi.getEnvelopeInternal())) {
                int jj = j as int
                if (jj > i && shareEdge(gi, pieces[jj])) pairs.add([i, jj] as int[])
            }
        }
        int[] parent = new int[n]
        for (int i = 0; i < n; i++) parent[i] = i
        def find = { int i -> while (parent[i] != i) { parent[i] = parent[parent[i]]; i = parent[i] }; i }
        pairs.each { int[] pr -> int a = find(pr[0]), b = find(pr[1]); if (a != b) parent[a] = b }
        def groups = [:].withDefault { [] }
        for (int i = 0; i < n; i++) groups[find(i)] << i
        groups.values().toList()
    }
    def boundary = clipGeom.getBoundary()
    def borderDist = new IndexedFacetDistance(boundary)
    borderDist.distance(gf.createPoint(clipGeom.getInteriorPoint().getCoordinate()))   // build the index before parallel use
    def pcIslet = PathClass.fromString(CLASS_ISLET)
    def pcCluster = PathClass.fromString('EndocrineCluster')
    def pcSingle = PathClass.fromString('EndocrineSingle')
    double minPiecePx2 = (res / pxUm) * (res / pxUm)   // one classifier pixel, for "has CD99 / has hormone"
    double minAreaPx2 = (P.MIN_AREA_UM2 as double) / pxArea
    double coreMinPx2 = (P.CORE_MIN_AREA_UM2 as double) / pxArea

    // ---- base objects for one (tC, tH, k, opening) combination: U = H + (C - H), components, holes filled ----
    def buildBase = { Map combo, double openHR ->
        long t0 = System.currentTimeMillis()
        def cPolys = polysOf(maskGeom(combo.keys.C))
        def hPolys = polysOf(maskGeom(combo.keys.H))
        def hTree = treeOf(hPolys)
        // C_only = C - H, per C part, against the H parts that touch it
        def cOnly = new ConcurrentLinkedQueue()
        cPolys.parallelStream().forEach { p ->
            def hs = hTree.query(p.getEnvelopeInternal()).findAll { it.intersects(p) }
            def d = hs.isEmpty() ? p : safeDiff(p, hs.size() == 1 ? hs[0] : safeUnion(hs))
            polysOf(d).each { if (it.getArea() > 0) cOnly.add(it) }
        }
        def cOnlyList = new ArrayList(cOnly)
        // optional: open the hormone pieces (OPEN_H_R_UM > 0); C_only above still uses the unopened H, so speckle
        // cannot turn CD99 pixels into CD99-only ones. Composition and route below also use the unopened masks.
        def hPieces = hPolys
        if (openHR > 0) {
            double rh = openHR / pxUm
            def hq = new ConcurrentLinkedQueue()
            hPolys.parallelStream().forEach { Polygon hp -> polysOf(openPoly(hp, rh)).each { if (it.getArea() > 0) hq.add(it) } }
            hPieces = new ArrayList(hq)
        }
        def pieces = new ArrayList(hPieces)
        pieces.addAll(cOnlyList)
        def groups = components(pieces)
        // union each group, fill holes (Heide), keep polygons
        def filled = new ConcurrentLinkedQueue()
        groups.parallelStream().forEach { List idx ->
            def g = idx.size() == 1 ? pieces[idx[0]] : safeUnion(idx.collect { pieces[it] })
            filledPolys(g).each { filled.add(it) }
        }
        def objsGeom = new ArrayList(filled)
        // hole filling can swallow a separate part that sat inside a hole: drop parts covered by another
        def fTree = treeOf(objsGeom)
        def swallowed = Collections.newSetFromMap(new IdentityHashMap())
        objsGeom.each { g ->
            for (def o : fTree.query(g.getEnvelopeInternal())) {
                if (!o.is(g) && o.getArea() > g.getArea() && o.covers(g)) { swallowed.add(g); break }
            }
        }
        objsGeom = objsGeom.findAll { !swallowed.contains(it) && it.getArea() >= minAreaPx2 }
        def coreList = hPolys.findAll { it.getArea() >= coreMinPx2 }
        [objsGeom: objsGeom, cPolys: cPolys, hPolys: hPolys, hTree: hTree, cOnlyList: cOnlyList, cOnlyTree: treeOf(cOnlyList),
         coreTree: treeOf(coreList), coreList: coreList, n_pieces: pieces.size(), n_hormone_parts_after_opening: hPieces.size(),
         n_hole_swallowed: swallowed.size(), base_seconds: (System.currentTimeMillis() - t0) / 1000.0]
    }

    // cores (hormone parts >= CORE_MIN_AREA_UM2) whose interior point lies in g
    def coresIn = { Geometry g, STRtree coreTree ->
        def prep = PreparedGeometryFactory.prepare(g)
        coreTree.query(g.getEnvelopeInternal()).findAll { prep.covers(it.getInteriorPoint()) }
    }

    // ---- optional targeted neck split of ONE object: cut thin CD99-only necks between >= 2 hormone cores ----
    // thin = the parts of the object's CD99-only pixels removed by an opening of radius NECK_R_UM. If removing
    // them leaves the cores in >= 2 separate components, the object is split there: each thin part touching
    // exactly one component is given back to it, thin parts touching >= 2 components (the necks) are dropped.
    // Otherwise the object is unchanged. CD99-only regions thicker than the neck width are never cut.
    def neckSplit = { Geometry g, Map base, List cores ->
        double rPx = (P.NECK_R_UM as double) / pxUm
        def prep = PreparedGeometryFactory.prepare(g)
        def thin = []
        base.cOnlyTree.query(g.getEnvelopeInternal()).each { Polygon q ->
            if (!prep.intersects(q)) return
            def o = openPoly(q, rPx)
            def t = (o == null || o.isEmpty()) ? q : safeDiff(q, o)
            polysOf(t).each { if (it.getArea() > 0) thin << it }
        }
        if (thin.isEmpty()) return null
        def rest = safeDiff(g, safeUnion(thin))
        def comps = polysOf(rest)
        if (comps.size() < 2) return null
        def compOfCore = cores.collect { core -> def ip = core.getInteriorPoint(); comps.findIndexOf { it.covers(ip) } }
        if (compOfCore.findAll { it >= 0 }.unique().size() < 2) return null
        def attached = comps.collect { [it] }
        int nNecks = 0
        thin.each { tp ->
            def touching = (0..<comps.size()).findAll { shareEdge(tp, comps[it]) }
            if (touching.size() == 1) attached[touching[0]] << tp
            else if (touching.size() >= 2) nNecks++
        }
        def out = []
        attached.each { parts -> filledPolys(parts.size() == 1 ? parts[0] : safeUnion(parts)).each { if (it.getArea() >= minAreaPx2) out << it } }
        out.size() >= 2 ? [geoms: out, n_necks: nNecks] : null
    }

    // ---- objects for one full combination (base + options): measurements, classes, optional filters ----
    def finishObjects = { Map combo, Map base, double openHR, boolean neck, boolean copos ->
        long t0 = System.currentTimeMillis()
        // cores per base object; optional neck split
        def geoms = []          // [geom, split flag]
        int nSplitObjects = 0, nNecksCut = 0
        def splitQ = new ConcurrentLinkedQueue()
        base.objsGeom.parallelStream().forEach { Geometry g ->
            def cores = coresIn(g, base.coreTree)
            def sp = (neck && cores.size() >= 2) ? neckSplit(g, base, cores) : null
            splitQ.add(sp == null ? [[g, 0]] : sp.geoms.collect { [it, 1] } + [[null, sp.n_necks]])
        }
        splitQ.each { List items ->
            items.each { it2 -> if (it2[0] == null) { nSplitObjects++; nNecksCut += it2[1] as int } else geoms << it2 }
        }
        geoms.sort { a, b ->   // deterministic object order (the parallel streams above are not)
            def ea = a[0].getEnvelopeInternal(), eb = b[0].getEnvelopeInternal()
            ea.getMinX() <=> eb.getMinX() ?: ea.getMinY() <=> eb.getMinY() ?: a[0].getArea() <=> b[0].getArea()
        }
        def insTree = treeOf(polysOf(maskGeom(combo.keys.INS))), gcgTree = treeOf(polysOf(maskGeom(combo.keys.GCG)))
        def sstTree = treeOf(polysOf(maskGeom(combo.keys.SST))), cTree = treeOf(base.cPolys)
        def objs = geoms.collect { PathObjects.createAnnotationObject(GeometryTools.geometryToROI(it[0], plane)) }
        def splitFlag = new IdentityHashMap()
        objs.eachWithIndex { o, i -> splitFlag[o] = geoms[i][1] }
        ObjectMeasurements.addShapeMeasurements(objs, cal,
                ObjectMeasurements.ShapeFeatures.AREA, ObjectMeasurements.ShapeFeatures.LENGTH,
                ObjectMeasurements.ShapeFeatures.CIRCULARITY, ObjectMeasurements.ShapeFeatures.SOLIDITY,
                ObjectMeasurements.ShapeFeatures.MAX_DIAMETER, ObjectMeasurements.ShapeFeatures.MIN_DIAMETER)
        String kMax = "Max diameter ${UM}", kMin = "Min diameter ${UM}"
        objs.parallelStream().forEach { o ->
            def g = o.getROI().getGeometry()
            double areaPx = g.getArea()
            double aC = areaIn(g, cTree), aH = areaIn(g, base.hTree), aCo = areaIn(g, base.cOnlyTree)
            double aI = areaIn(g, insTree), aG = areaIn(g, gcgTree), aS = areaIn(g, sstTree)
            def ml = o.getMeasurementList()
            double dMax = ml.get(kMax), dMin = ml.get(kMin)
            if (Double.isNaN(dMax) || Double.isNaN(dMin)) {
                // QuPath's own definitions (ObjectMeasurements, 0.7): 'Max diameter' = 2 x the radius of the MINIMUM
                // BOUNDING CIRCLE (not the longest chord), 'Min diameter' = the minimum caliper width
                dMax = 2.0 * new MinimumBoundingCircle(g).getRadius() * pxUm
                dMin = new MinimumDiameter(g).getLength() * pxUm
            }
            double meanCal = (dMax + dMin) / 2.0
            double areaUm2 = areaPx * pxArea
            double cd99OnlyFrac = areaPx > 0 ? aCo / areaPx : 0
            boolean hasC = aC >= minPiecePx2, hasH = aH >= minPiecePx2
            int route = (hasC && hasH) ? 3 : (hasH ? 2 : 1)
            def pos = [CD99: aC >= minPiecePx2, INS: aI >= minPiecePx2, GCG: aG >= minPiecePx2, SST: aS >= minPiecePx2]
            int nPos = pos.count { k2, v -> v }
            def onlyCh = nPos == 1 ? pos.find { k2, v -> v }.key : null
            boolean clampedOnly = onlyCh != null && onlyCh != 'CD99' && !Double.isNaN(clampFrac[onlyCh] as double) &&
                    (clampFrac[onlyCh] as double) >= (P.CLAMP_HEAVY_FRAC as double)
            int nCores = coresIn(g, base.coreTree).size()
            double dBorder = borderDist.distance(g) * pxUm
            ml.put("Mean caliper diameter ${UM}", meanCal)
            ml.put("Equivalent diameter ${UM}", 2.0 * Math.sqrt(areaUm2 / Math.PI))
            ml.put('Endocrine route', route)
            ml.put("CD99+ area ${UM}^2", aC * pxArea)
            ml.put("INS+ area ${UM}^2", aI * pxArea)
            ml.put("GCG+ area ${UM}^2", aG * pxArea)
            ml.put("SST+ area ${UM}^2", aS * pxArea)
            ml.put("Hormone+ area ${UM}^2", aH * pxArea)
            ml.put('CD99-only fraction', cd99OnlyFrac)
            ml.put('Hormone cores', nCores)
            ml.put('Multi-core flag', nCores >= 2 ? 1 : 0)
            ml.put('Neck split', splitFlag[o] as int)
            ml.put('Positive channels', nPos)
            ml.put('Clamped-only flag', clampedOnly ? 1 : 0)
            ml.put('Bridge flag', (areaUm2 >= (P.BRIDGE_AREA_UM2 as double) && cd99OnlyFrac >= (P.BRIDGE_CD99ONLY_FRAC as double)) ? 1 : 0)
            ml.put('Tissue border flag', dBorder <= (P.BORDER_UM as double) ? 1 : 0)
            ml.put("Distance to tissue border ${UM}", dBorder)
            ml.close()
            o.setPathClass(meanCal >= (P.ISLET_MIN_DIAM_UM as double) ? pcIslet :
                    (meanCal > (P.SINGLE_MAX_DIAM_UM as double) ? pcCluster : pcSingle))
        }
        // ---- optional: singles positive only in one heavily clamped hormone channel are dropped ----
        def cls = { o -> o.getPathClass().toString() }
        def m = { o, String k -> o.getMeasurementList().get(k) }
        def dropped = copos ? objs.findAll { cls(it) == 'EndocrineSingle' && m(it, 'Clamped-only flag') == 1 } : []
        if (!dropped.isEmpty()) {
            def dset = Collections.newSetFromMap(new IdentityHashMap()); dset.addAll(dropped)
            objs = objs.findAll { !dset.contains(it) }
        }
        // ---- summary ----
        double endoArea = objs.sum(0.0d) { it.getROI().getGeometry().getArea() * pxArea } as double
        def byClassRoute = new TreeMap()
        ENDOCRINE_CLASSES.each { c -> byClassRoute[c] = new TreeMap([1: 0, 2: 0, 3: 0]) }
        objs.each { o -> byClassRoute[cls(o)][(int) m(o, 'Endocrine route')] += 1 }
        def singleOnly = { String ch -> objs.count { o -> cls(o) == 'EndocrineSingle' && m(o, 'Positive channels') == 1 &&
                m(o, "${ch}+ area ${UM}^2".toString()) > 0 } }
        double mm2 = tissueAreaUm2 / 1e6
        def summary = new LinkedHashMap([
                tau_cd99: combo.tC, tau_h: combo.tH, floor_k: combo.k, neck_split: neck, single_copos: copos,
                open_h_r_um: openHR, sigma_um: P.SIGMA_UM, resolution_um: res,
                thresholds_counts: combo.thr_counts, thresholds_fold: combo.thr_fold, floor_binds: combo.floor_binds,
                n_objects: objs.size(),
                n_by_class: new LinkedHashMap(ENDOCRINE_CLASSES.collectEntries { c -> [c, objs.count { cls(it) == c }] }),
                n_by_route: new TreeMap([1: objs.count { m(it, 'Endocrine route') == 1 }, 2: objs.count { m(it, 'Endocrine route') == 2 },
                                         3: objs.count { m(it, 'Endocrine route') == 3 }]),
                n_by_class_route: byClassRoute,
                n_cd99_only_objects: objs.count { m(it, 'Endocrine route') == 1 },
                n_single_only_by_channel: new LinkedHashMap(CH.collectEntries { ch -> [ch, singleOnly(ch)] }),
                ins_only_singles_per_mm2: singleOnly('INS') / mm2,
                n_multi_core: objs.count { m(it, 'Multi-core flag') == 1 },
                n_multi_core_islet: objs.count { m(it, 'Multi-core flag') == 1 && cls(it) == 'Islet' },
                n_neck_split_objects_in: nSplitObjects, n_necks_cut: nNecksCut,
                n_neck_split_objects_out: objs.count { m(it, 'Neck split') == 1 },
                n_clamped_only: objs.count { m(it, 'Clamped-only flag') == 1 },
                n_single_copos_removed: dropped.size(),
                n_bridge_flag: objs.count { m(it, 'Bridge flag') == 1 },
                n_bridge_flag_islet: objs.count { m(it, 'Bridge flag') == 1 && cls(it) == 'Islet' },
                n_tissue_border_flag: objs.count { m(it, 'Tissue border flag') == 1 },
                n_pieces: base.n_pieces, n_hormone_parts: base.hPolys.size(), n_hormone_parts_after_opening: base.n_hormone_parts_after_opening,
                n_hormone_cores: base.coreList.size(), n_cd99_only_parts: base.cOnlyList.size(),
                n_hole_swallowed: base.n_hole_swallowed,
                endocrine_area_um2: endoArea, tissue_area_um2: tissueAreaUm2,
                endocrine_area_fraction: endoArea / tissueAreaUm2,
                endocrine_fraction_flag: (endoArea / tissueAreaUm2) > (P.MAX_ENDOCRINE_FRACTION as double),
                islet_area_um2: objs.findAll { cls(it) == CLASS_ISLET }.sum(0.0d) { m(it, "Area ${UM}^2") },
                largest_object_um2: objs.isEmpty() ? 0 : objs.collect { m(it, "Area ${UM}^2") }.max(),
                geometry_seconds: base.base_seconds + (System.currentTimeMillis() - t0) / 1000.0])
        [objs: objs, summary: summary]
    }

    // JSON-safe copy (NaN/Inf -> null; GStrings -> String)
    def plain
    plain = { v ->
        if (v == null) return null
        if (v instanceof CharSequence) return v.toString()
        if (v instanceof Double || v instanceof Float) return (Double.isNaN(v as double) || Double.isInfinite(v as double)) ? null : v
        if (v instanceof Map) { def mm = new LinkedHashMap(); v.each { k, x -> mm[String.valueOf(k)] = plain(x) }; return mm }
        if (v instanceof Collection) return v.collect { plain(it) }
        return v
    }
    def writeJson = { String path, Map content ->
        def f = new File(path); f.getParentFile()?.mkdirs()
        f.text = GsonTools.getInstance(true).toJson(plain(content))
    }
    def baseOut = (P.WRITE_GEOJSON as String) ? (P.WRITE_GEOJSON as String).replaceAll(/(?i)\.geojson$/, '') : null
    def common = new LinkedHashMap([image: imageName, tag: tag, method: 'drotar_cd99', anchors: anchors,
                                    clamp_mu: clampMu, clamp_sigma: clampSd, clamp_fraction: clampFrac,
                                    anchors_csv: P.ANCHORS_CSV, anchors_method_version: anchorVersion,
                                    annotation_hash_csv: csvHash, annotation_hash_now: nowHash, anchor_check_issues: checkIssues,
                                    params: P, trace_seconds: traceSec, n_traced_masks: layout.size()])

    if (sweep) {
        // --------------------------- sweep: write, never touch the hierarchy ---------------------------
        def results = []
        combos.each { Map combo -> openHs.each { double rH ->
            def base = buildBase(combo, rH)
            necks.each { boolean neck -> coposs.each { boolean copos ->
                def res1 = finishObjects(combo, base, rH, neck, copos)
                String f = "${baseOut}__tC${fmt(combo.tC)}_tH${fmt(combo.tH)}_k${fmt(combo.k)}" + (rH > 0 ? "_rH${fmt(rH)}" : '') +
                        (neck ? '_neck' : '') + (copos ? '_copos' : '') + '.geojson'
                if (!res1.objs.isEmpty()) exportObjectsToGeoJson(res1.objs, f, 'FEATURE_COLLECTION')
                res1.summary.geojson = res1.objs.isEmpty() ? null : f
                results << res1.summary
                println "[01v3] sweep tC=${fmt(combo.tC)} tH=${fmt(combo.tH)} k=${fmt(combo.k)} rH=${fmt(rH)} neck=${neck} copos=${copos}: " +
                        "${res1.summary.n_by_class} routes ${res1.summary.n_by_route} multi-core ${res1.summary.n_multi_core} " +
                        "INS-only singles ${String.format('%.2f', res1.summary.ins_only_singles_per_mm2)}/mm^2 " +
                        "endocrine ${String.format('%.2f%%', 100 * res1.summary.endocrine_area_fraction)} (${String.format('%.1f', res1.summary.geometry_seconds)} s)"
            } }
        } }
        common.sweep = results
        common.total_seconds = (System.currentTimeMillis() - tStart) / 1000.0
        writeJson("${baseOut}.summary.json", common)
        // tell a run_v3 footer there is nothing to save
        binding.setVariable('V3_RESULT', [ok: true, changed: false, message: "sweep of ${results.size()} combinations; hierarchy untouched"])
        println "[01v3] sweep of ${results.size()} combination(s) written to ${baseOut}.summary.json; hierarchy untouched"
    } else {
        // ---------------------------------- production ----------------------------------
        def combo = combos[0]
        double rH = openHs[0]
        def built = finishObjects(combo, buildBase(combo, rH), rH, necks[0], coposs[0])
        def objs = built.objs
        def summary = built.summary
        // swap only now, so a failure above leaves the hierarchy as it was
        def old = getAnnotationObjects().findAll { it.getPathClass()?.toString() in ENDOCRINE_CLASSES }
        def oldBands = getAnnotationObjects().findAll { it.getPathClass()?.toString() == CLASS_ISLET_EXPANDED }
        common.n_removed_previous = new LinkedHashMap(old.countBy { it.getPathClass().toString() })
        common.n_removed_previous[CLASS_ISLET_EXPANDED] = oldBands.size()
        if (!old.isEmpty() || !oldBands.isEmpty()) {
            println "[01v3] Removing previous objects ${common.n_removed_previous} (child detections move up)"
            removeObjects(old + oldBands, true)
        }
        addObjects(objs)
        // Detections of a removed object move to ITS parent; bands (and objects added by addObjects) live at the root,
        // so their cells would sit directly under the root until script 8. Give them to the region that contains
        // their nucleus centroid (Tissue first, then pLN; the test point script 8 uses).
        def root = hierarchy.getRootObject()
        def orphans = root.getChildObjects().findAll { it.isDetection() }
        def reparented = new LinkedHashMap([(CLASS_TISSUE): 0, pLN: 0, left_at_root: 0])
        if (!orphans.isEmpty()) {
            def claimed = Collections.newSetFromMap(new IdentityHashMap())
            [CLASS_TISSUE, 'pLN'].each { rc ->
                getAnnotationObjects().findAll { it.getPathClass()?.toString() == rc && it.getROI()?.isArea() }.each { region ->
                    def prep = PreparedGeometryFactory.prepare(region.getROI().getGeometry())
                    def kids = orphans.findAll { d -> if (claimed.contains(d)) return false
                        def r = PathObjectTools.getROI(d, true)
                        prep.contains(gf.createPoint(new Coordinate(r.getCentroidX(), r.getCentroidY()))) }
                    if (!kids.isEmpty()) { region.addChildObjects(kids); claimed.addAll(kids); reparented[rc] = (reparented[rc] as int) + kids.size() }
                }
            }
            reparented.left_at_root = orphans.size() - claimed.size()
            println "[01v3] Re-parented detections that were left under the root: ${reparented}"
        }
        common.reparented_root_detections = reparented
        fireHierarchyUpdate()
        // per-image classifiers, for visual review in the GUI (same thresholds as the traced server)
        if (P.SAVE_CLASSIFIERS as boolean) {
            def dir = Paths.get(P.OUT_CLASSIFIER_DIR as String, tag)
            dir.toFile().mkdirs()
            def ign = PathClass.fromString('Ignore*')
            def thr = { double t -> [ImageOps.Threshold.threshold(t), ImageOps.Core.ensureType(PixelType.UINT8)] }
            def save = { String name, op, Map labels ->
                PixelClassifiers.writeClassifier(PixelClassifiers.createClassifier(op, calRes, labels), dir.resolve(name + '.json'))
            }
            def tf = combo.thr_fold
            save('drotar_cd99_C', foldOp(['CD99']).appendOps(*thr(tf.CD99 as double)), [0: ign, 1: PathClass.fromString('mask_CD99')])
            save('drotar_cd99_H', foldOp(CH).appendOps(hRatioOp(tf.INS as double, tf.GCG as double, tf.SST as double), *thr(1.0d)),
                    [0: ign, 1: PathClass.fromString('mask_Hormone')])
            HORMONES.each { ch -> save("drotar_cd99_${ch}", foldOp([ch]).appendOps(*thr(tf[ch] as double)), [0: ign, 1: PathClass.fromString("mask_${ch}")]) }
            // one overlay for review: 1 CD99 only, 2 hormone only, 3 both
            save('drotar_cd99_routes', foldOp(CH).appendOps(
                    ImageOps.Core.splitMerge(ImageOps.Core.sequential(ImageOps.Channels.extract(0), ImageOps.Core.multiply(1.0d / (tf.CD99 as double))),
                            hRatioOp(tf.INS as double, tf.GCG as double, tf.SST as double)),
                    ImageOps.Threshold.threshold(1.0d, 1.0d), ImageOps.Core.multiply(1.0d, 2.0d), ImageOps.Channels.sum(),
                    ImageOps.Core.ensureType(PixelType.UINT8)),
                    [0: ign, 1: PathClass.fromString('route_CD99only'), 2: PathClass.fromString('route_HormoneOnly'),
                     3: PathClass.fromString('route_Both')])
            common.classifier_dir = dir.toString()
            println "[01v3] Saved per-image classifiers to ${dir}"
        }
        common.result = summary
        common.total_seconds = (System.currentTimeMillis() - tStart) / 1000.0
        if (baseOut) {
            if (!objs.isEmpty()) exportObjectsToGeoJson(objs, "${baseOut}.geojson", 'FEATURE_COLLECTION')
            writeJson("${baseOut}.summary.json", common)
            println "[01v3] Wrote ${baseOut}.geojson and ${baseOut}.summary.json"
        }
        if (summary.endocrine_fraction_flag)
            println String.format("[01v3] WARNING: endocrine area is %.1f%% of tissue (> %.0f%%): check thresholds/anchors",
                    100 * summary.endocrine_area_fraction, 100 * (P.MAX_ENDOCRINE_FRACTION as double))
        // run_v3 provenance picks this up; 'changed' stays true so the footer may save
        binding.setVariable('V3_RESULT', [ok: true, changed: true, message: "drotar_cd99 ${summary.n_by_class}",
                                          details: plain(common)])
        println "[01v3] ${imageName}: ${summary.n_by_class} routes ${summary.n_by_route} multi-core=${summary.n_multi_core} " +
                String.format("endocrine=%.2f%% of %.1f mm^2, %.1f s", 100 * summary.endocrine_area_fraction,
                        tissueAreaUm2 / 1e6, common.total_seconds)
    }
}
