import qupath.lib.analysis.heatmaps.DensityMaps
import qupath.lib.objects.PathObjectFilter
import qupath.lib.objects.PathObjectPredicates
import qupath.lib.images.servers.PixelCalibration

// --- Parameters ---
String  targetClassName    = "Other: Gal3Norm"
double  densityRadiusUm    = 50.0    // count target cells within this radius (merge/smoothing scale)
double  densityPixelSizeUm = 8.0    // density-map resolution; pinned to avoid null-pixel-size full-res contouring
double  countThreshold     = 2.0     // min target cells within radius for a pixel to be "dense"
boolean deleteExistingTargetAnnotations = true   // clear target-class annotations from earlier runs

def imageData = getCurrentImageData()
def hierarchy = imageData.getHierarchy()
def cal = imageData.getServer().getPixelCalibration()
if (!cal.hasPixelSizeMicrons()) { println "ERROR: image has no pixel size"; return }

def targetClass = getPathClass(targetClassName)

int nTarget = getDetectionObjects().count { it.getPathClass() == targetClass }
println "Target detections (${targetClassName}): ${nTarget}"
if (nTarget == 0) { println "No detections of target class; check the exact name."; return }

// Clear target-class annotations left over from prior runs, so re-runs are clean
if (deleteExistingTargetAnnotations) {
    def old = getAnnotationObjects().findAll { it.getPathClass() == targetClass }
    if (!old.isEmpty()) {
        old.each { it.setLocked(false) }
        removeObjects(old, false)
        println "Removed ${old.size()} existing ${targetClassName} annotations"
    }
}

// Predicate: detections classified exactly as the target class
def predicate = PathObjectPredicates.filter(PathObjectFilter.DETECTIONS_ALL)
        .and(PathObjectPredicates.exactClassification(targetClass))

// Pin the density-map pixel size explicitly
def densityCal = new PixelCalibration.Builder()
        .pixelSizeMicrons(densityPixelSizeUm, densityPixelSizeUm)
        .build()

// SUM = raw count of target cells within radius, per density-map pixel
def builder = DensityMaps.builder(predicate)
        .type(DensityMaps.DensityMapType.SUM)
        .radius(densityRadiusUm)
        .pixelSize(densityCal)

def densityServer = builder.buildServer(imageData)
println "Density map: ${densityServer.getWidth()} x ${densityServer.getHeight()} px @ ${densityPixelSizeUm} µm/px, ${densityServer.nChannels()} channel(s)"
println "  channel 0: ${densityServer.getMetadata().getChannels()[0].getName()}"

// Threshold channel 0 -> annotations. No SPLIT: one multipolygon covering all dense regions.
deselectAll()
def before = getAnnotationObjects() as Set
DensityMaps.threshold(hierarchy, densityServer, 0, countThreshold)
def created = (getAnnotationObjects() as Set) - before
println "Created ${created.size()} annotation(s) at count >= ${countThreshold}"

// Assign the target class and lock; report area
double totalAreaUm2 = 0
double pxAreaUm2 = cal.getPixelWidthMicrons() * cal.getPixelHeightMicrons()
created.each {
    it.setPathClass(targetClass)
    it.setLocked(true)
    totalAreaUm2 += it.getROI().getArea() * pxAreaUm2
}
fireHierarchyUpdate()
println String.format("Total annotated area: %.0f µm² (%.2f mm²)", totalAreaUm2, totalAreaUm2 / 1e6)