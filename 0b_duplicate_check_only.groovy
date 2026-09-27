/*
 * 0b_duplicate_check_only.groovy
 *
 * Re-runs the post-detection duplicate check from 0_cell_detection.groovy, against whatever
 * detections are already present on the currently open image. It removes redundant detections,
 * retaining one object from each connected group of matching near-duplicates. Use this after a
 * run of 0_cell_detection.groovy that got through InstanSeg but crashed (or was interrupted),
 * so InstanSeg does not need to be repeated.
 *
 * Saves the image data at the end so the detections are persisted this time.
 */

import java.util.IdentityHashMap

double DUP_RADIUS_UM = 1.0
double DUP_AREA_TOL  = 0.15

def server = getCurrentServer()
if (server == null) { println "ERROR: no image open."; return }

def cells = getDetectionObjects()
if (cells.isEmpty()) {
    println "ERROR: no detections found on this image."
    return
}
println String.format("[dup-check] %,d existing detections", cells.size())

def cal = server.getPixelCalibration()
double pxSize = cal.hasPixelSizeMicrons() ? cal.getAveragedPixelSizeMicrons() : 1.0
double rPx = DUP_RADIUS_UM / pxSize

// Plain map (not withDefault): withDefault's get() inserts missing keys as a side effect, which
// causes a ConcurrentModificationException if a neighbor lookup mutates the map mid-iteration.
def buckets = [:]
for (c in cells) {
    def roi = c.getROI()
    if (roi == null) continue
    def key = [(int) (roi.getCentroidX() / rPx), (int) (roi.getCentroidY() / rPx)]
    buckets.computeIfAbsent(key) { [] } << c
}

def parent = new IdentityHashMap()
for (c in cells) parent.put(c, c)

def findRoot
findRoot = { object ->
    def root = parent.get(object)
    if (!root.is(object)) {
        root = findRoot(root)
        parent.put(object, root)
    }
    root
}
def joinDuplicates = { first, second ->
    def firstRoot = findRoot(first)
    def secondRoot = findRoot(second)
    if (!firstRoot.is(secondRoot)) parent.put(secondRoot, firstRoot)
}

int dupPairs = 0
for (entry in new ArrayList(buckets.entrySet())) {
    def (bx, by) = entry.key
    def near = []
    for (dx in -1..1) for (dy in -1..1) near.addAll(buckets.getOrDefault([bx + dx, by + dy], []))
    for (a in entry.value) {
        def ra = a.getROI()
        for (b in near) {
            if (System.identityHashCode(b) <= System.identityHashCode(a)) continue
            def rb = b.getROI()
            double d = Math.hypot(ra.getCentroidX() - rb.getCentroidX(),
                                  ra.getCentroidY() - rb.getCentroidY())
            if (d > rPx) continue
            double aa = ra.getArea(), ab = rb.getArea()
            if (Math.abs(aa - ab) / Math.max(Math.max(aa, ab), 1e-9) < DUP_AREA_TOL) {
                dupPairs++
                joinDuplicates(a, b)
            }
        }
    }
}

if (dupPairs > 0) {
    def duplicates = cells.findAll { !findRoot(it).is(it) }
    println "[dup-check] ${dupPairs} duplicate detection pair(s) found."
    println "[dup-check] Removing ${duplicates.size()} redundant detection(s); retaining ${cells.size() - duplicates.size()}."
    removeObjects(duplicates, true)
    fireHierarchyUpdate()

    def remaining = getDetectionObjects()
    println "[dup-check] ${remaining.size()} detections remain after deduplication."
} else {
    println "[dup-check] PASS (0 pairs within ${DUP_RADIUS_UM} um)"
}

getProjectEntry().saveImageData(getCurrentImageData())
println "[dup-check] Saved image data."

