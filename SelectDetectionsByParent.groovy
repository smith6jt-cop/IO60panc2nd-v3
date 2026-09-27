def rois = getAnnotationObjects()
        .findAll { it.getPathClass()?.toString() == "pLN" && it.getROI() != null }
        .collect { it.getROI() }

def inside = getDetectionObjects().findAll { d ->
    def r = d.getROI()
    r != null && rois.any { it.contains(r.getCentroidX(), r.getCentroidY()) }
}
resetSelection()
selectObjects(inside)