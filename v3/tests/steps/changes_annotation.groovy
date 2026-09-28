// Test step (stands in for 'preflight'): a real in-memory change that the footer wants to save. Run with a
// huge --save-margin-min, the footer must skip the save (too close to the step's deadline), record why in a
// failing sentinel, and leave data.qpdata untouched.
def roi = qupath.lib.roi.ROIs.createRectangleROI(100, 100, 50, 50, qupath.lib.regions.ImagePlane.getDefaultPlane())
addObject(qupath.lib.objects.PathObjects.createAnnotationObject(roi, getPathClass('Other')))
V3_RESULT = [ok: true, changed: true, message: 'added one annotation', details: [:]]
