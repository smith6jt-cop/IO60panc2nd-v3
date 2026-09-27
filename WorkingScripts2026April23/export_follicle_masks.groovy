/**
 * Export binary masks of Follicle annotations for each image.
 * Run via: QuPath script -p /path/to/project scripts/export_follicle_masks.groovy
 */
import javax.imageio.ImageIO
import java.awt.Color
import java.awt.image.BufferedImage
import java.awt.RenderingHints

def imageData = getCurrentImageData()
if (imageData == null) {
    print("No image data available")
    return
}

def server = imageData.getServer()
def hierarchy = imageData.getHierarchy()

// Downsample for manageable file size
double downsample = 16.0
int w = (int) Math.ceil(server.getWidth() / downsample)
int h = (int) Math.ceil(server.getHeight() / downsample)

// Get Follicle annotations
def follicles = hierarchy.getAnnotationObjects().findAll {
    it.getPathClass() != null && it.getPathClass().getName() == 'Follicle'
}

if (follicles.isEmpty()) {
    print("No Follicle annotations found in " + server.getMetadata().getName())
    return
}

// Create binary mask
def mask = new BufferedImage(w, h, BufferedImage.TYPE_BYTE_GRAY)
def g2d = mask.createGraphics()
g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_OFF)
g2d.setColor(Color.BLACK)
g2d.fillRect(0, 0, w, h)
g2d.setColor(Color.WHITE)

// Scale coordinates to match downsampled image
def transform = java.awt.geom.AffineTransform.getScaleInstance(1.0 / downsample, 1.0 / downsample)
g2d.setTransform(transform)

// Fill each Follicle ROI
for (def ann : follicles) {
    def roi = ann.getROI()
    def shape = roi.getShape()
    g2d.fill(shape)
}
g2d.dispose()

// Save
def imageName = server.getMetadata().getName().replaceAll(/\.[^.]+$/, '').replaceAll(/ .*$/, '')
def outputDir = new File(buildFilePath(PROJECT_BASE_DIR, 'analysis', 'masks'))
outputDir.mkdirs()
def outputFile = new File(outputDir, imageName + '_follicle_mask.png')
ImageIO.write(mask, 'PNG', outputFile)
print("Saved: " + outputFile.getAbsolutePath() + " (" + follicles.size() + " follicles, " + w + "x" + h + ")")
