/**
 * Rebuild project thumbnails for multichannel (non-RGB) images. QuPath v0.7.0.
 *
 * Why the stock thumbnails are white:
 *   ProjectCommands.getThumbnailRGB(server) delegates to
 *   ProjectImportImagesCommand.getThumbnailRGB(server, null) with a null ImageDisplay.
 *   That builds a fresh ImageDisplay, which selects EVERY DirectServerChannelInfo
 *   (ImageDisplay.updateChannelOptions), auto-stretches each one independently to
 *   PathPrefs.autoBrightnessContrastSaturationPercentProperty() (default 0.1%), then
 *   blends them with updateRGBAdditive(). Summing N stretched channels clips to 255,255,255.
 *
 * Run with the plain Run button (this script loops over the project itself).
 * Do NOT use "Run for project".
 */

import qupath.lib.display.ChannelDisplayMode
import qupath.lib.display.DirectServerChannelInfo
import qupath.lib.display.ImageDisplay
import qupath.lib.images.ImageData
import java.awt.RenderingHints
import java.awt.image.BufferedImage

// ------------------------------- config -------------------------------
// Channel names as they appear in the QPTIFF metadata, WITHOUT the " (Cn)" postfix
// that ChannelDisplayInfo.getName() appends. Every extra channel is summed into the
// RGB output, so keep this to 1-3. Run once with a wrong name to print the list.
def CHANNEL_NAMES      = ['DAPI', 'Vimentin']
def SATURATION_PERCENT = 0.5                       // % saturated pixels per channel (import default: 0.1)
def MODE               = ChannelDisplayMode.COLOR  // COLOR | INVERTED_COLOR | GRAYSCALE | INVERTED_GRAYSCALE
def MAX_W              = 1000                      // QuPath's own thumbnail bounds
def MAX_H              = 600
// ----------------------------------------------------------------------

def project = getProject()
int nDone = 0

for (entry in project.getImageList()) {
    def server = entry.getServerBuilder().build()
    try {
        if (server.isRGB()) {
            println "Skipped (RGB): ${entry.getImageName()}"
            continue
        }

        // Histograms are built from a pyramid level of <= 2^20 px, not the thumbnail
        def imageData = new ImageData<BufferedImage>(server, ImageData.ImageType.FLUORESCENCE)
        def display = ImageDisplay.create(imageData)

        def direct = display.availableChannels().findAll { it instanceof DirectServerChannelInfo }
        def channels = direct.findAll { it.getOriginalChannelName() in CHANNEL_NAMES }
        if (channels.isEmpty()) {
            println "No match in ${entry.getImageName()}; available: ${direct*.getOriginalChannelName()}"
            continue
        }
        channels.each { display.autoSetDisplayRange(it, SATURATION_PERCENT / 100.0) }
        // Fixed range instead of percentile stretch:
        // channels.each { display.setMinMaxDisplay(it, 0f, 3000f) }

        def img = server.getDefaultThumbnail((int)(server.nZSlices() / 2), 0)
        def rgb = ImageDisplay.applyTransforms(img, null, channels, MODE)

        double scale = Math.min((double)MAX_W / rgb.getWidth(), (double)MAX_H / rgb.getHeight())
        if (scale < 1) {
            def out = new BufferedImage((int)(rgb.getWidth() * scale), (int)(rgb.getHeight() * scale), BufferedImage.TYPE_INT_RGB)
            def g2d = out.createGraphics()
            g2d.setRenderingHint(RenderingHints.KEY_INTERPOLATION, RenderingHints.VALUE_INTERPOLATION_BILINEAR)
            g2d.drawImage(rgb, 0, 0, out.getWidth(), out.getHeight(), null)
            g2d.dispose()
            rgb = out
        }

        entry.setThumbnail(rgb)   // writes <project>/data/<entryID>/thumbnail.jpg immediately
        nDone++
        println "Rebuilt: ${entry.getImageName()}  ${channels*.getOriginalChannelName()}"
    } finally {
        server.close()
    }
}

getQuPath().refreshProject()
println "Done: ${nDone} thumbnail(s)"