/*
 * fixture_crop.groovy -- write a small full-resolution, all-channel OME-TIFF crop of one project
 * image, for testing run_v3.sh without touching real data.
 *
 * Read-only on the real project: it loads a PRIVATE COPY of project.qpproj (never the master),
 * builds the entry's image server from that copy (so channel names and calibration come from the
 * project entry, exactly as the pipeline sees them) and reads pixels. It never opens the entry's
 * .qpdata and saves nothing.
 *
 * Headless, no -p/-i (the script opens the project copy itself):
 *   QuPath script -a <project copy .qpproj> -a "<image name>" -a x -a y -a w -a h \
 *                 -a <out.ome.tif> -a <pixel size um or 'keep'> fixture_crop.groovy
 *
 * The pixel size argument exists because 6450panc's entry says 0.5077663810243286 um while every
 * other image (and scripts 0/1) use 0.5078; the fixture should mimic the common case.
 */

import qupath.lib.images.servers.ImageServerMetadata
import qupath.lib.images.writers.ome.OMEPyramidWriter
import qupath.lib.projects.ProjectIO
import java.awt.image.BufferedImage

if (args == null || args.length < 8)
    throw new IllegalArgumentException('usage: <project copy> <image name> <x> <y> <w> <h> <out.ome.tif> <pixel um|keep>')
def (projPath, imageName) = [args[0], args[1]]
int x = args[2] as int, y = args[3] as int, w = args[4] as int, h = args[5] as int
def outPath = args[6]
def pixelArg = args[7]

if (new File(projPath).getCanonicalPath() == new File('/home/smith6jt/IO60panc2nd/project.qpproj').getCanonicalPath())
    throw new IllegalArgumentException('refusing to open the master project; pass a private copy')

def project = ProjectIO.loadProject(new File(projPath), BufferedImage)
def entry = project.getImageList().find { it.getImageName() == imageName }
if (entry == null) throw new IllegalArgumentException("no entry named '${imageName}'")

def server = entry.getServerBuilder().build()
def md = server.getMetadata()
def cal = md.getPixelCalibration()
println "[fixture_crop] ${imageName}: ${server.getWidth()}x${server.getHeight()}, ${md.getChannels().size()} channels, " +
        "pixel ${cal.getPixelWidth()} ${cal.getPixelWidthUnit()}, type ${md.getPixelType()}"
if (x < 0 || y < 0 || x + w > server.getWidth() || y + h > server.getHeight())
    throw new IllegalArgumentException("crop ${x},${y},${w}x${h} is outside the image")

if (pixelArg != 'keep') {
    double px = pixelArg as double
    server.setMetadata(new ImageServerMetadata.Builder(md).pixelSizeMicrons(px, px).build())
    println "[fixture_crop] pixel size set to ${px} um for the fixture (entry value ${cal.getPixelWidth()})"
}

long t0 = System.currentTimeMillis()
new File(outPath).getParentFile().mkdirs()
new OMEPyramidWriter.Builder(server)
        .region(x, y, w, h)
        .tileSize(512)
        .channelsPlanar()
        .losslessCompression()
        .dyadicDownsampling()
        .parallelize()
        .build()
        .writePyramid(outPath)
server.close()
println String.format("[fixture_crop] wrote %s (%.1f MB) in %.1f s", outPath,
        new File(outPath).length() / 1e6, (System.currentTimeMillis() - t0) / 1000.0)
