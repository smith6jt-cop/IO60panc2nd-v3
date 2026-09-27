/*
 * fixture_project.groovy -- create a one-image QuPath project for testing run_v3.sh.
 *
 * Creates <project dir>/project.qpproj with the fixture OME-TIFF as its only entry (image type
 * FLUORESCENCE), imports the fixture annotations (Tissue, pLN, a stale 'Other: Gal3Norm', and the
 * v2 islets that fall inside the crop, all in crop pixel coordinates) and saves. Refuses to write
 * anywhere but a fresh directory under /scratch/io60/fixture.
 *
 * Headless, no -p/-i:
 *   QuPath script -a <project dir> -a <fixture.ome.tif> -a <annotations.geojson> fixture_project.groovy
 */

import qupath.lib.images.ImageData
import qupath.lib.images.servers.ImageServers
import qupath.lib.io.PathIO
import qupath.lib.projects.Projects
import java.awt.image.BufferedImage

if (args == null || args.length < 3)
    throw new IllegalArgumentException('usage: <project dir> <fixture.ome.tif> <annotations.geojson>')
def projDir = new File(args[0]).getCanonicalFile()
def imagePath = new File(args[1]).getCanonicalFile()
def annPath = new File(args[2]).getCanonicalFile()

if (!projDir.getPath().startsWith('/scratch/io60/fixture/'))
    throw new IllegalArgumentException("fixture projects live under /scratch/io60/fixture/, not ${projDir}")
if (projDir.exists() && projDir.list().length > 0)
    throw new IllegalArgumentException("${projDir} is not empty; delete it first")
projDir.mkdirs()

def project = Projects.createProject(projDir, BufferedImage)
def server = ImageServers.buildServer(imagePath.toURI())
def entry = project.addImage(server.getBuilder())
entry.setImageName(imagePath.getName())
server.close()

def imageData = entry.readImageData()
imageData.setImageType(ImageData.ImageType.FLUORESCENCE)
def md = imageData.getServer().getMetadata()
println "[fixture_project] entry '${entry.getImageName()}' (id ${entry.getID()}): ${md.getWidth()}x${md.getHeight()}, " +
        "${md.getChannels().size()} channels, pixel ${md.getPixelCalibration().getPixelWidth()} " +
        "${md.getPixelCalibration().getPixelWidthUnit()}"

def objects = PathIO.readObjects(annPath)
imageData.getHierarchy().insertPathObjects(objects)   // placed by ROI containment: islets and pLN nest under Tissue, as in the real project
def byClass = imageData.getHierarchy().getAnnotationObjects().countBy { it.getPathClass()?.toString() ?: 'Unclassified' }
def tissue = imageData.getHierarchy().getAnnotationObjects().find { it.getPathClass()?.toString() == 'Tissue' }
println "[fixture_project] annotations ${byClass}; children of Tissue: ${tissue?.getChildObjects()?.size()}"

entry.saveImageData(imageData)
project.syncChanges()
println "[fixture_project] saved ${new File(projDir, 'project.qpproj')}"
