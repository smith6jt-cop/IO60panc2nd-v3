/*
 * islet_benchmark_tiles.groovy -- coded, blinded 1 mm x 1 mm tiles for calibrating the v3 islet detector.
 *
 * For the current image, three kinds of tile (full resolution, TILE_UM x TILE_UM):
 *   random            N_RANDOM tiles at random positions inside Tissue (>= RANDOM_MIN_TISSUE covered): the
 *                     unbiased sample that estimates prevalence and misses;
 *   object:<class>:r<route>
 *                     N_OBJECT tiles centred on detector objects, round-robin over class x route strata
 *                     (Islet / EndocrineCluster / EndocrineSingle x route 3, 2, 1), so rare strata (route 1 =
 *                     CD99 only, singles) are represented;
 *   exo_cd99          N_EXO exocrine tiles with the most CD99-bright pixels and no hormone signal: no Islet
 *                     object and hormone-positive area <= EXO_MAX_HORMONE_FRAC (scattered single endocrine cells
 *                     and tiny clusters are everywhere, ~20-40 per mm^2, so "none at all" selects almost nothing).
 *                     The bridge failure case: CD99-bright exocrine cells that could join islets into one object.
 *                     Scored at SCORE_UM (2 um keeps single CD99-bright cells above TAU_CD99; 4 um dilutes them).
 * Tiles of one image overlap by at most MAX_OVERLAP_FRAC of a tile.
 *
 * Output (OUT_DIR, default /scratch/io60/benchmark):
 *   tiles/<code>.ome.tif    DAPI, CD99, INS, GCG, SST crops, UINT16 counts as scanned (no anchors applied),
 *                           OME image name = <code>; the file is checked for the source image name/tag, so
 *                           annotators see neither donor nor disease status
 *   pred/<code>.geojson     detector objects clipped to the tile, in TILE pixel coordinates, with class,
 *                           route and flags (plus the Tissue/Artifact outline clipped to the tile), for scoring
 *   key.tsv                 code, image, x, y, w, h, tile_kind   -- the ONLY link back to the source image
 *   runs/<tag>.json         parameters and per-tile selection details (e.g. exo_cd99 scores)
 * Codes are random (SecureRandom), assigned after shuffling, and unique across key.tsv.
 *
 * Detector objects come from the hierarchy (after 1_panc_objects_v3.groovy), or from OBJECTS_GEOJSON (a
 * WRITE_GEOJSON output of script 1, so tiles can be cut without saving any image data). Read-only: the image
 * data is never modified or saved. An image already in key.tsv is skipped unless REPLACE=true (which removes
 * its old rows and files first).
 *
 * Headless:  QuPath script -p <project> -i "<image>" [-a KEY=VALUE ...] islet_benchmark_tiles.groovy
 */

import qupath.lib.images.ImageData
import qupath.lib.images.servers.ColorTransforms
import qupath.lib.images.servers.PixelCalibration
import qupath.lib.images.servers.PixelType
import qupath.lib.images.writers.ome.OMEPyramidWriter
import qupath.lib.io.GsonTools
import qupath.lib.io.PathIO
import qupath.lib.objects.PathObjects
import qupath.lib.objects.hierarchy.PathObjectHierarchy
import qupath.lib.regions.RegionRequest
import qupath.lib.roi.GeometryTools
import qupath.lib.roi.RoiTools
import qupath.opencv.ops.ImageOps
import org.locationtech.jts.geom.Envelope
import org.locationtech.jts.geom.util.AffineTransformation
import org.locationtech.jts.geom.util.PolygonExtracter
import org.locationtech.jts.geom.prep.PreparedGeometryFactory
import org.locationtech.jts.index.strtree.STRtree
import org.locationtech.jts.operation.union.UnaryUnionOp
import java.awt.geom.AffineTransform
import java.awt.image.BufferedImage
import java.nio.channels.FileChannel
import java.nio.file.Files
import java.nio.file.Paths
import java.nio.file.StandardOpenOption
import java.security.SecureRandom

// ==================== PARAMETERS ====================
String OUT_DIR              = '/scratch/io60/benchmark'
int N                       = 8        // tiles per kind, unless overridden below
int N_RANDOM                = -1       // -1 = N
int N_OBJECT                = -1
int N_EXO                   = -1
double TILE_UM              = 1000.0
double RANDOM_MIN_TISSUE    = 0.5      // fraction of the tile inside Tissue
double EXO_MIN_TISSUE       = 0.9
double MAX_OVERLAP_FRAC     = 0.25     // between tiles of this image
double EXO_MAX_HORMONE_FRAC = 0.005    // "no hormone signal": hormone-positive area fraction at SCORE_UM
double SCORE_UM             = 2.0      // scale for exo_cd99 scoring (~pyramid level 2)
double TAU_CD99             = 3.0      // folds for exo_cd99 scoring; match the detector's provisional values
double TAU_H                = 5.0
String ANCHORS_CSV          = '/home/smith6jt/IO60panc2nd/results/v3/islet_anchors.csv'
String OBJECTS_GEOJSON      = ''       // '' = detector objects from the hierarchy
long SEED                   = -1L      // -1 = derived from the image name (reproducible positions)
boolean REPLACE             = false
String CHANNELS             = 'DAPI,CD99,INS,GCG,SST'   // never IAPP (failed QC) or HLA-A
// ====================================================

def P = new LinkedHashMap([OUT_DIR: OUT_DIR, N: N, N_RANDOM: N_RANDOM, N_OBJECT: N_OBJECT, N_EXO: N_EXO, TILE_UM: TILE_UM,
        RANDOM_MIN_TISSUE: RANDOM_MIN_TISSUE, EXO_MIN_TISSUE: EXO_MIN_TISSUE, MAX_OVERLAP_FRAC: MAX_OVERLAP_FRAC,
        EXO_MAX_HORMONE_FRAC: EXO_MAX_HORMONE_FRAC, SCORE_UM: SCORE_UM, TAU_CD99: TAU_CD99, TAU_H: TAU_H,
        ANCHORS_CSV: ANCHORS_CSV, OBJECTS_GEOJSON: OBJECTS_GEOJSON, SEED: SEED, REPLACE: REPLACE, CHANNELS: CHANNELS])
def scriptArgs = (binding.hasVariable('args') && args != null) ? (args as List).collect { it.toString() } : []
scriptArgs.each { String a ->
    int eq = a.indexOf('=')
    if (eq <= 0) throw new IllegalArgumentException("script argument '${a}' is not KEY=VALUE")
    String k = a.substring(0, eq).trim(), v = a.substring(eq + 1).trim()
    if (!P.containsKey(k)) throw new IllegalArgumentException("unknown parameter '${k}'; known: ${P.keySet()}")
    def old = P[k]
    P[k] = (old instanceof Double) ? Double.parseDouble(v) : (old instanceof Integer) ? Integer.parseInt(v) :
            (old instanceof Long) ? Long.parseLong(v) : (old instanceof Boolean) ? Boolean.parseBoolean(v) : v
}
int nRandom = (P.N_RANDOM as int) >= 0 ? P.N_RANDOM as int : P.N as int
int nObject = (P.N_OBJECT as int) >= 0 ? P.N_OBJECT as int : P.N as int
int nExo    = (P.N_EXO as int) >= 0 ? P.N_EXO as int : P.N as int

def imageData = getCurrentImageData()
def server = imageData.getServer()
def entry = getProjectEntry()
String imageName = entry?.getImageName() ?: server.getMetadata().getName()
String tag = imageName.replaceAll(/[^A-Za-z0-9._-]+/, '_')
String stem = imageName.replaceAll(/(?i)\.(er\.)?qptiff.*$/, '').replaceAll(/(?i)\.ome\.tif+.*$/, '')
def cal = server.getPixelCalibration()
double pxUm = cal.getAveragedPixelSizeMicrons()
int W = server.getWidth(), H = server.getHeight()
int tilePx = (int) Math.round((P.TILE_UM as double) / pxUm)
if (tilePx > W || tilePx > H) throw new IllegalStateException("image ${W}x${H} px is smaller than a ${tilePx} px tile")
def chNames = server.getMetadata().getChannels().collect { it.getName() }
def outChannels = (P.CHANNELS as String).split(',').collect { it.trim() }
if (outChannels.any { it in ['IAPP', 'HLA-A'] }) throw new IllegalArgumentException('IAPP and HLA-A are never used')
def outIdx = outChannels.collect { c -> int i = chNames.indexOf(c); if (i < 0 || chNames.count { it == c } != 1) throw new IllegalStateException("channel ${c} not found exactly once"); i }
long seed = (P.SEED as long) >= 0 ? (P.SEED as long) : (imageName.hashCode() as long) & 0x7fffffffL
def rng = new Random(seed)
println "[bench] ${imageName}: tile ${tilePx} px (${P.TILE_UM} um), N random/object/exo = ${nRandom}/${nObject}/${nExo}, seed ${seed}"

// ---- output folders and key ----
def outDir = new File(P.OUT_DIR as String)
def tilesDir = new File(outDir, 'tiles'), predDir = new File(outDir, 'pred'), runsDir = new File(outDir, 'runs')
[tilesDir, predDir, runsDir].each { it.mkdirs() }
def keyFile = new File(outDir, 'key.tsv')
String KEY_HEADER = 'code\timage\tx\ty\tw\th\ttile_kind'
// all key.tsv reads/writes under an exclusive lock: several lanes may run this script at once
def withKeyLock = { Closure body ->
    def ch = FileChannel.open(Paths.get(outDir.getPath(), 'key.tsv.lock'), StandardOpenOption.CREATE, StandardOpenOption.WRITE)
    def lock = ch.lock()
    try { return body() } finally { lock.release(); ch.close() }
}
def readKey = { -> keyFile.exists() ? keyFile.readLines('UTF-8').findAll { it.trim() && !it.startsWith('code\t') }.collect { it.split('\t', -1) as List } : [] }
boolean skip = false
withKeyLock {
    def rows = readKey()
    def mine = rows.findAll { it[1] == imageName }
    if (!mine.isEmpty() && !(P.REPLACE as boolean)) {
        println "[bench] ${imageName} already has ${mine.size()} tiles in ${keyFile}; nothing to do (REPLACE=true regenerates)"
        skip = true
    } else if (!mine.isEmpty()) {
        mine.each { r -> new File(tilesDir, "${r[0]}.ome.tif").delete(); new File(predDir, "${r[0]}.geojson").delete() }
        keyFile.text = ([KEY_HEADER] + rows.findAll { it[1] != imageName }.collect { it.join('\t') }).join('\n') + '\n'
        println "[bench] REPLACE: removed ${mine.size()} previous tiles of ${imageName}"
    }
}

if (!skip) {
    // ---- tissue and detector objects ----
    def ENDO = ['Islet', 'EndocrineCluster', 'EndocrineSingle']
    def cls = { o -> o.getPathClass()?.toString() }
    def anns = getAnnotationObjects()
    def tissueAnns = anns.findAll { cls(it) == 'Tissue' && it.getROI()?.isArea() }
    if (tissueAnns.isEmpty()) throw new IllegalStateException('no Tissue annotation')
    def plane = tissueAnns[0].getROI().getImagePlane()
    def tissueGeom = UnaryUnionOp.union(tissueAnns.collect { it.getROI().getGeometry() })
    def tissuePrep = PreparedGeometryFactory.prepare(tissueGeom)
    def maskAnns = anns.findAll { cls(it) in ['Tissue', 'Artifact'] && it.getROI()?.isArea() }
    def objs = (P.OBJECTS_GEOJSON as String) ? PathIO.readObjects(Paths.get(P.OBJECTS_GEOJSON as String)).findAll { cls(it) in ENDO } :
            anns.findAll { cls(it) in ENDO }
    if (objs.isEmpty()) println "[bench] WARN: no detector objects (run 1_panc_objects_v3.groovy or set OBJECTS_GEOJSON); object tiles and preds will be empty"
    def route = { o -> double r = o.getMeasurementList().get('Endocrine route'); Double.isNaN(r) ? 0 : (int) r }
    def objTree = new STRtree()
    objs.each { objTree.insert(it.getROI().getGeometry().getEnvelopeInternal(), it) }
    objTree.build()
    println "[bench] ${objs.size()} detector objects " + objs.countBy { "${cls(it)}:r${route(it)}" }

    // ---- coarse maps for tissue coverage and exo_cd99 scoring: tissue, CD99 fold > tau, max hormone fold > tau ----
    def anchors = [:]
    def csv = new File(P.ANCHORS_CSV as String)
    if (nExo > 0) {
        if (!csv.exists()) throw new IOException("anchors CSV not found: ${csv}")
        def lines = csv.readLines('UTF-8')
        def hdr = lines[0].split(',', -1) as List
        int iImg = hdr.indexOf('image'), iTag = hdr.indexOf('tag'), iCh = hdr.indexOf('channel'), iA = hdr.indexOf('anchor')
        def rows = lines.drop(1).collect { it.split(',', -1) as List }   // no quoted fields in this CSV's key columns
        def mine = rows.findAll { it[iImg] == imageName }
        if (mine.isEmpty()) mine = rows.findAll { it[iTag] == tag }
        ['CD99', 'INS', 'GCG', 'SST'].each { ch ->
            def r = mine.find { it[iCh] == ch }
            if (r == null || !(r[iA] as double > 0)) throw new IllegalStateException("no anchor for ${imageName} ${ch}")
            anchors[ch] = r[iA] as double
        }
    }
    double ds = (P.SCORE_UM as double) / pxUm
    int Wc = (int) Math.ceil(W / ds), Hc = (int) Math.ceil(H / ds)
    // per-cell masks at SCORE_UM: tissue, CD99 fold > TAU_CD99, max hormone fold > TAU_H
    def tImg = new BufferedImage(Wc, Hc, BufferedImage.TYPE_BYTE_GRAY)
    def g2 = tImg.createGraphics()
    g2.setTransform(AffineTransform.getScaleInstance(1.0 / ds, 1.0 / ds))
    g2.setColor(java.awt.Color.WHITE)
    tissueAnns.each { g2.fill(RoiTools.getShape(it.getROI())) }
    g2.dispose()
    byte[] tBytes = ((java.awt.image.DataBufferByte) tImg.getRaster().getDataBuffer()).getData()
    byte[] cMask = new byte[Wc * Hc], hMask = new byte[Wc * Hc]
    long tS = System.currentTimeMillis()
    if (nExo > 0) {
        def op = ImageOps.buildImageDataOp(['CD99', 'INS', 'GCG', 'SST'].collect { ColorTransforms.createChannelExtractor(it) })
                .appendOps(ImageOps.Core.multiply(['CD99', 'INS', 'GCG', 'SST'].collect { 1.0d / (anchors[it] as double) } as double[]),
                           ImageOps.Core.splitMerge(ImageOps.Channels.extract(0),
                                   ImageOps.Core.sequential(ImageOps.Channels.extract(1, 2, 3), ImageOps.Channels.maximum())),
                           ImageOps.Core.ensureType(PixelType.FLOAT32))
        def privateServer = imageData.getServerBuilder().build()   // own reader pool, as in script 1
        try {
            privateServer.setMetadata(server.getMetadata())
            def privateData = new ImageData(privateServer, new PathObjectHierarchy(), imageData.getImageType())
            def calS = new PixelCalibration.Builder().pixelSizeMicrons(P.SCORE_UM as double, P.SCORE_UM as double).build()
            int TS = 1024   // score-scale pixels per block; blocks are read in parallel (one whole-image read is serial)
            def dsrv = ImageOps.buildServer(privateData, op, calS, TS, TS)
            double dsr = dsrv.getDownsampleForResolution(0)
            double tauC0 = P.TAU_CD99 as double, tauH0 = P.TAU_H as double
            def blocks = []
            for (int by = 0; by < Hc; by += TS) for (int bx = 0; bx < Wc; bx += TS) {
                def blk = GeometryTools.getDefaultFactory().toGeometry(new Envelope(bx * ds, Math.min(W, (bx + TS) * ds), by * ds, Math.min(H, (by + TS) * ds)))
                if (tissuePrep.intersects(blk)) blocks << [bx, by]
            }
            blocks.parallelStream().forEach { b ->
                int bx = b[0] as int, by = b[1] as int
                int x = (int) Math.round(bx * ds), y = (int) Math.round(by * ds)
                int w = Math.min(W - x, (int) Math.round(TS * ds)), h = Math.min(H - y, (int) Math.round(TS * ds))
                if (w <= 0 || h <= 0) return
                def r = dsrv.readRegion(RegionRequest.createInstance(dsrv.getPath(), dsr, x, y, w, h)).getRaster()
                int rw = Math.min(r.getWidth(), Wc - bx), rh = Math.min(r.getHeight(), Hc - by)
                float[] cr = new float[rw], hr = new float[rw]
                for (int yy = 0; yy < rh; yy++) {
                    r.getSamples(0, yy, rw, 1, 0, cr); r.getSamples(0, yy, rw, 1, 1, hr)
                    int base = (by + yy) * Wc + bx
                    for (int xx = 0; xx < rw; xx++) {
                        if (cr[xx] > tauC0) cMask[base + xx] = 1
                        if (hr[xx] > tauH0) hMask[base + xx] = 1
                    }
                }
            }
        } finally {
            privateServer.close()
        }
    }
    // summed-area tables (counts fit an int: <= Wc*Hc cells)
    int[] iT = new int[(Wc + 1) * (Hc + 1)], iC = new int[(Wc + 1) * (Hc + 1)], iH = new int[(Wc + 1) * (Hc + 1)]
    for (int y = 0; y < Hc; y++) {
        int rt = 0, rc = 0, rh = 0
        int base0 = y * (Wc + 1), base1 = (y + 1) * (Wc + 1), m = y * Wc
        for (int x = 0; x < Wc; x++) {
            boolean inT = tBytes[m + x] != 0
            if (inT) { rt++; if (cMask[m + x] != 0) rc++; if (hMask[m + x] != 0) rh++ }
            iT[base1 + x + 1] = iT[base0 + x + 1] + rt
            iC[base1 + x + 1] = iC[base0 + x + 1] + rc
            iH[base1 + x + 1] = iH[base0 + x + 1] + rh
        }
    }
    def boxSum = { int[] I, int x0, int y0, int x1, int y1 ->   // coarse cells [x0,x1) x [y0,y1)
        I[y1 * (Wc + 1) + x1] - I[y0 * (Wc + 1) + x1] - I[y1 * (Wc + 1) + x0] + I[y0 * (Wc + 1) + x0]
    }
    def coarse = { int x, int y ->
        int x0 = Math.max(0, (int) Math.floor(x / ds)), y0 = Math.max(0, (int) Math.floor(y / ds))
        int x1 = Math.min(Wc, (int) Math.ceil((x + tilePx) / ds)), y1 = Math.min(Hc, (int) Math.ceil((y + tilePx) / ds))
        [x0, y0, x1, y1, Math.max(1, (x1 - x0) * (y1 - y0))]
    }
    def tissueFrac = { int x, int y -> def c = coarse(x, y); boxSum(iT, c[0], c[1], c[2], c[3]) / (double) c[4] }
    println String.format("[bench] coarse maps %dx%d at %.1f um in %.1f s", Wc, Hc, P.SCORE_UM as double, (System.currentTimeMillis() - tS) / 1000.0)

    // ---- tile selection ----
    def chosen = []   // [x, y, kind, details]
    def clampXY = { double cx, double cy ->
        [(int) Math.max(0, Math.min(W - tilePx, Math.round(cx - tilePx / 2.0))), (int) Math.max(0, Math.min(H - tilePx, Math.round(cy - tilePx / 2.0)))]
    }
    def overlapOk = { int x, int y ->
        chosen.every { t ->
            double ox = Math.max(0, Math.min(x + tilePx, t[0] + tilePx) - Math.max(x, t[0]))
            double oy = Math.max(0, Math.min(y + tilePx, t[1] + tilePx) - Math.max(y, t[1]))
            ox * oy <= (P.MAX_OVERLAP_FRAC as double) * tilePx * tilePx
        }
    }
    // random inside Tissue
    def env = tissueGeom.getEnvelopeInternal()
    int tries = 0, nR = 0
    while (nR < nRandom && tries < 20000) {
        tries++
        double px = env.getMinX() + rng.nextDouble() * env.getWidth(), py = env.getMinY() + rng.nextDouble() * env.getHeight()
        if (!tissuePrep.contains(GeometryTools.getDefaultFactory().createPoint(new org.locationtech.jts.geom.Coordinate(px, py)))) continue
        def (int x, int y) = clampXY(px, py)
        double tf = tissueFrac(x, y)
        if (tf < (P.RANDOM_MIN_TISSUE as double) || !overlapOk(x, y)) continue
        chosen << [x, y, 'random', [tissue_frac: tf]]
        nR++
    }
    // centred on detector objects, round-robin over class x route strata
    def strata = new LinkedHashMap()
    ENDO.each { c -> [3, 2, 1].each { r -> def L = objs.findAll { cls(it) == c && route(it) == r }; if (L) { Collections.shuffle(L, rng); strata["object:${c}:r${r}".toString()] = L } } }
    int nO = 0
    def cursors = strata.collectEntries { k, v -> [k, 0] }
    while (nO < nObject && strata.any { k, v -> cursors[k] < v.size() }) {
        for (def k : strata.keySet()) {
            if (nO >= nObject) break
            def L = strata[k]
            while (cursors[k] < L.size()) {
                def o = L[cursors[k]]
                cursors[k] = cursors[k] + 1
                def (int x, int y) = clampXY(o.getROI().getCentroidX(), o.getROI().getCentroidY())
                if (!overlapOk(x, y)) continue
                chosen << [x, y, k, [object_area_um2: o.getMeasurementList().get("Area ${cal.getPixelWidthUnit()}^2".toString()), tissue_frac: tissueFrac(x, y)]]
                nO++
                break
            }
        }
    }
    // exocrine with scattered CD99-bright pixels and no hormone signal
    int nE = 0
    if (nExo > 0) {
        def cands = []
        int step = Math.max(1, (int) (tilePx / 4))
        for (int y = 0; y + tilePx <= H; y += step) for (int x = 0; x + tilePx <= W; x += step) {
            def c = coarse(x, y)
            double tf = boxSum(iT, c[0], c[1], c[2], c[3]) / (double) c[4]
            if (tf < (P.EXO_MIN_TISSUE as double)) continue
            double hf = boxSum(iH, c[0], c[1], c[2], c[3]) / (double) c[4]
            if (hf > (P.EXO_MAX_HORMONE_FRAC as double)) continue
            double cf = boxSum(iC, c[0], c[1], c[2], c[3]) / (double) c[4]
            if (cf <= 0) continue
            def box = GeometryTools.getDefaultFactory().toGeometry(new Envelope(x, x + tilePx, y, y + tilePx))
            boolean hasIslet = objTree.query(box.getEnvelopeInternal()).any { o -> cls(o) == 'Islet' && o.getROI().getGeometry().intersects(box) }
            if (hasIslet) continue
            cands << [x, y, cf, hf, tf]
        }
        cands.sort { a, b -> b[2] <=> a[2] }
        for (def c : cands) {
            if (nE >= nExo) break
            if (!overlapOk(c[0] as int, c[1] as int)) continue
            chosen << [c[0], c[1], 'exo_cd99', [cd99_pos_frac: c[2], hormone_pos_frac: c[3], tissue_frac: c[4]]]
            nE++
        }
        println "[bench] exo_cd99: ${cands.size()} eligible positions, best CD99+ fraction ${cands ? String.format('%.4f', cands[0][2]) : 'n/a'}"
    }
    println "[bench] selected ${nR} random, ${nO} object, ${nE} exo_cd99 tiles"

    // ---- codes (random, unique across key.tsv), shuffled so code order says nothing ----
    Collections.shuffle(chosen, rng)
    def sr = new SecureRandom()
    String ALPH = 'ABCDEFGHJKLMNPQRSTUVWXYZ23456789'
    def existing = withKeyLock { readKey().collect { it[0] } as Set }
    def newCode = { ->
        String c
        while (true) {
            c = 'T' + (1..8).collect { ALPH.charAt(sr.nextInt(ALPH.length())) }.join('')
            if (!existing.contains(c)) { existing.add(c); return c }
        }
    }

    // ---- write tiles, preds, key ----
    def gfac = GeometryTools.getDefaultFactory()
    def forbidden = [imageName, tag, stem].findAll { it && it.length() >= 6 }.unique()
    // ISO-8859-1 maps bytes 1:1 to chars, so String.indexOf is an exact (and fast) byte search
    def containsBytes = { String hay, byte[] needle -> hay.indexOf(new String(needle, 'ISO-8859-1')) >= 0 }
    def keyRows = []
    def details = []
    long tW = System.currentTimeMillis()
    def readServer = imageData.getServerBuilder().build()   // own reader pool, as in script 1
    readServer.setMetadata(server.getMetadata())
    try { chosen.each { t ->
        int x = t[0] as int, y = t[1] as int
        String code = newCode()
        String tilePath = new File(tilesDir, "${code}.ome.tif").getPath()
        new OMEPyramidWriter.Builder(readServer)
                .region(x, y, tilePx, tilePx)
                .channels(outIdx as int[])
                .name(code)
                .downsamples(1.0d)
                .tileSize(512)
                .losslessCompression()
                .parallelize()
                .build()
                .writePyramid(tilePath)
        String hay = new String(Files.readAllBytes(Paths.get(tilePath)), 'ISO-8859-1')
        def leak = forbidden.find { containsBytes(hay, it.getBytes('UTF-8')) || containsBytes(hay, it.getBytes('UTF-16LE')) }
        if (leak != null) { new File(tilePath).delete(); throw new IllegalStateException("tile ${code} contains the source name '${leak}'; blinding failed") }
        // detector objects (and Tissue/Artifact outlines) clipped to the tile, in tile coordinates
        def box = gfac.toGeometry(new Envelope(x, x + tilePx, y, y + tilePx))
        def shift = AffineTransformation.translationInstance(-x, -y)
        def preds = []
        def addClipped = { src, boolean withMeas ->
            def g = src.getROI().getGeometry()
            if (!g.intersects(box)) return
            def inter = g.intersection(box)
            def polys = PolygonExtracter.getPolygons(inter)
            if (polys.isEmpty()) return
            def gg = shift.transform(polys.size() == 1 ? polys[0] : gfac.buildGeometry(polys))
            def o = PathObjects.createAnnotationObject(GeometryTools.geometryToROI(gg, plane), src.getPathClass())
            if (withMeas) {
                def ml = o.getMeasurementList()
                src.getMeasurementList().getNames().each { n -> ml.put(n, src.getMeasurementList().get(n)) }
                ml.put('Clipped to tile', inter.getArea() < g.getArea() - 1e-6 ? 1 : 0)
                ml.close()
            }
            preds << o
        }
        objTree.query(box.getEnvelopeInternal()).each { addClipped(it, true) }
        maskAnns.each { addClipped(it, false) }
        String predPath = new File(predDir, "${code}.geojson").getPath()
        if (preds.isEmpty()) new File(predPath).text = '{"type":"FeatureCollection","features":[]}'
        else exportObjectsToGeoJson(preds, predPath, 'FEATURE_COLLECTION')
        keyRows << [code, imageName, x, y, tilePx, tilePx, t[2]].join('\t')
        details << [code: code, x: x, y: y, w: tilePx, h: tilePx, tile_kind: t[2], n_pred_objects: preds.count { cls(it) in ENDO }] + (t[3] as Map)
    } } finally {
        readServer.close()
    }
    withKeyLock {
        if (!keyFile.exists() || keyFile.length() == 0) keyFile.text = KEY_HEADER + '\n'
        keyFile.append(keyRows.join('\n') + '\n', 'UTF-8')
    }
    def plain
    plain = { v -> v instanceof Map ? v.collectEntries { k, x -> [String.valueOf(k), plain(x)] } : v instanceof Collection ? v.collect { plain(it) } :
            (v instanceof Double && (v.isNaN() || v.isInfinite())) ? null : (v instanceof CharSequence ? v.toString() : v) }
    new File(runsDir, "${tag}.json").text = GsonTools.getInstance(true).toJson(plain([image: imageName, tag: tag, params: P, seed: seed,
            tile_px: tilePx, n_objects: objs.size(), objects_source: (P.OBJECTS_GEOJSON ?: 'hierarchy'), tiles: details,
            written_at: java.time.OffsetDateTime.now().toString()]))
    println String.format("[bench] wrote %d tiles (%s) + preds + key rows in %.1f s -> %s", chosen.size(),
            chosen.countBy { it[2] }, (System.currentTimeMillis() - tW) / 1000.0, outDir)
}
