/*
 * islet_benchmark_tiles.groovy -- blinded benchmark tiles for calibrating the v3 islet detector (design v2).
 *
 * Goal (user ruling 2026-09-28): tiles that include enough REPRESENTATIVE POSITIVE AND NEGATIVE selections in every
 * stratum = acquisition group x disease (E1..E4 x ND / Aab+ / T1D), not a tile count. Each tile is drawn for one
 * CATEGORY from a candidate pool, with a quota per category and stratum (QUOTA_PER_STRATUM, split over the
 * stratum's images, at least MIN_PER_IMAGE when the pool is not empty):
 *
 *   positives (object pools; a pool element is a candidate GROUP, see below)
 *     pos_islet_small / _medium / _large   Islet-sized (mean caliper >= 35 um): < 60, 60-120, >= 120 um
 *     pos_ins_rich                         Islet/cluster, INS+ area >= 50% of hormone+ area
 *     pos_ins_deficient                    Islet/cluster with hormone signal but INS+ < 5% of it (GCG/SST-dominant)
 *     pos_infiltrated                      Islet whose 20 um periphery is immune-dense (>= INFILTRATED_MIN of its
 *                                          coarse pixels immune-bright: CD3e / CD8 / CD20 above the image's own p99)
 *     pos_cluster, pos_single              EndocrineCluster (10-35 um), EndocrineSingle (<= 10 um)
 *     pos_cd99_only, pos_hormone_only      Endocrine route 1 (CD99 only) / route 2 (hormone only) candidates
 *     pos_legacy_only                      a legacy Islets2 islet that no v3 setting detected
 *   hard cases
 *     hard_near_threshold                  detected by some but not all of the v3 settings given in CANDIDATES
 *   negatives (region pools over a window grid; top TOP_FRAC of eligible windows by the rule's score)
 *     neg_exocrine                         random exocrine window (tissue >= EXO_MIN_TISSUE) without islet/cluster
 *                                          candidates of any source
 *     neg_exo_cd99                         ... with the most CD99+ single candidates (scattered CD99-bright cells)
 *     neg_speckle                          ... with the most hormone-only single candidates (SST speckle in panel 2,
 *                                          INS in the metabolic panel; the dominant channel is recorded)
 *     neg_duct / neg_vessel                Pan-Cytokeratin- / CD34-bright windows (islets allowed: hard negatives)
 *     neg_immune                           immune-dense windows without islet/cluster candidates (aggregates)
 *     neg_edge                             tissue fraction 0.3-0.7 (the Tissue outline crosses the tile)
 *     neg_fold                             "bright in every channel" (DAPI >= p99 and Pan-Cytokeratin, CD34, CD3e
 *                                          >= p90 of the image): fold-like
 *
 * Candidates are the UNION of several detector settings (CANDIDATES: label=path of 1_panc_objects_v3 outputs, e.g.
 * one sweep), the legacy Islets2 islets (LEGACY_GEOJSON), the structure rules above and random draws. Overlapping
 * candidates of different sources form one GROUP (so an islet found by 5 settings is one pool element); a group's
 * attributes come from the DEFAULT_LABEL object when there is one. Every tile records its pool size, the number
 * drawn and the inclusion probability (n drawn / pool size; pools are shuffled with a seeded RNG). Rare categories
 * are drawn first so common ones cannot crowd them out (tiles of one image overlap <= MAX_OVERLAP_FRAC).
 *
 * Blinding (identity must not leak to annotators):
 *   - OUT_DIR/annotate/tiles/<code>.ome.tif is the ONLY annotator-visible output: DAPI, CD99, INS, GCG, SST as
 *     UINT16 counts, fixed channel names and a FIXED display palette, pixel size written as PIXEL_SIZE_OUT (0.5078)
 *     for every image, OME name = code, file mtime set to one fixed date; each file is checked for the source
 *     name/tag/stem (bytes) before it is kept. v3/qupath/islet_benchmark_report.py scans every tile's metadata.
 *   - tile centres are JITTERED (target + U(-JITTER_FRAC, JITTER_FRAC) x tile, target kept >= 10% inside), so
 *     an object is never exactly centred; region windows are jittered by half a grid step.
 *   - OUT_DIR/private/ (never share): key.tsv (code -> image, x, y, w, h, tile_kind = category, stratum,
 *     inclusion probability...), runs/<tag>.json (parameters, pools), pred/<label>/<code>.geojson (each source's
 *     objects clipped to the tile, tile coordinates: the scorer's --pred-dir label=private/pred/<label>) and
 *     regions/<code>.geojson (Tissue/pLN/Artifact outlines; kept OUT of pred/, which the scorer reads whole).
 * Codes are random (SecureRandom), assigned after shuffling, and unique across key.tsv.
 *
 * Read-only: the image data is never modified or saved. An image already in key.tsv is skipped unless REPLACE=true
 * (which removes its old rows and files first). Images marked excluded in the strata table (6457) are skipped.
 *
 * Headless:  QuPath script -p <project> -i "<image>" [-a KEY=VALUE ...] islet_benchmark_tiles.groovy
 * Coverage (category x stratum, gaps) and the leak check:  islet_benchmark_report.py --out-dir OUT_DIR
 */

import qupath.lib.common.ColorTools
import qupath.lib.images.ImageData
import qupath.lib.images.servers.ColorTransforms
import qupath.lib.images.servers.ImageChannel
import qupath.lib.images.servers.ImageServerMetadata
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
import org.locationtech.jts.algorithm.MinimumBoundingCircle
import org.locationtech.jts.algorithm.MinimumDiameter
import org.locationtech.jts.geom.Coordinate
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
String OUT_DIR              = '/scratch/io60/benchmark_v2'   // annotate/ = tiles only (share this); private/ = the rest
double TILE_UM              = 400.0
int QUOTA_PER_STRATUM       = 3        // tiles per category per stratum (group x disease), split over its images
int MIN_PER_IMAGE           = 1        // ... at least this many per image when its pool is not empty
double MAX_OVERLAP_FRAC     = 0.25     // between tiles of this image
double JITTER_FRAC          = 0.3
double MIN_TISSUE           = 0.5      // tissue fraction of every tile (except neg_edge)
double EXO_MIN_TISSUE       = 0.9      // ... of the exocrine negatives
double SCORE_UM             = 4.0      // coarse maps for the structure rules
double BRIGHT_Q             = 0.99     // "bright" = above this tissue quantile of the image's own coarse channel
double FOLD_Q               = 0.90
double TOP_FRAC             = 0.10     // region pools: the top fraction of eligible windows by the rule's score
double INFILTRATED_MIN      = 0.10     // pos_infiltrated: immune-bright fraction of the islet's 20 um periphery
String CANDIDATES           = ''       // label=path[,label=path...] detector outputs (1_panc_objects_v3 WRITE_GEOJSON)
String DEFAULT_LABEL        = 'default'
String LEGACY_GEOJSON       = ''       // '' = results/v3/annotations_all/<tag>.geojson (its 'Islet' = legacy Islets2); 'none'
String STRATA_TSV           = ''       // '' = the built-in donor table below; else TSV: donor, group, disease[, excluded]
String DUCT_CHANNEL         = 'Pan-Cytokeratin'
String VESSEL_CHANNEL       = 'CD34'
String IMMUNE_CHANNELS      = 'CD3e,CD8,CD20'
String CHANNELS             = 'DAPI,CD99,INS,GCG,SST'   // never IAPP (failed QC) or HLA-A
double PIXEL_SIZE_OUT       = 0.5078   // written into every tile (no per-image calibration digits leak)
long SEED                   = -1L      // -1 = derived from the image name (reproducible)
boolean REPLACE             = false
// ====================================================

def P = new LinkedHashMap([OUT_DIR: OUT_DIR, TILE_UM: TILE_UM, QUOTA_PER_STRATUM: QUOTA_PER_STRATUM, MIN_PER_IMAGE: MIN_PER_IMAGE,
        MAX_OVERLAP_FRAC: MAX_OVERLAP_FRAC, JITTER_FRAC: JITTER_FRAC, MIN_TISSUE: MIN_TISSUE, EXO_MIN_TISSUE: EXO_MIN_TISSUE,
        SCORE_UM: SCORE_UM, BRIGHT_Q: BRIGHT_Q, FOLD_Q: FOLD_Q, TOP_FRAC: TOP_FRAC, INFILTRATED_MIN: INFILTRATED_MIN,
        CANDIDATES: CANDIDATES, DEFAULT_LABEL: DEFAULT_LABEL, LEGACY_GEOJSON: LEGACY_GEOJSON, STRATA_TSV: STRATA_TSV,
        DUCT_CHANNEL: DUCT_CHANNEL, VESSEL_CHANNEL: VESSEL_CHANNEL, IMMUNE_CHANNELS: IMMUNE_CHANNELS, CHANNELS: CHANNELS,
        PIXEL_SIZE_OUT: PIXEL_SIZE_OUT, SEED: SEED, REPLACE: REPLACE])
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

// ---- strata: acquisition group x disease (user ruling 2026-09-28; the qptiff Responsivity date is NOT a batch) ----
def STRATA = [
        '6476': ['E1', 'ND'], 'HDL115': ['E1', 'ND'], '6521': ['E1', 'Aab+'], '6579': ['E1', 'T1D'], '6380': ['E1', 'T1D'], '6551': ['E1', 'T1D'],
        '6482': ['E2', 'ND'], '6374': ['E2', 'ND'], '6493': ['E2', 'ND'], '6539': ['E2', 'ND'], '6442': ['E2', 'ND'], '6461': ['E2', 'ND'],
        '6593': ['E2', 'T1D'], '6436': ['E2', 'T1D'], '6457': ['E2', 'T1D', 'excluded'],
        '6505': ['E3', 'Aab+'], '6450': ['E3', 'Aab+'], '6623': ['E3', 'T1D'], '6534': ['E3', 'T1D'], '6566': ['E3', 'T1D'], '6591': ['E3', 'T1D'], '6523': ['E3', 'T1D'],
        'HDL112': ['E4', 'ND'], '6347': ['E4', 'Aab+'], '6533': ['E4', 'T1D'], '6580': ['E4', 'T1D']]
if (P.STRATA_TSV) {
    STRATA = [:]
    new File(P.STRATA_TSV as String).readLines('UTF-8').findAll { it.trim() && !it.startsWith('#') && !it.startsWith('donor\t') }.each { l ->
        def f = l.split('\t', -1); STRATA[f[0].trim()] = [f[1].trim(), f[2].trim()] + ((f.size() > 3 && f[3].trim() in ['1', 'true', 'excluded']) ? ['excluded'] : [])
    }
}

def imageData = getCurrentImageData()
def server = imageData.getServer()
def entry = getProjectEntry()
String imageName = entry?.getImageName() ?: server.getMetadata().getName()
String tag = imageName.replaceAll(/[^A-Za-z0-9._-]+/, '_')
String stem = imageName.replaceAll(/(?i)\.(er\.)?qptiff.*$/, '').replaceAll(/(?i)\.ome\.tif+.*$/, '')
def dm = (imageName =~ /(HDL\d+)|(?<!\d)(\d{4})(?!\d)/)
String donor = dm.find() ? (dm.group(1) ?: dm.group(2)) : null
def st = STRATA[donor]
if (st == null) throw new IllegalStateException("donor '${donor}' of ${imageName} is not in the strata table (STRATA_TSV)")
String group = st[0], disease = st[1], stratum = "${st[0]}_${st[1]}".toString()
boolean excludedImg = st.size() > 2 && st[2] == 'excluded'
int nInStratum = STRATA.count { k, v -> v[0] == group && v[1] == disease && !(v.size() > 2 && v[2] == 'excluded') }
int nPerImage = Math.max(P.MIN_PER_IMAGE as int, (int) Math.ceil((P.QUOTA_PER_STRATUM as int) / (double) Math.max(1, nInStratum)))

def cal = server.getPixelCalibration()
double pxUm = cal.getAveragedPixelSizeMicrons()
int W = server.getWidth(), H = server.getHeight()
int tilePx = (int) Math.round((P.TILE_UM as double) / pxUm)
if (tilePx > W || tilePx > H) throw new IllegalStateException("image ${W}x${H} px is smaller than a ${tilePx} px tile")
def chNames = server.getMetadata().getChannels().collect { it.getName() }
def outChannels = (P.CHANNELS as String).split(',').collect { it.trim() }
if (outChannels.any { it in ['IAPP', 'HLA-A'] }) throw new IllegalArgumentException('IAPP and HLA-A are never used')
def idxOf = { String c -> int i = chNames.indexOf(c); if (i < 0 || chNames.count { it == c } != 1) throw new IllegalStateException("channel ${c} not found exactly once"); i }
def outIdx = outChannels.collect { idxOf(it) }
long seed = (P.SEED as long) >= 0 ? (P.SEED as long) : (imageName.hashCode() as long) & 0x7fffffffL
def rng = new Random(seed)
println "[bench] ${imageName}: donor ${donor}, stratum ${stratum} (${nInStratum} image(s)) -> ${nPerImage} tile(s) per category; " +
        "tile ${tilePx} px (${P.TILE_UM} um), seed ${seed}"

// ---- output folders and key: annotate/ holds tiles ONLY; everything that identifies a tile lives in private/ ----
def outDir = new File(P.OUT_DIR as String)
def annotateDir = new File(outDir, 'annotate'), privDir = new File(outDir, 'private')
def tilesDir = new File(annotateDir, 'tiles'), predRoot = new File(privDir, 'pred'), regionsDir = new File(privDir, 'regions'), runsDir = new File(privDir, 'runs')
[tilesDir, predRoot, regionsDir, runsDir].each { it.mkdirs() }
def keyFile = new File(privDir, 'key.tsv')
def KEY_COLS = ['code', 'image', 'x', 'y', 'w', 'h', 'tile_kind', 'pixel_size_um', 'donor', 'group', 'disease', 'stratum', 'target',
                'target_sources', 'target_x', 'target_y', 'pool_size', 'n_drawn', 'inclusion_prob', 'pool_rule', 'seed']
String KEY_HEADER = KEY_COLS.join('\t')
def withKeyLock = { Closure body ->   // several lanes may run this script at once
    def ch = FileChannel.open(Paths.get(privDir.getPath(), 'key.tsv.lock'), StandardOpenOption.CREATE, StandardOpenOption.WRITE)
    def lock = ch.lock()
    try { return body() } finally { lock.release(); ch.close() }
}
def readKey = { -> keyFile.exists() ? keyFile.readLines('UTF-8').findAll { it.trim() && !it.startsWith('code\t') }.collect { it.split('\t', -1) as List } : [] }
boolean skip = excludedImg
if (excludedImg) println "[bench] ${imageName} is excluded from analysis (strata table); no tiles"
withKeyLock {
    if (skip) return
    def rows = readKey()
    def mine = rows.findAll { it[1] == imageName }
    if (!mine.isEmpty() && !(P.REPLACE as boolean)) {
        println "[bench] ${imageName} already has ${mine.size()} tiles in ${keyFile}; nothing to do (REPLACE=true regenerates)"
        skip = true
    } else if (!mine.isEmpty()) {
        mine.each { r -> new File(tilesDir, "${r[0]}.ome.tif").delete(); new File(regionsDir, "${r[0]}.geojson").delete()
            predRoot.listFiles()?.each { d -> new File(d, "${r[0]}.geojson").delete() } }
        keyFile.text = ([KEY_HEADER] + rows.findAll { it[1] != imageName }.collect { it.join('\t') }).join('\n') + '\n'
        println "[bench] REPLACE: removed ${mine.size()} previous tiles of ${imageName}"
    }
}

if (!skip) {
    long t0 = System.currentTimeMillis()
    def gfac = GeometryTools.getDefaultFactory()
    def ENDO = ['Islet', 'EndocrineCluster', 'EndocrineSingle']
    def cls = { o -> o.getPathClass()?.toString() }
    String UM = cal.getPixelWidthUnit()
    // ---- analysis area: Tissue minus pLN and Artifact ----
    def anns = getAnnotationObjects()
    def unionOf = { List<String> cs -> def gs = anns.findAll { cls(it) in cs && it.getROI()?.isArea() }.collect { it.getROI().getGeometry() }; gs ? UnaryUnionOp.union(gs) : null }
    def tissueU = unionOf(['Tissue'])
    if (tissueU == null) throw new IllegalStateException('no Tissue annotation')
    def exclU = unionOf(['pLN', 'Artifact'])
    def areaGeom = exclU == null ? tissueU : tissueU.difference(exclU)
    def plane = anns.find { cls(it) == 'Tissue' }.getROI().getImagePlane()
    def regionAnns = anns.findAll { cls(it) in ['Tissue', 'pLN', 'Artifact'] && it.getROI()?.isArea() }

    // ---- candidates: detector settings + legacy Islets2, grouped by overlap ----
    def meanCaliperUm = { g -> ((2.0 * new MinimumBoundingCircle(g).getRadius()) + new MinimumDiameter(g).getLength()) / 2.0 * pxUm }
    def classOf = { double d -> d >= 35 ? 'Islet' : (d > 10 ? 'EndocrineCluster' : 'EndocrineSingle') }
    def mv = { o, String k -> double v = o.getMeasurementList().get(k); v }
    def sources = new LinkedHashMap()   // label -> objects
    (P.CANDIDATES as String).split(',').collect { it.trim() }.findAll { it }.each { spec ->
        int eq = spec.indexOf('=')
        if (eq <= 0) throw new IllegalArgumentException("CANDIDATES entry '${spec}' is not label=path")
        String label = spec.substring(0, eq).trim(), path = spec.substring(eq + 1).trim()
        if (label == 'legacy') throw new IllegalArgumentException("the label 'legacy' is reserved")
        sources[label] = new File(path).exists() ? PathIO.readObjects(Paths.get(path)).findAll { cls(it) in ENDO } : []
        if (!new File(path).exists()) println "[bench] WARN: ${path} does not exist (no objects: an empty detector output)"
    }
    String legacyPath = (P.LEGACY_GEOJSON as String) ?: "/home/smith6jt/IO60panc2nd/results/v3/annotations_all/${tag}.geojson".toString()
    if (legacyPath != 'none') {
        if (new File(legacyPath).exists()) sources['legacy'] = PathIO.readObjects(Paths.get(legacyPath)).findAll { cls(it) == 'Islet' }
        else println "[bench] WARN: no legacy Islets2 export ${legacyPath}"
    }
    if (!sources.containsKey(P.DEFAULT_LABEL as String)) println "[bench] WARN: DEFAULT_LABEL '${P.DEFAULT_LABEL}' is not among the candidate sources ${sources.keySet()}"
    def v3Labels = sources.keySet().findAll { it != 'legacy' }
    def cands = []   // [obj, label, geom]
    sources.each { label, objs -> objs.each { o -> def g = o.getROI().getGeometry(); if (!g.isEmpty() && g.intersects(areaGeom)) cands << [obj: o, label: label, geom: g] } }
    def cTree = new STRtree()
    cands.eachWithIndex { c, i -> cTree.insert(c.geom.getEnvelopeInternal(), i) }
    cTree.build()
    int[] parent = new int[cands.size()]
    for (int i = 0; i < parent.length; i++) parent[i] = i
    def find = { int i -> while (parent[i] != i) { parent[i] = parent[parent[i]]; i = parent[i] }; i }
    cands.eachWithIndex { c, i ->
        def prep = PreparedGeometryFactory.prepare(c.geom)
        cTree.query(c.geom.getEnvelopeInternal()).each { j -> int jj = j as int
            if (jj > i && prep.intersects(cands[jj].geom)) { int a = find(i), b = find(jj); if (a != b) parent[a] = b } }
    }
    def groupsIdx = [:].withDefault { [] }
    cands.eachWithIndex { c, i -> groupsIdx[find(i)] << i }
    def groups = groupsIdx.values().collect { List idx ->
        def members = idx.collect { cands[it] }
        def pick = { List ms -> ms.max { it.geom.getArea() } }
        def rep = pick(members.findAll { it.label == P.DEFAULT_LABEL }) ?: pick(members.findAll { it.label != 'legacy' }) ?: pick(members)
        def o = rep.obj
        double d = mv(o, "Mean caliper diameter ${UM}".toString())
        if (Double.isNaN(d)) d = meanCaliperUm(rep.geom)
        def c = rep.geom.getCentroid()
        [rep: rep, labels: members.collect { it.label }.toSet(), cls: classOf(d), diam: d,
         route: rep.label == 'legacy' ? 0 : (Double.isNaN(mv(o, 'Endocrine route')) ? 0 : (int) mv(o, 'Endocrine route')),
         ins: mv(o, "INS+ area ${UM}^2".toString()), gcg: mv(o, "GCG+ area ${UM}^2".toString()), sst: mv(o, "SST+ area ${UM}^2".toString()),
         hormone: mv(o, "Hormone+ area ${UM}^2".toString()), cx: c.getX(), cy: c.getY(), env: rep.geom.getEnvelopeInternal()]
    }
    def gTree = new STRtree()
    groups.each { gTree.insert(new Envelope(it.env), it) }
    gTree.build()
    println "[bench] candidates: " + sources.collectEntries { k, v -> [k, v.size()] } + " -> ${groups.size()} groups " + groups.countBy { it.cls }

    // ---- coarse maps at SCORE_UM: analysis area, DAPI, duct, vessel, immune (max of channel / own p99) ----
    double ds = (P.SCORE_UM as double) / pxUm
    int Wc = (int) Math.ceil(W / ds), Hc = (int) Math.ceil(H / ds)
    def tImg = new BufferedImage(Wc, Hc, BufferedImage.TYPE_BYTE_GRAY)
    def g2 = tImg.createGraphics()
    g2.setTransform(AffineTransform.getScaleInstance(1.0 / ds, 1.0 / ds))
    g2.setColor(java.awt.Color.WHITE)
    g2.fill(GeometryTools.geometryToShape(areaGeom))
    g2.dispose()
    byte[] tBytes = ((java.awt.image.DataBufferByte) tImg.getRaster().getDataBuffer()).getData()
    def immuneNames = (P.IMMUNE_CHANNELS as String).split(',').collect { it.trim() }.findAll { it && chNames.contains(it) }
    def mapNames = ['DAPI', P.DUCT_CHANNEL as String, P.VESSEL_CHANNEL as String] + immuneNames
    mapNames.each { idxOf(it) }
    int nMap = mapNames.size()
    float[][] maps = new float[nMap][Wc * Hc]
    long tS = System.currentTimeMillis()
    def op = ImageOps.buildImageDataOp(mapNames.collect { ColorTransforms.createChannelExtractor(it) }).appendOps(ImageOps.Core.ensureType(PixelType.FLOAT32))
    def privateServer = imageData.getServerBuilder().build()   // own reader pool, as in script 1
    try {
        privateServer.setMetadata(server.getMetadata())
        def privateData = new ImageData(privateServer, new PathObjectHierarchy(), imageData.getImageType())
        def calS = new PixelCalibration.Builder().pixelSizeMicrons(P.SCORE_UM as double, P.SCORE_UM as double).build()
        int TS = 1024
        def dsrv = ImageOps.buildServer(privateData, op, calS, TS, TS)
        double dsr = dsrv.getDownsampleForResolution(0)
        def areaPrep = PreparedGeometryFactory.prepare(areaGeom)
        def blocks = []
        for (int by = 0; by < Hc; by += TS) for (int bx = 0; bx < Wc; bx += TS) {
            def blk = gfac.toGeometry(new Envelope(bx * ds, Math.min(W, (bx + TS) * ds), by * ds, Math.min(H, (by + TS) * ds)))
            if (areaPrep.intersects(blk)) blocks << [bx, by]
        }
        blocks.parallelStream().forEach { b ->
            int bx = b[0] as int, by = b[1] as int
            int x = (int) Math.round(bx * ds), y = (int) Math.round(by * ds)
            int w = Math.min(W - x, (int) Math.round(TS * ds)), h = Math.min(H - y, (int) Math.round(TS * ds))
            if (w <= 0 || h <= 0) return
            def r = dsrv.readRegion(RegionRequest.createInstance(dsrv.getPath(), dsr, x, y, w, h)).getRaster()
            int rw = Math.min(r.getWidth(), Wc - bx), rh = Math.min(r.getHeight(), Hc - by)
            float[] row = new float[rw]
            for (int c = 0; c < nMap; c++) for (int yy = 0; yy < rh; yy++) {
                r.getSamples(0, yy, rw, 1, c, row)
                System.arraycopy(row, 0, maps[c], (by + yy) * Wc + bx, rw)
            }
        }
    } finally {
        privateServer.close()
    }
    // per-channel tissue quantiles (the image's own scale; tile selection only, never used for calling)
    int stepQ = Math.max(1, (int) ((Wc * (long) Hc) / 2_000_000))
    int nQ = 0
    for (int i = 0; i < Wc * Hc; i += stepQ) if (tBytes[i] != 0) nQ++
    int[] tIdx = new int[nQ]
    nQ = 0
    for (int i = 0; i < Wc * Hc; i += stepQ) if (tBytes[i] != 0) tIdx[nQ++] = i
    def sortedCache = [:]
    def quant = { int c, double q ->
        double[] v = sortedCache[c]
        if (v == null) { v = new double[tIdx.length]; for (int i = 0; i < tIdx.length; i++) v[i] = maps[c][tIdx[i]]; Arrays.sort(v); sortedCache[c] = v }
        v.length == 0 ? Double.NaN : v[(int) Math.min(v.length - 1, Math.floor(q * v.length))]
    }
    double dapiHi = quant(0, P.BRIGHT_Q as double), dapiQ = quant(0, P.FOLD_Q as double)
    double ductHi = quant(1, P.BRIGHT_Q as double), ductQ = quant(1, P.FOLD_Q as double)
    double vesHi = quant(2, P.BRIGHT_Q as double), vesQ = quant(2, P.FOLD_Q as double)
    def immHi = (3..<nMap).collect { quant(it, P.BRIGHT_Q as double) }, immQ = (3..<nMap).collect { quant(it, P.FOLD_Q as double) }
    // summed-area tables over tissue cells: tissue, duct-bright, vessel-bright, immune-bright, fold-like
    int SW = Wc + 1
    int[] iT = new int[SW * (Hc + 1)], iD = new int[SW * (Hc + 1)], iV = new int[SW * (Hc + 1)], iI = new int[SW * (Hc + 1)], iF = new int[SW * (Hc + 1)]
    byte[] immBright = new byte[Wc * Hc]
    for (int y = 0; y < Hc; y++) {
        int rt = 0, rd = 0, rv = 0, ri = 0, rf = 0
        int b0 = y * SW, b1 = (y + 1) * SW, m = y * Wc
        for (int x = 0; x < Wc; x++) {
            int i = m + x
            if (tBytes[i] != 0) {
                rt++
                if (maps[1][i] >= ductHi) rd++
                if (maps[2][i] >= vesHi) rv++
                boolean imm = false, immF = immQ.isEmpty()
                for (int c = 3; c < nMap; c++) { if (maps[c][i] >= immHi[c - 3]) imm = true; if (maps[c][i] >= immQ[c - 3]) immF = true }
                if (imm) { ri++; immBright[i] = 1 }
                if (maps[0][i] >= dapiHi && maps[1][i] >= ductQ && maps[2][i] >= vesQ && immF) rf++
            }
            iT[b1 + x + 1] = iT[b0 + x + 1] + rt; iD[b1 + x + 1] = iD[b0 + x + 1] + rd; iV[b1 + x + 1] = iV[b0 + x + 1] + rv
            iI[b1 + x + 1] = iI[b0 + x + 1] + ri; iF[b1 + x + 1] = iF[b0 + x + 1] + rf
        }
    }
    maps = null
    def boxSum = { int[] I, int x0, int y0, int x1, int y1 -> I[y1 * SW + x1] - I[y0 * SW + x1] - I[y1 * SW + x0] + I[y0 * SW + x0] }
    def coarse = { int x, int y ->
        int x0 = Math.max(0, (int) Math.floor(x / ds)), y0 = Math.max(0, (int) Math.floor(y / ds))
        int x1 = Math.min(Wc, (int) Math.ceil((x + tilePx) / ds)), y1 = Math.min(Hc, (int) Math.ceil((y + tilePx) / ds))
        [x0, y0, x1, y1, Math.max(1, (x1 - x0) * (y1 - y0))]
    }
    def frac = { int[] I, int x, int y -> def c = coarse(x, y); boxSum(I, c[0], c[1], c[2], c[3]) / (double) c[4] }
    println String.format("[bench] coarse maps %dx%d at %.1f um (%s) in %.1f s", Wc, Hc, P.SCORE_UM as double, mapNames, (System.currentTimeMillis() - tS) / 1000.0)

    // immune density of each Islet group's 20 um periphery (coarse cells in the ring)
    double ringPx = 20.0 / pxUm
    groups.findAll { it.cls == 'Islet' }.each { gr ->
        def ring = gr.rep.geom.buffer(ringPx).difference(gr.rep.geom)
        def e = ring.getEnvelopeInternal()
        def prep = PreparedGeometryFactory.prepare(ring)
        int n = 0, k = 0
        for (int cy = Math.max(0, (int) (e.getMinY() / ds)); cy <= Math.min(Hc - 1, (int) (e.getMaxY() / ds)); cy++)
            for (int cx = Math.max(0, (int) (e.getMinX() / ds)); cx <= Math.min(Wc - 1, (int) (e.getMaxX() / ds)); cx++) {
                int i = cy * Wc + cx
                if (tBytes[i] == 0 || !prep.contains(gfac.createPoint(new Coordinate((cx + 0.5) * ds, (cy + 0.5) * ds)))) continue
                n++; if (immBright[i] != 0) k++
            }
        gr.immuneRing = n > 0 ? k / (double) n : 0.0
    }

    // ---- category pools ----
    def hormoneShare = { gr, double v -> (Double.isNaN(v) || Double.isNaN(gr.hormone as double) || (gr.hormone as double) <= 0) ? Double.NaN : v / (gr.hormone as double) }
    def insShareIs = { gr, Closure test -> double sh = hormoneShare(gr, gr.ins as double); !Double.isNaN(sh) && test(sh) }
    def objCats = new LinkedHashMap([
            pos_islet_small    : { gr -> gr.cls == 'Islet' && gr.diam < 60 },
            pos_islet_medium   : { gr -> gr.cls == 'Islet' && gr.diam >= 60 && gr.diam < 120 },
            pos_islet_large    : { gr -> gr.cls == 'Islet' && gr.diam >= 120 },
            // NaN (no hormone area: CD99-only or legacy-only groups) must match NEITHER: in dynamic Groovy
            // `NaN >= 0.5` is true, which put CD99-only objects into pos_ins_rich (review 2026-09-29)
            pos_ins_rich       : { gr -> gr.cls in ['Islet', 'EndocrineCluster'] && insShareIs(gr, { double v -> v >= 0.5 }) },
            pos_ins_deficient  : { gr -> gr.cls in ['Islet', 'EndocrineCluster'] && insShareIs(gr, { double v -> v < 0.05 }) },
            pos_infiltrated    : { gr -> gr.cls == 'Islet' && (gr.immuneRing ?: 0) >= (P.INFILTRATED_MIN as double) },
            pos_cluster        : { gr -> gr.cls == 'EndocrineCluster' },
            pos_single         : { gr -> gr.cls == 'EndocrineSingle' },
            pos_cd99_only      : { gr -> gr.route == 1 },
            pos_hormone_only   : { gr -> gr.route == 2 },
            pos_legacy_only    : { gr -> gr.labels == (['legacy'] as Set) },
            hard_near_threshold: { gr -> def v3 = gr.labels.findAll { it != 'legacy' }; !v3.isEmpty() && v3.size() < v3Labels.size() }])
    def pools = new LinkedHashMap()
    objCats.each { name, pred -> pools[name] = [kind: 'object', rule: 'all matching candidate groups', items: groups.findAll { pred(it) }] }
    // windows on a grid (step = tile / 4) inside the image
    int step = Math.max(1, (int) (tilePx / 4))
    def wins = []
    for (int y = 0; y + tilePx <= H; y += step) for (int x = 0; x + tilePx <= W; x += step) {
        double tf = frac(iT, x, y)
        if (tf < 0.3) continue
        def box = new Envelope(x, x + tilePx, y, y + tilePx)
        int nIC = 0, nCd = 0, nHo = 0
        def hoCh = [INS: 0, GCG: 0, SST: 0]
        gTree.query(box).each { gr ->
            if (gr.cls in ['Islet', 'EndocrineCluster']) { if (gr.env.intersects(box)) nIC++ }
            else if (gr.cx >= x && gr.cx < x + tilePx && gr.cy >= y && gr.cy < y + tilePx) {
                if (gr.route in [1, 3]) nCd++
                if (gr.route == 2) { nHo++; def top = [INS: gr.ins, GCG: gr.gcg, SST: gr.sst].max { e -> Double.isNaN(e.value) ? -1 : e.value }; hoCh[top.key]++ }
            }
        }
        wins << [x: x, y: y, tf: tf, nIC: nIC, nCd: nCd, nHo: nHo, hoCh: hoCh, duct: frac(iD, x, y), ves: frac(iV, x, y),
                 imm: frac(iI, x, y), fold: frac(iF, x, y)]
    }
    def topPool = { List elig, Closure score -> // top TOP_FRAC of eligible windows by score (ties kept)
        if (elig.isEmpty()) return []
        def s = elig.sort(false) { a, b -> score(b) <=> score(a) }
        int n = Math.max(1, (int) Math.ceil((P.TOP_FRAC as double) * s.size()))
        double cut = score(s[n - 1])
        s.findAll { score(it) >= cut && score(it) > 0 }
    }
    double exoT = P.EXO_MIN_TISSUE as double, minT = P.MIN_TISSUE as double
    // these negatives must hold no islet/cluster candidate AFTER the tile is jittered, not only on the grid
    // (review 2026-09-29: 4 of 26 jittered negative tiles contained one)
    def NO_IC_CATS = ['neg_exocrine', 'neg_exo_cd99', 'neg_speckle', 'neg_immune']
    def hitsIC = { double x0, double y0 ->
        def b = new Envelope(x0, x0 + tilePx, y0, y0 + tilePx)
        gTree.query(b).any { gr -> gr.cls in ['Islet', 'EndocrineCluster'] && gr.env.intersects(b) }
    }
    def exo = wins.findAll { it.tf >= exoT && it.nIC == 0 }
    pools.neg_exocrine = [kind: 'window', rule: "tissue >= ${exoT}, no islet/cluster candidate", items: exo]
    pools.neg_exo_cd99 = [kind: 'window', rule: "exocrine, top ${P.TOP_FRAC} by CD99+ single candidates", items: topPool(exo.findAll { it.nCd > 0 }, { it.nCd })]
    pools.neg_speckle = [kind: 'window', rule: "exocrine, top ${P.TOP_FRAC} by hormone-only single candidates", items: topPool(exo.findAll { it.nHo > 0 }, { it.nHo })]
    pools.neg_duct = [kind: 'window', rule: "tissue >= ${minT}, top ${P.TOP_FRAC} by ${P.DUCT_CHANNEL} >= p${P.BRIGHT_Q}", items: topPool(wins.findAll { it.tf >= minT && it.duct > 0 }, { it.duct })]
    pools.neg_vessel = [kind: 'window', rule: "tissue >= ${minT}, top ${P.TOP_FRAC} by ${P.VESSEL_CHANNEL} >= p${P.BRIGHT_Q}", items: topPool(wins.findAll { it.tf >= minT && it.ves > 0 }, { it.ves })]
    pools.neg_immune = [kind: 'window', rule: "tissue >= ${minT}, no islet/cluster, top ${P.TOP_FRAC} by ${immuneNames} >= p${P.BRIGHT_Q}", items: topPool(wins.findAll { it.tf >= minT && it.nIC == 0 && it.imm > 0 }, { it.imm })]
    pools.neg_edge = [kind: 'window', rule: 'tissue fraction 0.3-0.7', items: wins.findAll { it.tf >= 0.3 && it.tf <= 0.7 }]
    pools.neg_fold = [kind: 'window', rule: "tissue >= ${minT}, top ${P.TOP_FRAC} by fold-like (DAPI >= p${P.BRIGHT_Q}, others >= p${P.FOLD_Q})", items: topPool(wins.findAll { it.tf >= minT && it.fold > 0 }, { it.fold })]

    // ---- draw: rare categories first; jittered positions; overlap limit ----
    def chosen = []   // [x, y, category, details]
    def overlapOk = { int x, int y ->
        chosen.every { t ->
            double ox = Math.max(0, Math.min(x + tilePx, t[0] + tilePx) - Math.max(x, t[0]))
            double oy = Math.max(0, Math.min(y + tilePx, t[1] + tilePx) - Math.max(y, t[1]))
            ox * oy <= (P.MAX_OVERLAP_FRAC as double) * tilePx * tilePx
        }
    }
    double jit = (P.JITTER_FRAC as double) * tilePx, inset = 0.1 * tilePx
    def clampI = { double v, int lo, int hi -> (int) Math.max(lo, Math.min(hi, Math.round(v))) }
    def poolStats = new LinkedHashMap()
    pools.keySet().sort { a, b -> pools[a].items.size() <=> pools[b].items.size() ?: a <=> b }.each { String cat ->
        def pool = pools[cat]
        def items = new ArrayList(pool.items)
        Collections.shuffle(items, rng)
        int got = 0, rejected = 0
        for (def item : items) {
            if (got >= nPerImage) break
            int x, y
            String target = '', srcs = ''
            double tx = Double.NaN, ty = Double.NaN
            if (pool.kind == 'object') {
                tx = item.cx; ty = item.cy
                double cx = tx + (rng.nextDouble() * 2 - 1) * jit, cy = ty + (rng.nextDouble() * 2 - 1) * jit
                x = clampI(cx - tilePx / 2.0, 0, W - tilePx); y = clampI(cy - tilePx / 2.0, 0, H - tilePx)
                // keep the target inside the tile (>= 10% from its edges), even after clamping at the image edge
                x = clampI(Math.min(Math.max(x, tx + inset - tilePx), tx - inset), 0, W - tilePx)
                y = clampI(Math.min(Math.max(y, ty + inset - tilePx), ty - inset), 0, H - tilePx)
                target = "${item.cls}:r${item.route}:d${String.format(Locale.ROOT, '%.1f', item.diam as double)}".toString()
                srcs = item.labels.sort().join('+')
            } else {
                x = clampI(item.x + (rng.nextDouble() - 0.5) * step, 0, W - tilePx); y = clampI(item.y + (rng.nextDouble() - 0.5) * step, 0, H - tilePx)
                target = cat == 'neg_speckle' ? "speckle:${item.hoCh.max { e -> e.value }.key}".toString() : 'window'
            }
            double tf = frac(iT, x, y)
            boolean tissueOk = cat == 'neg_edge' ? (tf >= 0.3 && tf <= 0.7) : tf >= (cat in ['neg_exocrine', 'neg_exo_cd99', 'neg_speckle'] ? exoT : minT)
            if (!tissueOk || !overlapOk(x, y) || (cat in NO_IC_CATS && hitsIC(x, y))) { rejected++; continue }
            chosen << [x, y, cat, [target: target, target_sources: srcs, target_x: tx, target_y: ty, tissue_frac: tf]]
            got++
        }
        poolStats[cat] = [kind: pool.kind, rule: pool.rule, pool_size: items.size(), n_target: items.isEmpty() ? 0 : nPerImage, n_drawn: got,
                          n_rejected: rejected, inclusion_prob: items.isEmpty() ? 0.0 : got / (double) items.size()]
        chosen.findAll { it[2] == cat }.each { it[3].pool_size = items.size(); it[3].n_drawn = got; it[3].inclusion_prob = got / (double) items.size(); it[3].pool_rule = pool.rule }
    }
    println "[bench] pools (size -> drawn): " + poolStats.collect { k, v -> "${k} ${v.pool_size}->${v.n_drawn}" }.join(', ')

    // ---- codes (random, unique across key.tsv), shuffled so code order says nothing ----
    Collections.shuffle(chosen, rng)
    def sr = new SecureRandom()
    String ALPH = 'ABCDEFGHJKLMNPQRSTUVWXYZ23456789'
    def existing = withKeyLock { readKey().collect { it[0] } as Set }
    def newCode = { ->
        while (true) {
            String c = 'T' + (1..8).collect { ALPH.charAt(sr.nextInt(ALPH.length())) }.join('')
            if (!existing.contains(c)) { existing.add(c); return c }
        }
    }

    // ---- tile server: fixed channel names, FIXED palette, fixed pixel size, neutral name (private instance) ----
    def PALETTE = [DAPI: ColorTools.packRGB(0, 0, 255), CD99: ColorTools.packRGB(255, 255, 255), INS: ColorTools.packRGB(255, 0, 0),
                   GCG: ColorTools.packRGB(0, 255, 0), SST: ColorTools.packRGB(255, 255, 0)]
    def fixedChannels = chNames.withIndex().collect { String n, int i ->
        int j = outIdx.indexOf(i)
        j >= 0 ? ImageChannel.getInstance(outChannels[j], PALETTE[outChannels[j]] ?: ColorTools.packRGB(200, 200, 200)) :
                ImageChannel.getInstance("C${i}".toString(), ColorTools.packRGB(128, 128, 128))
    }
    def readServer = imageData.getServerBuilder().build()   // own reader pool, as in script 1
    readServer.setMetadata(new ImageServerMetadata.Builder(server.getMetadata()).name('tile').channels(fixedChannels)
            .pixelSizeMicrons(P.PIXEL_SIZE_OUT as double, P.PIXEL_SIZE_OUT as double).build())
    long FIXED_MTIME = java.time.OffsetDateTime.parse('2026-01-01T00:00:00Z').toInstant().toEpochMilli()
    def forbidden = [imageName, tag, stem].findAll { it && it.length() >= 6 }.unique()
    def containsBytes = { String hay, byte[] needle -> hay.indexOf(new String(needle, 'ISO-8859-1')) >= 0 }
    def keyRows = [], details = []
    def labels = new ArrayList(sources.keySet())
    labels.each { new File(predRoot, it).mkdirs() }
    def clipTo = { src, box, shift, boolean withMeas ->
        def g = src.getROI().getGeometry()
        if (!g.intersects(box)) return null
        def inter = g.intersection(box)
        def polys = PolygonExtracter.getPolygons(inter)
        if (polys.isEmpty()) return null
        def gg = shift.transform(polys.size() == 1 ? polys[0] : gfac.buildGeometry(polys))
        def o = PathObjects.createAnnotationObject(GeometryTools.geometryToROI(gg, plane), src.getPathClass())
        if (withMeas) {
            def ml = o.getMeasurementList()
            src.getMeasurementList().getNames().each { n -> ml.put(n, src.getMeasurementList().get(n)) }
            ml.put('Clipped to tile', inter.getArea() < g.getArea() - 1e-6 ? 1 : 0)
            ml.close()
        }
        o
    }
    def srcTrees = sources.collectEntries { lab, objs -> def t = new STRtree(); objs.each { t.insert(it.getROI().getGeometry().getEnvelopeInternal(), it) }; t.build(); [lab, t] }
    long tW = System.currentTimeMillis()
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
        new File(tilePath).setLastModified(FIXED_MTIME)
        def box = gfac.toGeometry(new Envelope(x, x + tilePx, y, y + tilePx))
        def shift = AffineTransformation.translationInstance(-x, -y)
        labels.each { lab ->
            def preds = srcTrees[lab].query(box.getEnvelopeInternal()).collect { clipTo(it, box, shift, true) }.findAll { it != null }
            String pp = new File(new File(predRoot, lab), "${code}.geojson").getPath()
            if (preds.isEmpty()) new File(pp).text = '{"type":"FeatureCollection","features":[]}'
            else exportObjectsToGeoJson(preds, pp, 'FEATURE_COLLECTION')
        }
        def regs = regionAnns.collect { clipTo(it, box, shift, false) }.findAll { it != null }
        String rp = new File(regionsDir, "${code}.geojson").getPath()
        if (regs.isEmpty()) new File(rp).text = '{"type":"FeatureCollection","features":[]}'
        else exportObjectsToGeoJson(regs, rp, 'FEATURE_COLLECTION')
        def d = t[3] as Map
        def fmtD = { v -> (v instanceof Double && (v.isNaN() || v.isInfinite())) ? '' : (v instanceof Double ? String.format(Locale.ROOT, '%.6g', v) : String.valueOf(v)) }
        def row = [code: code, image: imageName, x: x, y: y, w: tilePx, h: tilePx, tile_kind: t[2], pixel_size_um: fmtD(pxUm), donor: donor,
                   group: group, disease: disease, stratum: stratum, target: d.target, target_sources: d.target_sources,
                   target_x: fmtD(d.target_x), target_y: fmtD(d.target_y), pool_size: d.pool_size, n_drawn: d.n_drawn,
                   inclusion_prob: fmtD(d.inclusion_prob), pool_rule: d.pool_rule, seed: seed]
        keyRows << KEY_COLS.collect { String.valueOf(row[it] ?: (row[it] == 0 ? 0 : '')).replace('\t', ' ') }.join('\t')
        details << [code: code] + row + [tissue_frac: d.tissue_frac]
    } } finally {
        readServer.close()
    }
    withKeyLock {
        if (!keyFile.exists() || keyFile.length() == 0) keyFile.text = KEY_HEADER + '\n'
        keyFile.append(keyRows.join('\n') + (keyRows ? '\n' : ''), 'UTF-8')
    }
    def plain
    plain = { v -> v instanceof Map ? v.collectEntries { k, x -> [String.valueOf(k), plain(x)] } : v instanceof Collection ? v.collect { plain(it) } :
            (v instanceof Double && (v.isNaN() || v.isInfinite())) ? null : (v instanceof CharSequence ? v.toString() : v) }
    new File(runsDir, "${tag}.json").text = GsonTools.getInstance(true).toJson(plain([image: imageName, tag: tag, donor: donor, group: group,
            disease: disease, stratum: stratum, n_images_in_stratum: nInStratum, n_per_category: nPerImage, params: P, seed: seed,
            tile_px: tilePx, sources: sources.collectEntries { k, v -> [k, v.size()] }, n_groups: groups.size(),
            groups_by_class: groups.countBy { it.cls }, pools: poolStats, n_windows: wins.size(),
            quantiles: [dapi_hi: dapiHi, duct_hi: ductHi, vessel_hi: vesHi, immune_hi: immHi], tiles: details,
            seconds: (System.currentTimeMillis() - t0) / 1000.0, written_at: java.time.OffsetDateTime.now().toString()]))
    println String.format("[bench] wrote %d tiles (%s) + preds (%s) + key rows in %.1f s (total %.1f s) -> %s", chosen.size(),
            new TreeMap(chosen.countBy { it[2] }), labels, (System.currentTimeMillis() - tW) / 1000.0, (System.currentTimeMillis() - t0) / 1000.0, outDir)
}
