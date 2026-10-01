/*
 * djl_cuda_check.groovy -- fast fail for GPU jobs: can QuPath's DJL PyTorch engine use CUDA here?
 *
 * HiPerGator kills a job whose GPU sits idle for an hour, and InstanSeg quietly falling back to the
 * CPU would look exactly like that. So c0_gpu.sbatch runs this first (seconds) and stops unless the
 * sentinel says ok. No image, no project, no download (the engine must already be in ~/.djl.ai).
 *
 *   QuPath script -a <sentinel.json> djl_cuda_check.groovy
 */

if (args == null || args.length < 1) throw new IllegalArgumentException('usage: djl_cuda_check.groovy <sentinel.json>')
def out = new File(args[0])
def result = [ok: false, time: java.time.OffsetDateTime.now().toString(), host: InetAddress.getLocalHost().getHostName(),
              cuda_visible_devices: System.getenv('CUDA_VISIBLE_DEVICES')]
try {
    def engine = qupath.ext.djl.DjlTools.getEngine('PyTorch', false)   // false: never download on a compute node
    result.engine_version = engine.getVersion()
    result.gpu_count = engine.getGpuCount()
    if (result.gpu_count > 0) {
        // a real CUDA op, not just a device count: 1024x1024 ones @ ones sums to 1024^3
        def manager = ai.djl.ndarray.NDManager.newBaseManager(ai.djl.Device.gpu(0), 'PyTorch')
        try {
            long t0 = System.nanoTime()
            def a = manager.ones(new ai.djl.ndarray.types.Shape(1024, 1024))
            double s = a.matMul(a).sum().toType(ai.djl.ndarray.types.DataType.FLOAT64, false).getDouble()
            result.matmul_sum = s
            result.matmul_ms = (System.nanoTime() - t0) / 1e6
            result.device = a.getDevice().toString()
            result.ok = Math.abs(s - Math.pow(1024, 3)) < 1 && a.getDevice().isGpu()
        } finally {
            manager.close()
        }
    } else {
        result.error = 'PyTorch engine loaded but sees no CUDA device'
    }
} catch (Throwable t) {
    result.error = t.toString()
}
// InstanSeg reads pixels through QuPath's OpenCV natives. On HiPerGator they load only with a newer
// libstdc++ on LD_LIBRARY_PATH (hpg_env.sh); without it step 0 dies after the model loads (2026-10-01).
try {
    def cv = Class.forName('org.bytedeco.opencv.global.opencv_core')
    org.bytedeco.javacpp.Loader.load(cv)
    result.opencv_version = cv.getField('CV_VERSION').get(null)
    result.opencv_ok = true
} catch (Throwable t) {
    result.opencv_ok = false
    result.opencv_error = t.toString()
    result.ok = false
}
out.getParentFile()?.mkdirs()
out.text = new com.google.gson.GsonBuilder().setPrettyPrinting().serializeNulls().create().toJson(result)
println "[djl_cuda_check] ${result.ok ? 'OK' : 'FAILED'}: ${result}"
