// Test step: leaves early with a top-level `return`, as the real scripts do on a missing precondition.
// The appended footer must never run, so no sentinel is written and the driver records a failure.
println '[test] body starts'
if (getAnnotationObjects().size() > 0) {
    println '[test] early top-level return'
    return
}
println '[test] not reached'
