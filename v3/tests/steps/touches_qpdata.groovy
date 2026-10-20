// Test step (stands in for 'preflight'): something else writes data.qpdata while the step runs -- here its
// mtime moves back an hour, as a GUI save would move it -- and the step then reports a change to save.
// The footer must refuse to save over the changed file and write NO sentinel, so that the next plan finds
// data.qpdata different from the last record and blocks the image until --accept-external.
def f = getProjectEntry().getEntryPath().resolve('data.qpdata').toFile()
f.setLastModified(f.lastModified() - 3_600_000L)
println "[test] moved the mtime of ${f} back by one hour"
V3_RESULT = [ok: true, changed: true, message: 'pretends to have changed something', details: [:]]
