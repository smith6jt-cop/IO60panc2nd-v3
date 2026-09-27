// Test step: renames one channel and reports success. saveImageData() then rewrites
// data/<id>/server.json, but QuPath reloads channel names from project.qpproj -- the driver must
// report this as drift and block the image. Stands in for 'preflight'.
def names = getCurrentServer().getMetadata().getChannels().collect { it.getName() }
names[1] = names[1] + '_renamed'
setChannelNames(*names as String[])
V3_RESULT = [ok: true, changed: true, message: "renamed channel 2 to ${names[1]}".toString(), details: [:]]
println "[test] ${V3_RESULT.message}"
