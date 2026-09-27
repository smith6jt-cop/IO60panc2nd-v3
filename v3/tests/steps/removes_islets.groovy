// Test step: finishes normally but leaves no Islet annotations, so the footer's islets_present
// postcondition fails -- the footer must then NOT save (data.qpdata keeps its islets).
def islets = getAnnotationObjects().findAll { it.getPathClass()?.toString() == 'Islet' }
removeObjects(islets, true)
println "[test] removed ${islets.size()} islets in memory"
