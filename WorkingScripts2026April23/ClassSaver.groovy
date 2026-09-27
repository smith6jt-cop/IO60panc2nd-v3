//def targetClass = getPathClass("Islet")
//
//// Save backup to file and strip
//def backupFile = new File(buildFilePath(PROJECT_BASE_DIR, "class_backup.csv"))
//backupFile.text = ""
//int i = 0
//for (def anno in getAnnotationObjects()) {
//    def pc = anno.getPathClass()
//    if (pc != null && pc != targetClass) {
//        backupFile << "${i},${pc}\n"
//        anno.setPathClass(null)
//    }
//    i++
//}
//fireHierarchyUpdate()
//print "Stripped. Run Signed Distance now."

def annotations = getAnnotationObjects().toList()
def backupFile = new File(buildFilePath(PROJECT_BASE_DIR, "class_backup.csv"))
for (def line in backupFile.readLines()) {
    def parts = line.split(",", 2)
    int idx = Integer.parseInt(parts[0])
    def cls = getPathClass(parts[1])
    annotations[idx].setPathClass(cls)
}
fireHierarchyUpdate()
print "Restored"