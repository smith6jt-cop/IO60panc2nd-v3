// Test step: throws after changing the hierarchy in memory; nothing may be saved.
removeObjects(getAnnotationObjects().findAll { it.getPathClass()?.toString() == 'Islet' }, true)
throw new IllegalStateException('[test] deliberate failure after an in-memory change')
