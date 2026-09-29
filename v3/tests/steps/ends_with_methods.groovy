// Test step: a body that ends with a method definition and a class (the footer is appended after them).
// Stands in for 'preflight' (postcondition v3_result/preflight reads V3_RESULT).
def n = twice(21)
println "[test] body ran, twice(21) = ${n}, Box.v = ${new Box().v}"
V3_RESULT = [ok: n == 42, changed: false, message: 'method-ending body ran', details: [n: n]]

int twice(int x) { return 2 * x }

class Box { int v = 7 }
