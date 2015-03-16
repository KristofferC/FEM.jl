import FEM.getmat
import FEM.getvec

facts("FEM.Tools") do

context("FEM.Tools.FloatArrayPool") do

matpool = MatPool()
vecpool = VecPool()

A = getmat(2, 2, "J", matpool)
B = getmat(2, 2, "J", matpool)
C = getmat(2, 2, "K", matpool)
D = getmat(3, 2, "J", matpool)

@fact is(A, B) => true
@fact is(A, D) => false
@fact is(B, C) => false

a = getvec(4, "N", vecpool)
b = getvec(4, "N", vecpool)
c = getvec(6, "N", vecpool)
d = getvec(4, "K", vecpool)

@fact is(a, b) => true
@fact is(b, c) => false
@fact is(b, d) => false

function function_creating_matrix(matpool::MatPool)
    return getmat(2,2, "a", matpool)
end

f = getmat(2,2, "a", matpool)
v = function_creating_matrix(matpool)
@fact is(f, v) => true


end # context

end # facts