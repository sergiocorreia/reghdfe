/*
Also see:
- https://juliamath.github.io/IterativeSolvers.jl/dev/
- https://github.com/JuliaMath/IterativeSolvers.jl/tree/master/src
*/

clear all
cls
discard

include "groupreg.mata", adopath

set obs 4

mata:

// Test hypot()
hypot(3, 4)
hypot(1e-200,1e-200)
hypot(1e154,1e154)
hypot(1e-15, 1e-15)
hypot(1e306,1e306)
hypot(., 10)
hypot(10, .)

A1 = (1,2 \ 4, 6 \ 10, 0 \ -5 , 4)
x = (10 \ 12)
z = (1\2\3\0)
w1 = (1\1\1\1)
w2 = (1\2\3\1)

idx = st_addvar("long", ("x1", "x2"))
st_store(., idx, A1)
idx = st_addvar("long", ("y"))
st_store(., idx, z)
idx = st_addvar("long", ("w1"))
st_store(., idx, w1)
idx = st_addvar("long", ("w2"))
st_store(., idx, w2)

A2 = mock_matrix_init(A1)
A2.weights = w2

// Test simple properties
rows(A2)
A2.rows()
A2.cols()
A2.size()

// Test methods
y = A2.mult(x)
assert(A1*x == y)
y


y = A2.mult_transpose(z)
assert(cross(A1, w2, z) == y)

// Benchmarks
qrsolve(A1 :* sqrt(w2), z :* sqrt(w2))

// Test LSMR
solution = lsmr(A2, z, maxiter=10, ., ., ., verbose=2)
solution.alphas
solution.stop_code

end

exit

reg y x* [fw=w2], nocons
reg y x* [aw=w2], nocons
reg y x* [pw=w2], nocons


exit

r
rc
assert(0)

svsolve(A1, z)

cholsolve(cross(A1,A1), cross(A1,z))
lusolve(cross(A1,A1), cross(A1,z))
