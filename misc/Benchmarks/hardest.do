* This is "the hardest case ever" for iterative methods
* In this case, just use a direct one (areg/res2fe/etc)

clear all
cls
set more off

local N 100

set obs `N'
gen id1 = 1 + int((_n-1)/2)
gen id2 = 1 + int((_n)/2)
gen n = 1 + (_n==1 | _n==_N)
expand n
drop n

sort id1 id2
li id1 id2, sep(0)

cap drop y
gen y = (_n==1) // (_n==_N)
gen x = uniform()

do twowayreg.ado
tic
di "twowayset"
twowayset id1 id2
di "projvar"
projvar y x, p(w_)
reg w*, noc robust
drop w_*
toc, report

tic
reghdfe y x, absorb(id1 id2) dof(none) tol(1e-10) fast
// v(3) accel(cg) transf(sym) keepsing clear // fast
toc, report


exit

* alternatives:
reg2hdfe y, id1(id1) id2(id2) tol(1e-8) // verbose // 5000 iters
exit


*hdfe  y, absorb(id1 id2) dof(none) tol(1e-10) v(3) accel(a) transf(kac) keepsing clear // fast
hdfe  y, absorb(id1 id2) dof(none) tol(1e-10) v(3) accel(cg) transf(sym) keepsing clear // fast

* CG+SYM requires 52 iterations (~100 effective iterations)
* CG+CIM 100
* SD+KAC 500
* A+KAC 7000+ (!!)


exit

