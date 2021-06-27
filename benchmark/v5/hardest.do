* This is "the hardest case ever" for iterative methods
* In this case, we need to apply pruning, or use a direct method (areg/res2fe/etc)

clear all
cls
set more off

// Reload
reghdfe, reload
pr drop _all

set matastrict off
do twowayreg.ado

// Warmup
sysuse auto
twowayset turn trunk
projvar price weight, p(w_)
reghdfe price, a(turn trunk)
cls


// -------------------------------------------------------------------------------------------------
// Zig-zag dataset
// -------------------------------------------------------------------------------------------------

clear
local N 2000 // avoid more than 1,000 obs as twowayreg hangs
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

// -------------------------------------------------------------------------------------------------
// Benchmark: Somaini-Wolak
// -------------------------------------------------------------------------------------------------

timer clear
timer on 1
di "twowayset"
twowayset id1 id2
di "projvar"
projvar y x, p(w_)
reg w*, noc robust
drop w_*
timer off 1

timer on 2
reghdfe y x, absorb(id1 id2) dof(none)
timer off 2

timer on 3
reghdfe y x, absorb(id1 id2) dof(none) noprune
timer off 3

timer on 4
reghdfe y x, absorb(id1 id2) dof(none) noprune accel(cg) transf(cim)
timer off 4

timer on 5
reghdfe y x, absorb(id1 id2) dof(none) noprune accel(sd) transf(kac)
timer off 5

timer list
timer clear

!del last_*

exit

/*
1) projvar
2) reghdfe prune
3) reghdfe noprune
4) reghdfe cg+cim noprune
5) reghdfe sd+kac noprune

. timer list
   1:      5.08 /        1 =       5.0780
   2:      0.13 /        1 =       0.1250
   3:      0.68 /        1 =       0.6830
   4:      1.09 /        1 =       1.0940
   5:      1.25 /        1 =       1.2500
*/

* CG+SYM requires N/2 iterations (N effective iterations)
* CG+CIM requires N iterations (2x the amount of work)
* SD+KAC same (?)
* A+KAC does not converge (acceleration is too aggressive!)
