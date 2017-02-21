cls
set more off
clear all
*set processors 1

cap log close _all
log using areg_xtreg.log, replace
version

set obs 2000000
gen y = uniform()
gen x1 = uniform()
gen x2 = uniform()
gen x3 = uniform()
gen id = 1 + int((_n-1)/1000)
bys id: gen long t = _n
compress
sort id t
xtset id t

local vce vce(cluster id)

set rmsg on

reghdfe y in 1/10, noabsorb

****
cls
set trace off
cd ../test
*do setup
cd ../benchmark

reghdfe y x* , absorb(id) v(1) nosample
exit
****

timer clear

loc i 0
timeit `++i': areg y x* , absorb(id) `vce'
timeit `++i': xtreg y x*, fe `vce'
timeit `++i': reghdfe y x* , absorb(id) `vce' old // v3
timeit `++i': reghdfe y x* , absorb(id) `vce' old fast // v3-fast
timeit `++i': reghdfe y x* , absorb(id) `vce' old fast keepsing dof(none) // v3-fastest
timeit `++i': reghdfe y x* , absorb(id) `vce' // v4

timer list
log close _all

set rmsg off
exit
