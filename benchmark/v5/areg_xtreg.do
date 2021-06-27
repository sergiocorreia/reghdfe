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

timer clear

loc i 0
timeit `++i': areg y x* , absorb(id) `vce'
timeit `++i': xtreg y x*, fe `vce'
timeit `++i': reghdfe y x* , absorb(id) `vce' old // v3
timeit `++i': reghdfe y x* , absorb(id) `vce' old fast // v3-fast
timeit `++i': reghdfe y x* , absorb(id) `vce' old fast keepsing dof(none) // v3-fastest
timeit `++i': reghdfe y x* , absorb(id) `vce' // v4
timeit `++i': reghdfe y x* , absorb(id) `vce' nosample keepsing // v4

timer list

// TLDR for a 2mm sample:
// areg			2.6s
// xtreg		6.2s
// old reghdfe	2.4-3s
// new reghdfe	1.6-1.7s

// 2.4/2.6 = 92% runtime of areg before
// 1.6/2.6 = 62% now (!)

log close _all

set rmsg off
exit
