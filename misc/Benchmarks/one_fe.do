rebuild_git reghdfe
clear all
set more off
cls

set obs 5000000


gen y = uniform()
gen x1 = uniform()
gen x2 = uniform()
gen x3 = uniform()
gen id = 1 + int((_n-1)/1000)

set rmsg on
*areg y x* , absorb(id)
*reghdfe y x* , absorb(id) fast old
*reghdfe y x* , absorb(id) keepsingletons fast timeit

sort id
*areg y x* , absorb(id)
*reghdfe y x* , absorb(id) fast old
reghdfe y x* , absorb(id) keepsingletons fast timeit

set rmsg off
exit



cls
set more off
reghdfe y x* , absorb(id) keepsingletons fast timeit v(2)

exit


reghdfe y x* , absorb(id) fast old v(4)
