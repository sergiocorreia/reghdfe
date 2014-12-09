clear all
set more off
pr drop _all
discard
cap cls

* Many variables
set varabbrev off

set obs 1000
gen fe = _n<500
gen y = uniform()

forv i=1/100 {
	gen xxxxxxxxxx`i' = uniform()
}

forv i=1/50 {
	gen zzzz`i' = uniform()
	gen wwww`i' = uniform()
}

reghdfe y x*, a(fe)
ds x*
local xs = r(varlist)
reghdfe y `xs', a(fe)

set trace off
set tracedepth 4
reghdfe y xxxxxxxxxx1 (zzzz1=wwww1), a(fe) verbose(3) nested
reghdfe y x* (z*=w*), a(fe) verbose(3)
