clear all
cls
set more off
set type double


loc n = 1e6

set obs `n'
gen y = rnormal()
forv i=1/5 {
	gen x`i' = rnormal()
}

gen long i = 1 + int((_n-1)/10)
gen long j = int(ceil(runiform()*1e5))
gdistinct i j


set rmsg on
reghdfe y x*, a(i j) dof(none) nosample v(1) pool(3) // prune()


* 0:		6.8 6.9 7.0
* 1: 		9.8 9.9
* 2:		9.0
* 3:		7.9
