* BENCHMARK WRT RES2FE
* PAPER: http://web.stanford.edu/group/fwolak/cgi-bin/sites/default/files/jem-2014-0008.pdf
* CODE: http://economics.mit.edu/files/8662
*  - This is the Example.do file)
*  - Note: Timers (tic/toc) from https://github.com/sergiocorreia/stata-misc
clear all
cls

cap log close _all
log using res2fe.log, replace
version

set matastrict off
do twowayreg.ado

// Warmup
sysuse auto
twowayset turn trunk
projvar price weight, p(w_)
reghdfe price, a(turn trunk)
cls
clear

*** 0) Preliminaries

forvalues lo = 3/3 {
di `lo'
forvalues wo = 2/2 {
di `wo'
foreach vars of numlist 2 10 {

di `vars'

loc long = 10^`lo'
loc wide = 10^`wo'
*loc vars = 2
loc lout = 0.1
loc reps = 1

loc toto = `long'*`wide'
set more off

forvalues rep = 1/`reps' {


*** 1) Generate Data
drop _all
set obs `toto'
** Variables
forvalues var = 1/`vars' {
	gen x`var'= rnormal(0)
	}
** Fixed Effects
* Indicators
gen hhid = floor((_n-1)/`wide')
gen ttid = _n-1-hhid*`wide'
** Drop a fraction of observations;
gen out= uniform()
sort out
drop if _n<`lout'*`toto'
* Effects
gen hhef = rnormal(0)
gen ttef = rnormal(0)
bysort hhid: replace hhef = hhef[1]
gen hid = 1
replace hid = hid[_n-1] + 1*(hhid[_n-1]~=hhid[_n]) if _n>1
bysort ttid: replace ttef = ttef[1]
gen tid = 1
replace tid = tid[_n-1] + 1*(ttid[_n-1]~=ttid[_n]) if _n>1


** Dependent Variable
gen y = hhef + ttef + rnormal(0)
forvalues var = 1/`vars' {
	qui replace y= y + x`var'
	}

sort hid tid

timer clear

*** 2) Run Our procedure
timer on 1
qui {
di "twowayset"
twowayset hid tid // try both and pick the fastest permutation
di "projvar"
projvar y x*, p(w_)
reg w_y w_x*, noc robust
drop w_*
}
timer off 1

* reghdfe
timer on 2
qui reghdfe y x*, vce(robust) absorb(hid tid) nosample dof(none) tol(1e-6) keepsingletons
timer off 2

timer list
timer clear

}
}
}
}

// results: for N=90k, |FE1|=1000, |FE2|=100, speed is .65 for reghdfe vs .86 for resfe

log close _all

!del last_*
exit
