* BENCHMARK WRT RES2FE
* PAPER: http://web.stanford.edu/group/fwolak/cgi-bin/sites/default/files/jem-2014-0008.pdf
* CODE: http://economics.mit.edu/files/8662
*  - This is the Example.do file)
*  - Note: Timers (tic/toc) from https://github.com/sergiocorreia/stata-misc

which tic.ado
cap log close _all
log using res2fe.log, replace
version

do twowayreg.ado
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

*** 2) Run Our procedure
tic
di "twowayset"
twowayset hid tid
di "projvar"
projvar y x*, p(w_)
reg w_y w_x*, noc robust
drop w_*
toc, report

* Old and Slow
tic
reghdfe y x*, vce(robust) absorb(tid hid) old
toc, report

* Slow
tic
reghdfe y x*, vce(robust) absorb(tid hid)
toc, report

* Fast
tic
reghdfe y x*, vce(robust) absorb(tid hid) fast dof(none) tol(1e-6) keepsingletons group(20) // v(3) timeit 
toc, report

}
}
}
}

log close _all
