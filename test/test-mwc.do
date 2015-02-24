* Test multi-way clustering with -regress- and -ivreg2-
* Maybe we can also do it with -ivregress-

* Note: We don-t benchmark against -cgmreg- but against -fixed_cgmreg-
* This is due to a bug in the former, where they count K=rows(X'X)
* But the rows are also included for omitted variables (e.g. the base variable in a factor expansion i.turn)

cd "D:/Github/reghdfe/source"
cscript "reghdfe with multi-way clustering" adofile reghdfe

* Setup
	cap cls
	discard
	clear all
	set more off
	qui adopath + "D:/Github/reghdfe/test"
	set more off

* Convenience
	capture program drop TrimMatrix
	program define TrimMatrix, eclass
	args size
		assert `size'>0
		matrix trim_b = e(b)
		matrix trim_V = e(V)
		matrix trim_b = trim_b[1, 1..`size']
		matrix trim_V = trim_V[1..`size',1..`size']
		ereturn matrix trim_b = trim_b
		ereturn matrix trim_V = trim_V
	end
	
* Dataset
	sysuse auto
	drop if missing(rep)

* [TEST] Verify that it gives the same results as -areg- with one-way-clustering
	noi di as text " - two-way ols"
	areg price weight disp, absorb(rep) cluster(turn)
	TrimMatrix 2
	storedresults save benchmark e()
	
	reghdfe price weight disp, a(rep) tol(1e-10) vce(cluster turn)
	TrimMatrix 2
	storedresults compare benchmark e(), tol(1e-10) include(scalar: N r2 matrix: trim_b trim_V macros: wexp wtype clustvar)
	storedresults drop benchmark
	
* [TEST] Verify that it gives the same results as -cgmreg- with one-way-clustering
	noi di as text " - two-way ols"
	fixed_cgmreg price weight disp i.rep, cluster(turn)
	TrimMatrix 2
	storedresults save benchmark e()
	
	reghdfe price weight disp, a(rep) tol(1e-10) vce(cluster turn)
	TrimMatrix 2
	storedresults compare benchmark e(), tol(1e-10) include(scalar: N r2 matrix: trim_b trim_V macros: wexp wtype clustvar)
	storedresults drop benchmark

exit
* [TEST] Verify that the e() unique to reghdfe works
	assert e(N_clustvars)==1
	
exit

* [TEST] Verify that it gives the same results as -cgmreg- with two-way-clustering
	noi di as text " - two-way ols"
	fixed_cgmreg price weight disp i.rep, cluster(turn foreign)
	TrimMatrix 2
	storedresults save benchmark e()
	
	reghdfe price weight disp, a(rep) tol(1e-10) vce(cluster turn foreign)
	assert e(N_clustvars)==1
	TrimMatrix 2
	storedresults compare benchmark e(), tol(1e-10) include(scalar: N r2 N_clustvars matrix: trim_b trim_V macros: wexp wtype clustvar)
	storedresults drop benchmark
	
	
exit
* [TEST] Verify that it gives the same results as -cgmreg- when 
* [TEST] Verify that it gives the same results as -ivreg2- 
* [TEST] Nested within clusters
	
	
* [TEST] [fw] should fail with non-integer weights
	gen float m = 1 + uniform()
	rcof "reghdfe price weight disp [fw=m], a(turn#foreign)" == 401 // may not use noninteger freq weights
	drop m
	
* [TEST] Freq. weight with huge N (to test numeric stability)
	noi di as text " - freq. weight with large N"
	cap drop n
	gen double n = 100000
	
	areg price weight disp turn#foreign [fw=n], a(turn)
	TrimMatrix 2
	storedresults save areg e()
	
	reghdfe price weight disp [fw=n], a(turn#foreign) tol(1e-10)
	TrimMatrix 2
	storedresults compare areg e(), tol(1e-10) include(scalar: N r2 matrix: trim_b trim_V macros: wexp wtype)
	storedresults drop areg
	
* [TEST] Analitic weights
	noi di as text " - analitic weight"
	cap drop n
	gen n = 1 + (10 * uniform())
	reg price weight disp trunk turn##foreign [aw=n] // , a(turn)
	TrimMatrix 3
	storedresults save areg e()
	
	reghdfe price weight disp trunk [aw=n], a(turn#foreign) tol(1e-10) dofmethod(bounds)
	TrimMatrix 3
	storedresults compare areg e(), tol(1e-10) include(scalar: N r2 matrix: trim_b trim_V macros: wexp wtype)
	storedresults drop areg

* [TEST] [aw] should fail with negative weights
	gen float m = 1
	replace m = -2 in 1/2
	rcof "reghdfe price weight disp [aw=m], a(turn#foreign)" == 402 // negative weights
	drop m
	
* [TEST] Probability/sampling weights
	noi di as text " - probability weight"
	cap drop n
	gen n = 1 + (10 * uniform())
	reg price weight disp trunk turn##foreign [pw=n] // , a(turn)
	TrimMatrix 3
	storedresults save areg e()
	
	reghdfe price weight disp trunk [pw=n], a(turn#foreign) tol(1e-10) dofmethod(bounds)
	TrimMatrix 3
	storedresults compare areg e(), tol(1e-10) include(scalar: N r2 matrix: trim_b trim_V macros: wexp wtype)
	storedresults drop areg
	
	
	
	exit
/*
clear all
discard
cls
set more off

sysuse auto
keep if rep<.

replace turn = int(turn/14)
collapse (count) n=price (mean) price, by(turn rep foreign)

reghdfe price foreign [aw=n] , a(i.turn##c.rep)
local b = _b[foreign]
reg price foreign i.turn i.turn#c.rep [aw=n]
assert abs(`b'-_b[foreign])<0.01


exit
*/

clear all
discard
cls
set more off
sysuse auto

gen n = 100000 // 0000 // int(uniform()*10+3)
reghdfe price weight [fw=n], a(turn##c.length)
exit
areg price weight turn#c.length [fw=n], a(turn)





* THIS IS THE SIMPLEST ILLUSTRATION OF THE BUG

clear all
discard
cls
set more off

sysuse auto
keep if rep<.

collapse (count) n=price (mean) price, by(foreign turn)

gen c = 1
reghdfe price foreign [fw=n] , a(turn turn#c.turn)
reghdfe price foreign [fw=n] , a(turn##c.turn)
reg price foreign i.turn##c.turn [fw=n]



discard
clear all
set obs 10
gen c = 1
set seed 32453245
gen y = uniform()
gen x = uniform()
gen v = uniform()
gen n = 1 + int(10*uniform())
reghdfe y x [fw=n] , a(i.c i.c#c.v)
reg y x v [fw=n]





