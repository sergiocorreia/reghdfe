cd "D:/Github/reghdfe" // /source
cscript "reghdfe with weights" adofile reghdfe

* Setup
	discard
	clear all
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
	
* Create fake dataset
	sysuse auto
	gen n = int(uniform()*10+3)
	
* [TEST] Freq. weight
	noi di as text " - freq. weight"
	qui areg price weight disp turn#foreign [fw=n], a(turn)
	TrimMatrix 2
	storedresults save areg e()

	reghdfe price weight disp [fw=n], a(turn#foreign) tol(1e-10)
	TrimMatrix 2
	storedresults compare areg e(), tol(1e-10) include(scalar: N r2 matrix: trim_b trim_V macros: wexp wtype)

foreach suite in default avar {
	di as text "SUITE=<`suite'>"
	reghdfe price weight disp [fw=n], a(turn#foreign) tol(1e-10) vce(ols, suite(`suite'))
	TrimMatrix 2
	storedresults compare areg e(), tol(1e-10) include(scalar: N r2 matrix: trim_b trim_V macros: wexp wtype)	
}

	storedresults drop areg

* [TEST] [fw] should fail with non-integer weights
	gen float m = 1 + uniform()
	rcof "reghdfe price weight disp [fw=m], a(turn#foreign)" == 401 // may not use noninteger freq weights
	drop m
	
* [TEST] Freq. weight with huge N (to test numeric stability)
	noi di as text " - freq. weight with large N"
	cap drop n
	gen double n = 100000
	
	qui areg price weight disp turn#foreign [fw=n], a(turn)
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
	qui reg price weight disp trunk turn##foreign [aw=n] , robust
	TrimMatrix 3
	storedresults save areg e()
	
	reghdfe price weight disp trunk [aw=n], a(turn#foreign) tol(1e-10) dof(none) vce(robust)
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
	qui reg price weight disp trunk turn##foreign [pw=n] , robust
	TrimMatrix 3
	storedresults save areg e()
	
	reghdfe price weight disp trunk [pw=n], a(turn#foreign) tol(1e-10) dof(none) vce(robust)
	TrimMatrix 3
	storedresults compare areg e(), tol(1e-10) include(scalar: N r2 matrix: trim_b trim_V macros: wexp wtype)
	storedresults drop areg
	
	
cd "D:/Github/reghdfe/test"
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





