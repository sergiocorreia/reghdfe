* Test that it doesn't crash
set varabbrev on
cd "D:/Github/reghdfe/source"
*cscript "minimal test" adofile reghdfe

* Setup
	cap cls
	discard
	clear all
	set more off
	qui adopath + "D:/Github/reghdfe/test"
	set more off


set trace off
pr drop _all
	
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

	replace length = 0 if rep==3
	replace length = 5 if rep==1

bys foreign: gen t = _n
tsset foreign t

*reghdfe price weight (disp=gear), a(rep#foreign##c.length turn) vce(cluster foreign t, bw(2) kernel(tru))
set varabbrev off
set trace on
reghdfe price weight (displacement=gear_ratio), a(rep78#foreign##c.length turn) vce(, kiefer) verbose(3)
set varabbrev on
asd
reghdfe price weight (disp=gear), a(rep#foreign##c.length turn) vce(cluster foreign t, bw(2) kernel(tru))
reghdfe price weight (disp=gear), a(rep#foreign##c.length turn) vce(cluster foreign t, bw(2) kernel(tru))

asd
	
*reghdfe price weight disp, a(foreign#rep turn#rep) tol(1e-10) vce(cluster rep#foreign) verbose(3) dof(all)	
*reghdfe price weight disp, a(rep#i.turn##c.length foreign) tol(1e-10) vce(cluster turn) verbose(3) dof(all)
reghdfe price weight (disp=gear), a(rep#foreign##c.length turn) ///
	vce(cluster turn foreign#head, bw(2) kernel(tru)) ///
	tol(1e-10) verbose(3) dof(all)
	
ahora_con_bw_y_avar_y_default_y_todo
reghdfe price weight disp, a(rep foreign) tol(1e-10) vce(cluster turn#rep turn) verbose(3) dof(all)

asd
	
* [TEST] Verify that it gives the same results as -areg- with one-way-clustering
	reghdfe price weight disp, a(rep) tol(1e-10) vce(cluster turn#rep) verbose(3)
	reghdfe price weight (disp=length), a(rep) tol(1e-10) vce(cluster turn foreign) verbose(3)
	reghdfe price weight (disp=length), a(rep) tol(1e-10) vce(robust, bw(2)) // not passing it
	asd
	reghdfe price weight disp, a(rep) tol(1e-10) vce(cluster turn foreign) verbose(3)
	reghdfe price weight disp, a(rep) tol(1e-10) vce(cluster turn rep) verbose(3)

	
exit

* TODO
* Fix EstimateDoF
* Fix Estimate -> the call to the wrappers
* Fix ivreg2 wrapper
* Fix regress wrapper, with calls depending on the situation
* Create the custom mwc

