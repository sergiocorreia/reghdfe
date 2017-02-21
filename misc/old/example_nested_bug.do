clear all
set more off
cls

* Params
	local N 100
	local reps 1000

	local dtasize = max(`N', `reps')

* Create dataset
	set obs `N'
	gen id = _n in 1/`N'
	replace id = (id - 1) if mod(id,2)==0 & id<=10 in 1/`N'
	bys id: gen t = _n if id<.
	xtset id t
	by id: gen N = _N if id<.
	tab N
	gen byte singleton = N==1
	compress

	* X = xtreg_fe A = areg B = xtreg_fe
	* C=cluster O=ols R=robust B=bootstrap J=Jacknife
	* F=Full Sample S = Small Sample
	gen pvalue_AOF = .
	gen pvalue_XOF = .
	gen pvalue_AOS = .
	gen pvalue_XOS = .
	gen pvalue_ARF = .
	gen pvalue_XRF = . // --> This would be the same as clustered
	gen pvalue_ARS = .
	gen pvalue_XRS = . // --> This would be the same as clustered
	gen pvalue_ACF = .
	gen pvalue_XCF = .
	gen pvalue_ACS = .
	gen pvalue_XCS = .
	gen pvalue_XBF = .
	gen pvalue_XBS = .
	gen pvalue_XJF = .

	gen y = .
	gen x = .
	gen alpha = .
	gen e = .

	set seed 123456

forval rep = 1/`reps' {
	di as text "[rep=`rep']"

	replace x = rnormal() in 1/`N'
	replace alpha = rnormal() in 1/`N'
	by id: replace alpha = alpha[1] if id<.
	replace e = rnormal() // can add AR/MA terms within cluster
	replace y = alpha + 0 * x + e in 1/`N'

// -------------------------------------------------------------------------------------------------

	* AREG OLS all obs
	qui areg y x in 1/`N', absorb(id)
	test x // avoids writing the formula for the pvalue
	replace pvalue_AOF = r(p) in `rep'

	* XTREG OLS all obs
	qui xtreg y x in 1/`N', fe vce(conventional)
	test x // avoids writing the formula for the pvalue
	replace pvalue_XOF = r(p) in `rep'

	* AREG OLS useful obs
	qui areg y x if !singleton in 1/`N', absorb(id)
	test x // avoids writing the formula for the pvalue
	replace pvalue_AOS = r(p) in `rep'

	* XTREG OLS useful obs
	qui xtreg y x if !singleton in 1/`N', fe vce(conventional)
	test x // avoids writing the formula for the pvalue
	replace pvalue_XOS = r(p) in `rep'

// -------------------------------------------------------------------------------------------------

	* AREG ROBUST all obs
	qui areg y x in 1/`N', absorb(id) robust
	test x // avoids writing the formula for the pvalue
	replace pvalue_ARF = r(p) in `rep'

	* XTREG ROBUST all obs
	* N/A
	
	* AREG ROBUST useful obs
	qui areg y x if !singleton in 1/`N', absorb(id) robust
	test x // avoids writing the formula for the pvalue
	replace pvalue_ARS = r(p) in `rep'

	* XTREG ROBUST useful obs
	* N/A

// -------------------------------------------------------------------------------------------------

	* AREG CLUSTER all obs
	qui areg y x in 1/`N', absorb(id) vce(cluster id)
	test x // avoids writing the formula for the pvalue
	replace pvalue_ACF = r(p) in `rep'

	* AREG CLUSTER useful obs
	qui areg y x if !singleton in 1/`N', absorb(id) vce(cluster id)
	test x // avoids writing the formula for the pvalue
	replace pvalue_ACS = r(p) in `rep'

	* XTREG CLUSTER all obs
	qui xtreg y x in 1/`N', fe vce(cluster id)
	test x // avoids writing the formula for the pvalue
	replace pvalue_XCF = r(p) in `rep'

	* XTREG CLUSTER useful obs
	qui xtreg y x if !singleton in 1/`N', fe vce(cluster id)
	test x // avoids writing the formula for the pvalue
	replace pvalue_XCS = r(p) in `rep'

// -------------------------------------------------------------------------------------------------

	*** XTREG Bootstrap all obs
	**qui xtreg y x in 1/`N', fe vce(boot, reps(400))
	**test x // avoids writing the formula for the pvalue
	**replace pvalue_XBF = r(p) in `rep'

	*** XTREG Bootstrap useful obs
	**qui xtreg y x if !singleton in 1/`N', fe vce(boot, reps(400))
	**test x // avoids writing the formula for the pvalue
	**replace pvalue_XBS = r(p) in `rep'

	* XTREG Jacknife all obs
	qui xtreg y x in 1/`N', fe vce(jack)
	test x // avoids writing the formula for the pvalue
	replace pvalue_XBF = r(p) in `rep'

// -------------------------------------------------------------------------------------------------
	if mod(`rep',10)==0 {
		save "example_nested_bug.dta", replace
	}
}

// -------------------------------------------------------------------------------------------------

	local suffixes AOF XOF AOS XOS ARF XRF ARS XRS ACF XCF ACS XCS XBF XBS
	local levels 05 10
	foreach suffix of local suffixes {
		foreach level of local levels {
			gen accept`level'_`suffix' = (pvalue_`suffix' <= `level'/100) if pvalue_`suffix'<.
		}
	}
	format %5.3f accept*
	save "example_nested_bug.dta", replace

	su accept05*, sep(2) format
	su accept10*, sep(2) format
	tw (kdensity pvalue_AOF) (kdensity pvalue_XCF) (kdensity pvalue_XCS) (kdensity pvalue_XBF)

	collapse (mean) accept*, fast

exit


/*
cls
clear all
set more off

sysuse auto, clear
gen id = _n - (id<=6 & mod(id,2)==0)
tab id
bys id: gen t = _n
xtset id t
bys id: gen is_singleton = (_N==1)

xtreg price weight length, fe vce(cluster id)
drop if is_singleton
xtreg price weight length, fe vce(cluster id)



global vce vce(cluster id) // vce(robust) vce(cluster t)
cap noi xtreg price weight length, fe $vce // dfadj
cap noi areg price weight length, absorb(id) $vce
cap noi reghdfe price weight length, absorb(id) $vce

exit




cap noi xtreg price weight length, fe vce(boot, reps(500) seed(10101))
*/



