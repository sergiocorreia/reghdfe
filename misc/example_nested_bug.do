clear all
set more off
cls

* Params
	local N 100
	local reps 100

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

	* X = xtreg_fe A = areg B = xtreg_fe
	* C=cluster O=ols B=bootstrap
	* F=Full Sample S = Small Sample
	gen pvalue_AO = .
	gen pvalue_XO = .
	gen pvalue_ACF = .
	gen pvalue_XCF = .
	gen pvalue_ACS = .
	gen pvalue_XCS = .
	gen pvalue_XBF = .
	gen y = .
	gen x = .
	gen alpha = .
	gen e = .

	set seed 12345

forval rep = 1/`reps' {

	* Update data	
	replace x = rnormal() in 1/`N'
	replace alpha = rnormal() in 1/`N'
	by id: replace alpha = alpha[1] if id<.
	replace e = rnormal() // can add AR/MA terms within cluster
	replace y = alpha + 0 * x + e in 1/`N'

	local vce  // Cluster by ID

	* AREG OLS all obs
	areg y x in 1/`N', absorb(id)
	test x // avoids writing the formula for the pvalue
	replace pvalue_AO = r(p)

	* XTREG OLS all obs
	xtreg y x in 1/`N', fe vce(conventional)
	test x // avoids writing the formula for the pvalue
	replace pvalue_XO = r(p)

	* AREG CLUSTER all obs
	areg y x in 1/`N', absorb(id) vce(cluster id)
	test x // avoids writing the formula for the pvalue
	replace pvalue_ACF = r(p)

	* AREG CLUSTER useful obs
	areg y x if !singleton in 1/`N', absorb(id) vce(cluster id)
	test x // avoids writing the formula for the pvalue
	replace pvalue_ACS = r(p)

	* XTREG CLUSTER all obs
	xtreg y x in 1/`N', fe vce(cluster id)
	test x // avoids writing the formula for the pvalue
	replace pvalue_XCF = r(p)

	* XTREG CLUSTER useful obs
	xtreg y x if !singleton in 1/`N', fe vce(cluster id)
	test x // avoids writing the formula for the pvalue
	replace pvalue_XCS = r(p)

	* XTREG Bootstrap all obs
	xtreg y x in 1/`N', fe vce(boot, reps(400))
	test x // avoids writing the formula for the pvalue
	replace pvalue_XBF = r(p)

	su pvalue*, sep(2)
}

collapse (mean) pvalue*, fast
save "example_nested_bug.dta", replace


exit







cls
clear all
set more off
sysuse auto

gen id = _n
replace id = 1 if id==2
replace id = 3 if id==4
replace id = 5 if id==6
bys id: gen t = _n


xtset id t
xtreg price weight length if id<=6, fe vce(cluster id)

global vce vce(cluster id) // vce(robust) vce(cluster t)
cap noi xtreg price weight length, fe $vce // dfadj
cap noi areg price weight length, absorb(id) $vce
cap noi reghdfe price weight length, absorb(id) $vce

exit




cap noi xtreg price weight length, fe vce(boot, reps(500) seed(10101))
