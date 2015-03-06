cd "D:/Github/reghdfe/source"
cscript "reghdfe comparison with ivreg2" adofile reghdfe

* Setup
	discard
	clear all
	set more off
	* cls

* Convenience: "Trim <size>" will trim e(b) and e(V)
	capture program drop TrimMatrix
	program define TrimMatrix, eclass
	args size scale
		assert `size'>0
		matrix trim_b = e(b)
		matrix trim_V = e(V)
		matrix trim_b = trim_b[1, 1..`size']
		matrix trim_V = trim_V[1..`size',1..`size']
		if ("`scale'"!="") matrix trim_V = trim_V * `scale'
		ereturn matrix trim_b = trim_b
		ereturn matrix trim_V = trim_V
	end
	
* Create fake dataset
	sysuse auto
	*gen n = int(uniform()*10+3) // used for weights
	*replace length = 0 if rep==3 // used for DoF adjustment of cont var
	*replace length = 5 if rep==1
	*gen byte one = 1
	bys turn: gen t = _n
	tsset turn t
	drop if missing(rep)

	local include_ols  ///
		scalar: N rmse rss r2 df_r F ///
		matrix: trim_b trim_V ///
		macros: wexp wtype 

* [TEST]
	
	noi di as text " -  Comparison with ivreg2"
	local depvar price
	local endogvars weight length
	local indepvars trunk
	local instruments gear displace head
	local absvars turn
	*local clustervars rep
	fvunab tmp : `indepvars' `endogvars'
	local K : list sizeof tmp

	qui tab `absvars', gen(ABS_)

	* (ols) Clustered-xtivreg2; clustervar==absvar
	xtivreg2 `depvar' `indepvars' `endogvars', cluster(`absvars') small i(`absvars') fe t(t)
	TrimMatrix `K'
	matrix list e(trim_V)
	storedresults save benchmark e()

	reghdfe `depvar' `indepvars' `endogvars', absorb(`absvars') vce(cluster `absvars', suite(avar)) nocon dropsingletons
	estadd scalar r2 = e(r2_within), replace
	local scale = e(unclustered_df_r) / (e(unclustered_df_r) + 1)
	TrimMatrix `K' `scale'
	estadd scalar F = e(F) / `scale', replace
	matrix list e(trim_V)
	storedresults compare benchmark e(), tol(1e-8) include(`include_ols')
	storedresults drop benchmark

// -------------------------------------------------------------------------------------------------

	* (iv) Unadjusted-ivreg2
	local clustervars // Empty
	ivreg2 `depvar' `indepvars' ABS_* (`endogvars'=`instruments'), partial(ABS_*) small
	TrimMatrix `K'
	storedresults save benchmark e()

	reghdfe `depvar' `indepvars' (`endogvars'=`instruments'), absorb(`absvars') vce(unadjusted, suite(avar)) tol(1e-12) nocon // dropsingletons
	estadd scalar r2 = e(r2_within), replace
	*local scale = 1 // e(unclustered_df_r) / (e(unclustered_df_r) + 1)
	TrimMatrix `K' // `scale'
	*estadd scalar F = e(F) / `scale', replace
	*matrix list e(trim_V)
	storedresults compare benchmark e(), tol(1e-6) include(`include_ols')
	storedresults drop benchmark

// -------------------------------------------------------------------------------------------------

	* (iv) Robust-ivreg2
	ivreg2 `depvar' `indepvars' ABS_* (`endogvars'=`instruments'), partial(ABS_*) small robust
	TrimMatrix `K'
	storedresults save benchmark e()

	reghdfe `depvar' `indepvars' (`endogvars'=`instruments'), absorb(`absvars') vce(robust, suite(avar)) tol(1e-12) nocon // dropsingletons
	estadd scalar r2 = e(r2_within), replace
	*local scale = 1 // e(unclustered_df_r) / (e(unclustered_df_r) + 1)
	TrimMatrix `K' // `scale'
	*estadd scalar F = e(F) / `scale', replace
	*matrix list e(trim_V)
	storedresults compare benchmark e(), tol(1e-6) include(`include_ols')
	storedresults drop benchmark

// -------------------------------------------------------------------------------------------------

	* (iv) Clustered-ivreg2
	local clustervars rep
	ivreg2 `depvar' `indepvars' ABS_* (`endogvars'=`instruments'), partial(ABS_*) cluster(`clustervars') small
	TrimMatrix `K'
	storedresults save benchmark e()

	reghdfe `depvar' `indepvars' (`endogvars'=`instruments'), absorb(`absvars') vce(cluster `clustervars', suite(avar)) tol(1e-12) nocon // dropsingletons
	estadd scalar r2 = e(r2_within), replace
	*local scale = 1 // e(unclustered_df_r) / (e(unclustered_df_r) + 1)
	TrimMatrix `K' // `scale'
	*estadd scalar F = e(F) / `scale', replace
	*matrix list e(trim_V)
	storedresults compare benchmark e(), tol(1e-6) include(`include_ols')
	storedresults drop benchmark

// -------------------------------------------------------------------------------------------------

cd "D:/Github/reghdfe/test"
exit
