cd "D:/Github/reghdfe" // /source
cscript "reghdfe with weights and MWC" adofile reghdfe

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
	use "D:/Github/reghdfe/test/data/nlswork"

* Variables
	local depvar ln_wage
	local indepvars wks_work ttl_exp
	local weightvar tenure
	local fweightvar ind_code
	local absvars age hours
	local clustervars age hours
	
	fvunab tmp : `indepvars' `endogvars'
	local K : list sizeof tmp

	local include ///
		scalar: N N_clust /// r2 df_r
		matrix: trim_b trim_V ///
		macros: wtype // wexp

	local i 0
	foreach absvar of local absvars {
		qui tab `absvar', gen(ABS`++i'_)
	}

	drop if `weightvar'==. | `weightvar'<=0

// -------------------------------------------------------------------------------------------------

* [TEST] Freq. weight
	local weight fweight
	noi di as text " - freq. weight"
	ivreg2 `depvar' `indepvars' ABS* [`weight'=`fweightvar'], cluster(`clustervars') partial(ABS*) small
	di e(df_r) // Missing?
	TrimMatrix `K'
	matrix list e(trim_V), format(%20.12e)
	storedresults save benchmark e()

	reghdfe `depvar' `indepvars' [`weight'=`fweightvar'], absorb(`absvars') vce(cluster `clustervars', suite(avar))
	TrimMatrix `K'
	matrix list e(trim_V), format(%20.12e)
	storedresults compare benchmark e(), tol(1e-8) include(`include')

	reghdfe `depvar' `indepvars' [`weight'=`fweightvar'], absorb(`absvars') vce(cluster `clustervars', suite(mwc))
	TrimMatrix `K'
	matrix list e(trim_V), format(%20.12e)
	storedresults compare benchmark e(), tol(1e-8) include(`include')

	storedresults drop benchmark

// -------------------------------------------------------------------------------------------------

* [TEST] Analytic. weight
	local weight aweight
	noi di as text " - analytic weight"
	ivreg2 `depvar' `indepvars' ABS* [`weight'=`weightvar'], cluster(`clustervars') partial(ABS*) small
	di e(df_r) // Missing?
	TrimMatrix `K'
	matrix list e(trim_V), format(%20.12e)
	storedresults save benchmark e()

	reghdfe `depvar' `indepvars' [`weight'=`weightvar'], absorb(`absvars') vce(cluster `clustervars', suite(avar))
	TrimMatrix `K'
	matrix list e(trim_V), format(%20.12e)
	storedresults compare benchmark e(), tol(1e-8) include(`include')

	reghdfe `depvar' `indepvars' [`weight'=`weightvar'], absorb(`absvars') vce(cluster `clustervars', suite(mwc))
	TrimMatrix `K'
	matrix list e(trim_V), format(%20.12e)
	storedresults compare benchmark e(), tol(1e-8) include(`include')

	storedresults drop benchmark

// -------------------------------------------------------------------------------------------------

* [TEST] Prob. weight
	local weight pweight
	noi di as text " - prob. weight"
	ivreg2 `depvar' `indepvars' ABS* [`weight'=`weightvar'], cluster(`clustervars') partial(ABS*) small
	di e(df_r) // Missing?
	TrimMatrix `K'
	matrix list e(trim_V), format(%20.12e)
	storedresults save benchmark e()

	reghdfe `depvar' `indepvars' [`weight'=`weightvar'], absorb(`absvars') vce(cluster `clustervars', suite(avar))
	TrimMatrix `K'
	matrix list e(trim_V), format(%20.12e)
	storedresults compare benchmark e(), tol(1e-8) include(`include')

	reghdfe `depvar' `indepvars' [`weight'=`weightvar'], absorb(`absvars') vce(cluster `clustervars', suite(mwc))
	TrimMatrix `K'
	matrix list e(trim_V), format(%20.12e)
	storedresults compare benchmark e(), tol(1e-8) include(`include')

	storedresults drop benchmark

// -------------------------------------------------------------------------------------------------

cd "D:/Github/reghdfe/test"
exit
