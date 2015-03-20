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
	local absvars idcode // age
	local clustervars idcode // age

	drop if `weightvar'==. | `weightvar'<=0
	drop if `fweightvar'==. | `fweightvar'<=0
	
	fvunab tmp : `indepvars' `endogvars'
	local K : list sizeof tmp

	local include ///
		scalar: N N_clust df_r /// r2
		matrix: trim_b trim_V ///
		macros: wtype // wexp

	local i 0

// -------------------------------------------------------------------------------------------------
rebuild_git reghdfe
* [TEST] Freq. weight
	local weight fweight
	noi di as text " - freq. weight"
	xtreg `depvar' `indepvars' [`weight'=`fweightvar'], vce(cluster `clustervars') fe // absorb(`absvars') 
	TrimMatrix `K'
	matrix list e(trim_b), format(%20.12e)
	matrix list e(trim_V), format(%20.12e)
	di e(df_r)
	storedresults save benchmark e()

foreach suite in default avar mwc {
	di as text "SUITE=<`suite'> WEIGHTEXP=[`weight'=`fweightvar']"
	reghdfe `depvar' `indepvars' [`weight'=`fweightvar'], absorb(`absvars') vce(cluster `clustervars', suite(`suite')) tol(1e-12)
	TrimMatrix `K'
	matrix list e(trim_b), format(%20.12e)
	matrix list e(trim_V), format(%20.12e)
	storedresults compare benchmark e(), tol(1e-10) include(`include')
}

	storedresults drop benchmark
asd
// -------------------------------------------------------------------------------------------------

* [TEST] Analytic. weight
	local weight aweight
	noi di as text " - analytic weight: [`weight'=`weightvar']"
	xtreg `depvar' `indepvars' [`weight'=`weightvar'], vce(cluster `clustervars') fe // absorb(`absvars') 
	di e(df_r) // Missing?
	TrimMatrix `K'
	matrix list e(trim_V), format(%20.12e)
	storedresults save benchmark e()

	reghdfe `depvar' `indepvars' [`weight'=`weightvar'], absorb(`absvars') vce(cluster `clustervars', suite(avar))
	TrimMatrix `K'
	matrix list e(trim_V), format(%20.12e)
	storedresults compare benchmark e(), tol(1e-10) include(`include')

	reghdfe `depvar' `indepvars' [`weight'=`weightvar'], absorb(`absvars') vce(cluster `clustervars', suite(mwc))
	TrimMatrix `K'
	matrix list e(trim_V), format(%20.12e)
	storedresults compare benchmark e(), tol(1e-10) include(`include')

	storedresults drop benchmark

// -------------------------------------------------------------------------------------------------

* [TEST] Prob. weight
	local weight pweight
	noi di as text " - prob. weight: [`weight'=`weightvar']"
	xtreg `depvar' `indepvars' [`weight'=`weightvar'], vce(cluster `clustervars') fe // absorb(`absvars') 
	di e(df_r) // Missing?
	TrimMatrix `K'
	matrix list e(trim_V), format(%20.12e)
	storedresults save benchmark e()

	reghdfe `depvar' `indepvars' [`weight'=`weightvar'], absorb(`absvars') vce(cluster `clustervars', suite(avar))
	TrimMatrix `K'
	matrix list e(trim_V), format(%20.12e)
	storedresults compare benchmark e(), tol(1e-10) include(`include')

	reghdfe `depvar' `indepvars' [`weight'=`weightvar'], absorb(`absvars') vce(cluster `clustervars', suite(mwc))
	TrimMatrix `K'
	matrix list e(trim_V), format(%20.12e)
	storedresults compare benchmark e(), tol(1e-10) include(`include')

	storedresults drop benchmark

// -------------------------------------------------------------------------------------------------

cd "D:/Github/reghdfe/test"
exit
