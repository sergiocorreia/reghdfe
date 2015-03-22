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

	* XTREG REQUIRES WEIGHTS TO BE CONSTANT WITHIN ABSVAR
	bys idcode (year): gen double cfw = `fweightvar'[1]
	local cfweightvar cfw

	
	fvunab tmp : `indepvars' `endogvars'
	local K : list sizeof tmp

	local include ///
		scalar: N N_clust df_r /// r2
		matrix: trim_b trim_V ///
		macros: wtype // wexp

	local i 0
rebuild_git reghdfe

// -------------------------------------------------------------------------------------------------

* [TEST] Freq. weight - Comparison vs XTREG
	local weight fweight
	noi di as text " - freq. weight vs xtreg"
	xtreg `depvar' `indepvars' [`weight'=`cfweightvar'], vce(cluster `clustervars') fe // absorb(`absvars') 
	TrimMatrix `K'
	matrix list e(trim_b), format(%20.12e)
	matrix list e(trim_V), format(%20.12e)
	di e(df_r)
	storedresults save benchmark e()

foreach suite in default avar mwc {
	di as text "SUITE=<`suite'> WEIGHTEXP=[`weight'=`cfweightvar']"
	reghdfe `depvar' `indepvars' [`weight'=`cfweightvar'], absorb(`absvars') vce(cluster `clustervars', suite(`suite'))
	TrimMatrix `K'
	matrix list e(trim_b), format(%20.12e)
	matrix list e(trim_V), format(%20.12e)
	storedresults compare benchmark e(), tol(1e-10) include(`include')
}

	storedresults drop benchmark

// -------------------------------------------------------------------------------------------------

* [TEST] Freq. weight - Comparison vs AREG
	local weight fweight
	noi di as text " - freq. weight vs areg"
	areg `depvar' `indepvars' [`weight'=`fweightvar'], vce(cluster `clustervars') absorb(`absvars') 
	TrimMatrix `K'
	matrix list e(trim_b), format(%20.12e)
	matrix list e(trim_V), format(%20.12e)
	di e(df_r)
	storedresults save benchmark e()

foreach suite in default avar mwc {
	di as text "SUITE=<`suite'> WEIGHTEXP=[`weight'=`fweightvar']"
	reghdfe `depvar' `indepvars' [`weight'=`fweightvar'], absorb(`absvars') vce(cluster `clustervars', suite(`suite')) dof(pairwise continuous) // exclude -clusters- from dof()
	TrimMatrix `K'
	matrix list e(trim_b), format(%20.12e)
	matrix list e(trim_V), format(%20.12e)
	storedresults compare benchmark e(), tol(1e-10) include(`include')
}

	storedresults drop benchmark

// -------------------------------------------------------------------------------------------------

* [TEST] Analytic. weight vs xtreg
	local weight aweight
	noi di as text " - analytic weight: [`weight'=`cfweightvar']"
	xtreg `depvar' `indepvars' [`weight'=`cfweightvar'], vce(cluster `clustervars') fe // absorb(`absvars') 
	di e(df_r) // Missing?
	TrimMatrix `K'
	matrix list e(trim_V), format(%20.12e)
	storedresults save benchmark e()

foreach suite in default avar mwc {
	di as text "SUITE=<`suite'> WEIGHTEXP=[`weight'=`cfweightvar']"
	reghdfe `depvar' `indepvars' [`weight'=`cfweightvar'], absorb(`absvars') vce(cluster `clustervars', suite(`suite'))
	TrimMatrix `K'
	matrix list e(trim_b), format(%20.12e)
	matrix list e(trim_V), format(%20.12e)
	storedresults compare benchmark e(), tol(1e-10) include(`include')
}

	storedresults drop benchmark

// -------------------------------------------------------------------------------------------------

* [TEST] Analytic. weight vs areg
	local weight aweight
	noi di as text " - analytic weight: [`weight'=`weightvar']"
	areg `depvar' `indepvars' [`weight'=`weightvar'], vce(cluster `clustervars') absorb(`absvars') 
	di e(df_r) // Missing?
	TrimMatrix `K'
	matrix list e(trim_V), format(%20.12e)
	storedresults save benchmark e()

foreach suite in default avar mwc {
	di as text "SUITE=<`suite'> WEIGHTEXP=[`weight'=`weightvar']"
	reghdfe `depvar' `indepvars' [`weight'=`weightvar'], absorb(`absvars') vce(cluster `clustervars', suite(`suite')) dof(pairwise continuous) // exclude -clusters- from dof()
	TrimMatrix `K'
	matrix list e(trim_b), format(%20.12e)
	matrix list e(trim_V), format(%20.12e)
	storedresults compare benchmark e(), tol(1e-10) include(`include')
}

	storedresults drop benchmark

// -------------------------------------------------------------------------------------------------

* [TEST] Prob. weight vs xtreg
	local weight pweight
	noi di as text " - prob. weight: [`weight'=`cfweightvar']"
	xtreg `depvar' `indepvars' [`weight'=`cfweightvar'], vce(cluster `clustervars') fe // absorb(`absvars') 
	di e(df_r) // Missing?
	TrimMatrix `K'
	matrix list e(trim_V), format(%20.12e)
	storedresults save benchmark e()

foreach suite in default avar mwc {
	di as text "SUITE=<`suite'> WEIGHTEXP=[`weight'=`cfweightvar']"
	reghdfe `depvar' `indepvars' [`weight'=`cfweightvar'], absorb(`absvars') vce(cluster `clustervars', suite(`suite'))
	TrimMatrix `K'
	matrix list e(trim_b), format(%20.12e)
	matrix list e(trim_V), format(%20.12e)
	storedresults compare benchmark e(), tol(1e-10) include(`include')
}

	storedresults drop benchmark

// -------------------------------------------------------------------------------------------------

* [TEST] Prob. weight vs areg
	local weight pweight
	noi di as text " - prob. weight: [`weight'=`weightvar']"
	areg `depvar' `indepvars' [`weight'=`weightvar'], vce(cluster `clustervars') absorb(`absvars') 
	di e(df_r) // Missing?
	TrimMatrix `K'
	matrix list e(trim_V), format(%20.12e)
	storedresults save benchmark e()

foreach suite in default avar mwc {
	di as text "SUITE=<`suite'> WEIGHTEXP=[`weight'=`weightvar']"
	reghdfe `depvar' `indepvars' [`weight'=`weightvar'], absorb(`absvars') vce(cluster `clustervars', suite(`suite')) dof(pairwise continuous) // exclude -clusters- from dof()
	TrimMatrix `K'
	matrix list e(trim_b), format(%20.12e)
	matrix list e(trim_V), format(%20.12e)
	storedresults compare benchmark e(), tol(1e-10) include(`include')
}

	storedresults drop benchmark

// -------------------------------------------------------------------------------------------------

cd "D:/Github/reghdfe/test"
exit
