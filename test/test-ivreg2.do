noi cscript "reghdfe comparison with ivreg2" adofile reghdfe

* Setup
	discard
	clear all
	set more off
	* cls

* Convenience: "Trim <size>" will trim e(b) and e(V)
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
	*gen n = int(uniform()*10+3) // used for weights
	*replace length = 0 if rep==3 // used for DoF adjustment of cont var
	*replace length = 5 if rep==1
	*gen byte one = 1
	bys turn: gen t = _n
	tsset turn t
	drop if missing(rep)

* [TEST] One absvar
	
	noi di as text " -  Comparison with ivreg2"
	local lhs price
	local rhs weight length
	local absvars turn
	local clustervars rep
	fvunab tmp : `rhs'
	local K : list sizeof tmp

	qui tab `absvars', gen(TURN_)

	* 1. Run benchmark
	cap ivreg2 `lhs' `rhs' TURN_*, cluster(`clustervars') small
	assert _rc==0
	TrimMatrix `K'
	storedresults save benchmark e()
	
	* 2. Run reghdfe
	reghdfe `lhs' `rhs', absorb(`absvars') vce(cluster `clustervars')
	TrimMatrix `K'
	
	* 3. Compare
	* Ignored:
	* - tss: ivreg2 doesn't set it
	* - F df_a df_m -> Because ivreg2 doesn't absorb the FEs
	storedresults compare benchmark e(), tol(1e-8) include( ///
		scalar: N rmse rss mss r2 r2_a df_r ll ///
		matrix: trim_b trim_V ///
		macros: wexp wtype )
	storedresults drop benchmark

* [TEST] Two absvars
	local absvars `absvars' foreign
	
	* 1. Run benchmark
	cap ivreg2 `lhs' `rhs' TURN_* foreign, cluster(`clustervars') small
	assert _rc==0
	TrimMatrix `K'
	storedresults save benchmark e()
	
	* 2. Run reghdfe
	reghdfe `lhs' `rhs', absorb(`absvars') vce(cluster `clustervars')
	TrimMatrix `K'
	
	* 3. Compare
	* Ignored:
	* - tss: ivreg2 doesn't set it
	* - F df_a df_m -> Because ivreg2 doesn't absorb the FEs
	storedresults compare benchmark e(), tol(1e-8) include( ///
		scalar: N rmse rss mss r2 r2_a df_r ///
		matrix: trim_b trim_V ///
		macros: wexp wtype )
	storedresults drop benchmark

* [TEST] One absvar, TWO clustervars
	
	noi di as text " -  Comparison with ivreg2"
	local lhs price
	local rhs weight length
	local absvars foreign
	local clustervars rep turn
	fvunab tmp : `rhs'
	local K : list sizeof tmp

	qui tab `absvars', gen(FOREIGN_)

	* 1. Run benchmark
	cap ivreg2 `lhs' `rhs' FOREIGN_*, cluster(`clustervars') small
	assert _rc==0
	TrimMatrix `K'
	storedresults save benchmark e()
	
	* 2. Run reghdfe
	reghdfe `lhs' `rhs', absorb(`absvars') vce(cluster `clustervars', suite(avar))
	
	TrimMatrix `K'
	
	* 3. Compare
	* Ignored:
	* - tss: ivreg2 doesn't set it
	* - F df_a df_m -> Because ivreg2 doesn't absorb the FEs
	storedresults compare benchmark e(), tol(1e-8) include( ///
		scalar: N rmse rss mss r2 r2_a df_r ///
		matrix: trim_b trim_V ///
		macros: wexp wtype )
	storedresults drop benchmark

	
	
cd "D:/Github/reghdfe/test"
exit
