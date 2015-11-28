cscript "reghdfe template" adofile reghdfe

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

* [TEST] Simple example
	
	noi di as text " -  description of test"
	local lhs price
	local rhs weight length
	local absvars turn
	fvunab tmp : `rhs'
	local K : list sizeof tmp

	* Custom adjustment: in this simple case we can compare _cons
	local K = `K' + 1

	* 1. Run benchmark
	areg `lhs' `rhs', absorb(`absvars')
	TrimMatrix `K'
	storedresults save benchmark e()
	
	* 2. Run reghdfe
	reghdfe `lhs' `rhs', absorb(`absvars')
	TrimMatrix `K'
	
	* 3. Compare
	storedresults compare benchmark e(), tol(1e-12) include( ///
		scalar: N rmse tss rss mss r2 r2_a F df_r df_a df_m F_absorb ///
		matrix: trim_b trim_V ///
		macros: wexp wtype )
	storedresults drop benchmark

cd "C:/Git/reghdfe/test"
exit
