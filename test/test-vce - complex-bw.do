cd "D:/Github/reghdfe/source"
cscript "reghdfe with weights" adofile reghdfe

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
	
	noi di as text " -  VCE with absvars==clustervars==tsset, and bw(2)"
	local lhs price
	local rhs weight length
	local absvars turn t
	fvunab tmp : `rhs'
	local K : list sizeof tmp

	* 1. Run benchmark
	qui tab turn, gen(TURN_)
	qui tab t, gen(T_)
	ivreg2 `lhs' `rhs' TURN_* T_*, partial(TURN_* T_*) cluster(`absvars') bw(2)
	TrimMatrix `K'
	storedresults save benchmark e()
	
	* 2. Run reghdfe
	reghdfe `lhs' `rhs', a(`absvars') vce(cluster `absvars', bw(2)) verbose(3)
	TrimMatrix `K'
	
	* 3. Compare
	storedresults compare benchmark e(), tol(1e-12) include( ///
		scalar: N rss r2 r2_a F df_r df_a df_m ///
		matrix: trim_b trim_V ///
		macros: wexp wtype )
	storedresults drop benchmark

cd "D:/Github/reghdfe/test"
exit
