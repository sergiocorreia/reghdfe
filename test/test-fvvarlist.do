noi cscript "reghdfe with factor variables" adofile reghdfe

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
	*gen byte spam = price>6000
	bys turn: gen t = _n
	tsset turn t
	drop if missing(rep)
	
	local included_e scalar: N rmse tss rss r2 r2_a F df_r df_m ll ll_0 /// F_absorb
		matrix: /// trim_b trim_V ///
		macros: wexp wtype

* [TEST] Undjusted
	local lhs price
	local rhs weight foreign##trunk##c.length
	local absvars turn
	fvunab tmp : `rhs'
	local K : list sizeof tmp
	
	* 1. Run benchmark
	areg `lhs' `rhs', absorb(`absvars')
	TrimMatrix `K'
	local bench_df_a = e(df_a)
	storedresults save benchmark e()
	
	* 2. Run reghdfe and compare
	reghdfe `lhs' `rhs', absorb(`absvars') keepsingletons tol(1e-10)
	TrimMatrix `K'
	storedresults compare benchmark e(), tol(1e-10) include(`included_e')
	assert `bench_df_a'==e(df_a)-1

	storedresults drop benchmark

cd "C:/Git/reghdfe/test"
exit
