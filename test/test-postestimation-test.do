cscript "reghdfe postestimation: test" adofile reghdfe

* Setup
	discard
	clear all
	set more off
	* cls

	local included_e scalar: N rmse tss rss r2 r2_a F df_r df_m ///
		matrix: trim_b trim_V ///
		macros: wexp wtype

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

	
	noi di as text " -  description of test"
	local lhs price
	local rhs weight length gear disp
	local absvars turn
	fvunab tmp : `rhs'
	local K : list sizeof tmp

	local testvars length gear

* [TEST] Unadjusted

	* 1. Run benchmark
	areg `lhs' `rhs', absorb(`absvars')
	TrimMatrix `K'
	local bench_df_a = e(df_a)
	test `testvars'
	storedresults save bench_r r()
	storedresults save bench_e e()

	* 2. Run reghdfe and compare
	
	* 2a) vce suite = avar
	reghdfe `lhs' `rhs', absorb(`absvars') vce(unadjusted, suite(avar)) keepsingletons
	TrimMatrix `K'
	assert `bench_df_a'==e(df_a)-1
	storedresults compare bench_e e(), tol(1e-12) include(`included_e')
	test `testvars'
	storedresults compare bench_r r(), tol(1e-12) // include(scalar: drop df_r F df p)

	* 2a) vce suite = default
	reghdfe `lhs' `rhs', absorb(`absvars') vce(unadjusted, suite(default)) keepsingletons
	TrimMatrix `K'
	assert `bench_df_a'==e(df_a)-1
	storedresults compare bench_e e(), tol(1e-12) include(`included_e')
	test `testvars'
	storedresults compare bench_r r(), tol(1e-12) // include(scalar: drop df_r F df p)

	storedresults drop bench_e bench_r


* [TEST] Robust

	* 1. Run benchmark
	areg `lhs' `rhs', absorb(`absvars') robust
	TrimMatrix `K'
	local bench_df_a = e(df_a)
	test `testvars'
	storedresults save bench_r r()
	storedresults save bench_e e()

	* 2. Run reghdfe and compare
	
	* 2a) vce suite = avar
	reghdfe `lhs' `rhs', absorb(`absvars') vce(robust, suite(avar)) keepsingletons
	TrimMatrix `K'
	assert `bench_df_a'==e(df_a)-1
	storedresults compare bench_e e(), tol(1e-10) include(`included_e') // Note the lowered tol
	test `testvars'
	storedresults compare bench_r r(), tol(1e-12) // include(scalar: drop df_r F df p)

	* 2a) vce suite = default
	reghdfe `lhs' `rhs', absorb(`absvars') vce(robust, suite(default)) keepsingletons
	TrimMatrix `K'
	assert `bench_df_a'==e(df_a)-1
	storedresults compare bench_e e(), tol(1e-10) include(`included_e') // Note the lowered tol
	test `testvars'
	storedresults compare bench_r r(), tol(1e-12) // include(scalar: drop df_r F df p)

	storedresults drop bench_e bench_r
	

* [TEST] Cluster

	* 1. Run benchmark
	areg `lhs' `rhs', absorb(`absvars') cluster(trunk)
	TrimMatrix `K'
	local bench_df_a = e(df_a)
	test `testvars'
	storedresults save bench_r r()
	storedresults save bench_e e()

	* 2. Run reghdfe and compare
	
	* 2a) vce suite = avar
	reghdfe `lhs' `rhs', absorb(`absvars') vce(cluster trunk, suite(avar)) keepsingletons
	TrimMatrix `K'
	assert `bench_df_a'==e(df_a)-1
	storedresults compare bench_e e(), tol(1e-10) include(`included_e') // Note the lowered tol
	test `testvars'
	storedresults compare bench_r r(), tol(1e-12) // include(scalar: drop df_r F df p)

	* 2a) vce suite = default
	reghdfe `lhs' `rhs', absorb(`absvars') vce(cluster trunk, suite(default)) keepsingletons
	TrimMatrix `K'
	assert `bench_df_a'==e(df_a)-1
	storedresults compare bench_e e(), tol(1e-10) include(`included_e') // Note the lowered tol
	test `testvars'
	storedresults compare bench_r r(), tol(1e-12) // include(scalar: drop df_r F df p)

	storedresults drop bench_e bench_r
cd "D:/Github/reghdfe/test"
exit
