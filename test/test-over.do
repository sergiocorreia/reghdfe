noi cscript "reghdfe with cache" adofile reghdfe

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
	replace foreign = foreign + 10
	label define origin 10 "Local" 11 "It's foreign", add
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
	local othervar gear head
	local absvars turn
	local over foreign
	fvunab tmp : `rhs'
	local K : list sizeof tmp

	local include ///
		scalar: N rmse tss rss mss r2 r2_a F df_r df_m F_absorb ///
		matrix: trim_b trim_V ///
		macros: wexp wtype 

	* 1. Run benchmarks
	areg `lhs' `rhs' if foreign==10, absorb(`absvars')
	TrimMatrix `K'
	local bench_df_a0 = e(df_a)
	storedresults save benchmark0 e()

	areg `lhs' `rhs' if foreign==11, absorb(`absvars')
	TrimMatrix `K'
	local bench_df_a1 = e(df_a)
	storedresults save benchmark1 e()
	
	* 2. Save cache
	local fn "D:/Github/tmp/thecache"
	
	reghdfe `lhs' `othervar' `rhs', absorb(`absvars') savecache("`fn'") over(foreign) cores(2)
	assert "`e(over_levels)'"!=""

	* 3. Use cache and compare
	reghdfe `lhs' `rhs' if foreign==10, absorb(`absvars') usecache("`fn'") over(foreign) cores(2)
	TrimMatrix `K'
	assert `bench_df_a0'==e(df_a)-1
	storedresults compare benchmark0 e(), tol(1e-12) include(`include')

	reghdfe `lhs' `rhs' if foreign==11, absorb(`absvars') usecache("`fn'") over(foreign) cores(2)
	TrimMatrix `K'
	assert `bench_df_a1'==e(df_a)-1
	storedresults compare benchmark1 e(), tol(1e-12) include(`include')

	storedresults drop benchmark0 benchmark1

	// -------------------------------------------------------------------------------------------------

	* 1. Run benchmarks
	areg `lhs' `rhs' if foreign==10, absorb(`absvars')
	local bench_df_a0 = e(df_a)
	TrimMatrix `K'
	storedresults save benchmark0 e()

	areg `lhs' `rhs' if foreign==11, absorb(`absvars')
	local bench_df_a1 = e(df_a)
	TrimMatrix `K'
	storedresults save benchmark1 e()
	
	* 2. Save cache
	local fn "D:/Github/tmp/thecache"
	set trace off
	reghdfe `lhs' `othervar' `rhs', absorb(`absvars') savecache("`fn'") over(foreign)

	* 3. Use cache and compare
	reghdfe `lhs' `rhs' if foreign==10, absorb(`absvars') usecache("`fn'") over(foreign) summarize
	TrimMatrix `K'
	assert `bench_df_a0'==e(df_a)-1
	storedresults compare benchmark0 e(), tol(1e-12) include(`include')

	reghdfe `lhs' `rhs' if foreign==11, absorb(`absvars') usecache("`fn'") over(foreign) summarize
	TrimMatrix `K'
	assert `bench_df_a1'==e(df_a)-1
	storedresults compare benchmark1 e(), tol(1e-12) include(`include')

	storedresults drop benchmark0 benchmark1

cd "D:/Github/reghdfe/test"
exit
