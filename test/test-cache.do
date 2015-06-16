cscript "reghdfe with cache" adofile reghdfe

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
	local lhs L.price
	local rhs weight length
	local othervar gear head
	local absvars turn
	fvunab tmp : `rhs'
	local K : list sizeof tmp

	local exclude scalar: F_absorb macros: cmdline


	*	scalar: N rmse tss rss mss r2 r2_a F df_r df_m F_absorb ///
	*	matrix: trim_b trim_V ///
	*	macros: wexp wtype 

	* 1. Run benchmark
	reghdfe `lhs' `rhs', absorb(`absvars')
	storedresults save benchmark e()
	
	* 2. Save cache
	preserve
	reghdfe `lhs' `othervar' `rhs', absorb(`absvars') cache(save)

	* 3. Use cache
	reghdfe `lhs' `othervar', absorb(`absvars') cache(use)
	reghdfe `lhs' `rhs', absorb(`absvars') cache(use)
	
	* 3. Compare
	storedresults compare benchmark e(), tol(1e-12) exclude(`exclude')

	* 4. Cleanup
	reghdfe, cache(clear)
	restore

	* Repeat with -fast- in cache(use), a keep(), and vce(cluster)
	reghdfe `lhs' `rhs', absorb(`absvars') vce(cluster `absvars')
	storedresults save benchmark e()

	reghdfe `lhs' `othervar' `rhs', absorb(`absvars') cache(save, keep(foreign)) fast vce(cluster `absvars')
	reghdfe `lhs' `rhs', absorb(`absvars') cache(use) fast vce(cluster `absvars')
	storedresults compare benchmark e(), tol(1e-12) exclude(`exclude')

	reghdfe `lhs' `rhs', absorb(`absvars') cache(use) fast vce(cluster `absvars') resid(resid)
	bys foreign: su resid

	* Cleanup
	reghdfe, cache(clear)

	storedresults drop benchmark
/*
* [TEST] Complex lags example
	sysuse auto, clear
	bys turn: gen t = _n
	tsset turn t

	noi di as text " -  complex factor/time vars"
	local lhs L.price
	local rhs ibn.foreign#c.(L(-1/1).length L(-1/1).weight)
	local absvars turn
	fvunab tmp : `rhs'
	local K : list sizeof tmp

	local include ///
		scalar: N rmse tss rss mss r2 r2_a F df_r df_m F_absorb ///
		matrix: trim_b trim_V ///
		macros: wexp wtype 

	* 1. Run benchmark
	areg `lhs' `rhs', absorb(`absvars')
	TrimMatrix `K'
	local bench_df_a = e(df_a)
	storedresults save benchmark e()
	
	* 2. Save cache
	local fn "D:/Github/tmp/thecache"
	reghdfe `lhs' `othervar' `rhs', absorb(`absvars') savecache("`fn'")

	* 3. Use cache
	reghdfe `lhs' `rhs', absorb(`absvars') usecache("`fn'")
	TrimMatrix `K'
	
	* 3. Compare
	assert `bench_df_a'==e(df_a)-1
	storedresults compare benchmark e(), tol(1e-10) include(`include')

*/
cd "D:/Github/reghdfe/test"
exit

