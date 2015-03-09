cd "D:/Github/reghdfe" // /source
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
	local lhs price
	local rhs weight length
	local othervar gear head
	local absvars turn
	fvunab tmp : `rhs'
	local K : list sizeof tmp

	local include ///
		scalar: N rmse tss rss mss r2 r2_a F df_r df_a df_m F_absorb ///
		matrix: trim_b trim_V ///
		macros: wexp wtype 

	* Custom adjustment: in this simple case we can compare _cons
	local K = `K' + 1

	* 1. Run benchmark
	areg `lhs' `rhs', absorb(`absvars')
	TrimMatrix `K'
	storedresults save benchmark e()
	
	* 2. Save cache
	local fn "D:/Github/tmp/thecache"
	reghdfe `lhs' `othervar' `rhs', absorb(`absvars') savecache("`fn'")

	* 3. Use cache
	reghdfe `lhs' `othervar', absorb(`absvars') usecache("`fn'")
	reghdfe `lhs' `rhs', absorb(`absvars') usecache("`fn'")
	TrimMatrix `K'
	
	* 3. Compare
	storedresults compare benchmark e(), tol(1e-12) include(`include')

	* Repeat with -fast- in usecache
	local fn "D:/Github/tmp/thecache"
	reghdfe `lhs' `othervar' `rhs', absorb(`absvars') savecache("`fn'")
	reghdfe `lhs' `rhs', absorb(`absvars') usecache("`fn'") fast
	TrimMatrix `K'
	storedresults compare benchmark e(), tol(1e-12) include(`include')

	* Repeat with multiple cores
	local fn "D:/Github/tmp/thecache"
	reghdfe `lhs' `othervar' `rhs', absorb(`absvars') savecache("`fn'") cores(3)
	reghdfe `lhs' `rhs', absorb(`absvars') usecache("`fn'") fast
	TrimMatrix `K'
	storedresults compare benchmark e(), tol(1e-12) include(`include')

	storedresults drop benchmark

cd "D:/Github/reghdfe/test"
exit
