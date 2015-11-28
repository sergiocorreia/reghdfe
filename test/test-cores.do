cscript "reghdfe with multicore" adofile reghdfe

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
	
	local included_e scalar: N rmse tss rss r2 r2_a F df_r df_a df_m F_absorb ///
		matrix: trim_b trim_V ///
		macros: wexp wtype

	local lhs price
	local rhs weight length gear disp
	local absvars turn trunk
	fvunab tmp : `rhs'
	local K : list sizeof tmp
	
	* 1. Run benchmark w/1 core
	reghdfe `lhs' `rhs', absorb(`absvars')
	TrimMatrix `K'
	storedresults save benchmark e()
	
	* 2. Run with 2 cores and compare
	reghdfe `lhs' `rhs', absorb(`absvars') cores(2)
	TrimMatrix `K'
	storedresults compare benchmark e(), tol(1e-12) exclude(macros: cmdline)

	* 3. Run with 4 cores and compare
	reghdfe `lhs' `rhs', absorb(`absvars') cores(4)
	TrimMatrix `K'
	storedresults compare benchmark e(), tol(1e-12) exclude(macros: cmdline)

	* 2. Run with 3 cores and -fast-
	reghdfe `lhs' `rhs', absorb(`absvars') cores(3) fast
	TrimMatrix `K'
	storedresults compare benchmark e(), tol(1e-12) exclude(macros: cmdline corr1 corr2)

	* 2. Run with 1 core and assert that -fast- works
	reghdfe `lhs' `rhs', absorb(`absvars') cores(1) fast
	TrimMatrix `K'
	storedresults compare benchmark e(), tol(1e-12) exclude(macros: cmdline corr1 corr2)


	storedresults drop benchmark

cd "C:/Git/reghdfe/test"
exit
