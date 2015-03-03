cd "D:/Github/reghdfe/source"
cscript "reghdfe clustering and absorbing by the same variable" adofile reghdfe

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

* [TEST]
	local lhs price
	local rhs weight length
	local absvars turn // make it the same as panelvar
	local clustervar `absvars'
	fvunab tmp : `rhs'
	local K : list sizeof tmp

	* Custom adjustment: in this simple case we can compare _cons
	drop if missing(rep)

// -------------------------------------------------------------------------------------------------
* Must match -areg- without the dof adjustments

	* 1. Run benchmark
	areg `lhs' `rhs', absorb(`absvars') cluster(`clustervar')
	TrimMatrix `K'
	storedresults save benchmark e()
	
	* 2. Run reghdfe
	* To match -areg- we need to have dof(none)!
	reghdfe `lhs' `rhs', absorb(`absvars') vce(cluster `clustervar') dof(none)
	TrimMatrix `K'

	* 3. Compare
	storedresults compare benchmark e(), tol(1e-12) include( ///
		scalar: N rmse tss rss r2 r2_a F df_r df_a df_m rmse ///  
		matrix: trim_b trim_V ///
		macros: wexp wtype )
	storedresults drop benchmark

// -------------------------------------------------------------------------------------------------
* Must match -xtreg- with the dof adjustments

	* 1. Run benchmark
	xtreg `lhs' `rhs', fe cluster(`clustervar')
	TrimMatrix `K'
	storedresults save benchmark e()
		di e(df_m)
	* 2. Run reghdfe
	* To match -areg- we need to have dof(none)!
	reghdfe `lhs' `rhs', absorb(`absvars') vce(cluster `clustervar')
	TrimMatrix `K'
	
	* 3. Compare
	* NOTE: I'm removing e(df_m) from the comparison b/c -xtreg,fe- acts crazy. See:
	* http://www.stata.com/statalist/archive/2010-03/msg00941.html
	* Also removing r2/r2a/df_a as expected, b/c xtreg has a different formula for the R2s
	storedresults compare benchmark e(), tol(1e-12) include( ///
		scalar: N tss rss F df_r /// df_m df_a r2 r2_a 
		matrix: trim_b trim_V ///
		macros: wexp wtype )
	storedresults drop benchmark

// -------------------------------------------------------------------------------------------------
cd "D:/Github/reghdfe/test"
exit
