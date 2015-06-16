cscript "reghdfe with clusters" adofile reghdfe

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

* [TEST] Cluster
	local lhs price
	local rhs weight length
	local absvars turn
	local clustervar rep
	fvunab tmp : `rhs'
	local K : list sizeof tmp

	drop if missing(rep)

	* 1. Run benchmark
	areg `lhs' `rhs', absorb(`absvars') cluster(`clustervar')
	di as result "areg `lhs' `rhs', absorb(`absvars') cluster(`clustervar')"
	matrix list e(V)
	TrimMatrix `K'
	local bench_df_a = e(df_a)
	storedresults save benchmark e()
	
	* 2. Run reghdfe
	di as result "reghdfe `lhs' `rhs', absorb(`absvars') vce(cluster `clustervar') keepsingletons"
	reghdfe `lhs' `rhs', absorb(`absvars') vce(cluster `clustervar') keepsingletons // dof(none)	
	matrix list e(V)
	TrimMatrix `K'
	
	* 3. Compare
	storedresults compare benchmark e(), tol(1e-12) include( ///
		scalar: N rmse tss rss r2 r2_a F df_r df_m ll ll_0 /// F_absorb 
		matrix: trim_b trim_V ///
		macros: wexp wtype )
	storedresults drop benchmark
	* NOTE: What should I use to build F_absorb in this case?
	assert `bench_df_a'==e(df_a)-1

* [TEST] Interacted cluster
	local lhs price
	local rhs weight length
	local absvars turn
	local clustervar turn#trunk
	fvunab tmp : `rhs'
	local K : list sizeof tmp

	drop if missing(rep)

	* 1. Run benchmark
	egen turn_trunk = group(turn trunk)
	areg `lhs' `rhs', absorb(`absvars') cluster(turn_trunk)
	TrimMatrix `K'
	local bench_df_a = e(df_a)
	storedresults save benchmark e()
	
	* 2. Run reghdfe
	reghdfe `lhs' `rhs', absorb(`absvars') vce(cluster `clustervar') keepsingletons // dof(none)
	TrimMatrix `K'
	
	* 3. Compare
	storedresults compare benchmark e(), tol(1e-12) include( ///
		scalar: N rmse tss rss r2 r2_a F df_r df_m ll ll_0 /// F_absorb 
		matrix: trim_b trim_V ///
		macros: wexp wtype )
	storedresults drop benchmark
	assert `bench_df_a'==e(df_a)-1

	
cd "D:/Github/reghdfe/test"
exit
