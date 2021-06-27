noi cscript "reghdfe: ols with cluster VCE" adofile reghdfe


* Dataset
	sysuse auto
	bys turn: gen t = _n
	tsset turn t
	drop if missing(rep)

	
	local included_e ///
		scalar: N rmse tss rss r2 r2_a F df_r df_m ll ll_0 /// mss is not reported by areg if clustervar!=absvar
		matrix: b V ///
		macros: wexp wtype

* [TEST] Cluster

	local lhs price
	local rhs weight length
	local absvars turn
	local clustervar rep78
	fvunab tmp : `rhs'
	local K : list sizeof tmp
	
	* 1. Run benchmark
	areg `lhs' `rhs', absorb(`absvars') cluster(`clustervar')
	local bench_df_a = e(df_a)
	storedresults save benchmark e()
	
	* 2. Run reghdfe
	reghdfe `lhs' `rhs', absorb(`absvars') keepsingletons verbose(-1) vce(cluster `clustervar')
	storedresults compare benchmark e(), tol(1e-12) include(`included_e')
	assert `bench_df_a'==e(df_a)-1

	* 2. Run shorthand
	reghdfe `lhs' `rhs', absorb(`absvars') keepsingletons verbose(-1) cluster(`clustervar')
	storedresults compare benchmark e(), tol(1e-12) include(`included_e')
	assert `bench_df_a'==e(df_a)-1

	storedresults drop benchmark

* [TEST] Interacted Cluster

	local lhs price
	local rhs weight length
	local absvars turn
	local clustervar rep78#trunk
	fvunab tmp : `rhs'
	local K : list sizeof tmp
	
	* 1. Run benchmark
	egen rep_trunk = group(rep78 trunk)
	areg `lhs' `rhs', absorb(`absvars') cluster(rep_trunk)
	local bench_df_a = e(df_a)
	storedresults save benchmark e()
	
	* 2. Run reghdfe
	reghdfe `lhs' `rhs', absorb(`absvars') keepsingletons verbose(-1) vce(cluster `clustervar')
	storedresults compare benchmark e(), tol(1e-12) include(`included_e')
	assert `bench_df_a'==e(df_a)-1

	* 2. Run shorthand
	reghdfe `lhs' `rhs', absorb(`absvars') keepsingletons verbose(-1) cluster(`clustervar')
	storedresults compare benchmark e(), tol(1e-12) include(`included_e')
	assert `bench_df_a'==e(df_a)-1

	storedresults drop benchmark

exit
