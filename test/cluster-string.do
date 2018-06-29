noi cscript "reghdfe: ols with a string cluster" adofile reghdfe


* Dataset
	sysuse auto
	bys turn: gen t = _n
	tsset turn t
	drop if missing(rep)
	gen REP = "X" + strofreal(rep) if foreign
	gen xrep = rep if foreign

	
	local included_e ///
		scalar: N rmse tss rss r2 r2_a F df_r df_m ll ll_0 /// mss is not reported by areg if clustervar!=absvar
		matrix: b V ///
		macros: wexp wtype

* [TEST] Cluster

	local lhs price
	local rhs weight length
	local absvars turn
	local clustervar REP
	fvunab tmp : `rhs'
	local K : list sizeof tmp
	
	* 1. Run benchmark
	areg `lhs' `rhs', absorb(`absvars') cluster(`clustervar')
	local bench_df_a = e(df_a)
	storedresults save benchmark e()
	
	* 2. Run reghdfe
	//cou
	//su `lhs' `rhs' `absvars'
	//tab `clustervar'
	reghdfe `lhs' `rhs', absorb(`absvars') keepsingletons verbose(-1) vce(cluster `clustervar')
	storedresults compare benchmark e(), tol(1e-11) include(`included_e')
	assert `bench_df_a'==e(df_a)-1

	storedresults drop benchmark

exit
