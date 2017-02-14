noi cscript "reghdfe: fixed effect strictly nested within cluster" adofile reghdfe

* Dataset
	sysuse auto
	egen turn_trunk = group(turn trunk)
	bys turn_trunk: gen t = _n
	xtset turn_trunk t
	
	local included_e ///
		scalar: N tss rss F df_r ll ll_0 /// rmse r2 r2_a df_m 
		matrix: trim_b trim_V ///
		macros: wexp wtype

* [TEST] Cluster is absvar

	local lhs price
	local rhs weight length
	local absvars turn#trunk
	local clustervar turn
	fvunab tmp : `rhs'
	local K : list sizeof tmp

	* 1. Run benchmark
	xtreg `lhs' `rhs', fe cluster(`clustervar')
	trim_cons
	local bench_df_m = e(df_m)
	local bench_df_a = e(df_a)
	local bench_within = e(r2_w)
	storedresults save benchmark e()
	
	* 2. Run reghdfe
	reghdfe `lhs' `rhs', absorb(`absvars') keepsingletons verbose(-1) vce(cluster `clustervar')

	* NOTE: See statalist post for discussion on df_m discrepancy
	* http://www.stata.com/statalist/archive/2010-03/msg00941.html	
	* "So I think some explanation is necessary. I see no reason, conceptually, why xtreg,fe with small-sample statistics should not be exactly equivalent to areg"
	
	//loc adj1 = e(N) - e(rank) - e(df_a) //+ e(M_due_to_nested)
	//local adjdof = `adj1' / (`adj1'-1)
	//local adjdof 1
	//notrim ?? `adjdof'
	notrim

	assert e(df_a)==0
	assert `bench_df_m'==e(df_m)-1
	// BUGBUG
	// assert `bench_df_a'==e(df_a)+1
	
	storedresults compare benchmark e(), tol(1e-12) include(`included_e')

	* WILL NOT HOLD B/C xtreg still include the nested ones in e(df_a)
	* assert `bench_df_a'==e(df_a)-1

	assert abs(`bench_within'-e(r2_within))<1e-6

	storedresults drop benchmark

exit
