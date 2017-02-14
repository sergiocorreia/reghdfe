noi cscript "reghdfe: ols with unadjusted VCE" adofile reghdfe

* Dataset
	sysuse auto
	bys turn: gen t = _n
	tsset turn t
	drop if missing(rep)
	
	local included_e ///
		scalar: N rmse tss rss mss r2 r2_a F df_r df_m ll ll_0 ///
		matrix: trim_b trim_V ///
		macros: wexp wtype

* [TEST] Unadjusted

	local lhs price
	local rhs weight length
	local absvars turn
	fvunab tmp : `rhs'
	local K : list sizeof tmp
	
	* 1. Run benchmark
	areg `lhs' `rhs', absorb(`absvars')
	trim_cons
	local bench_df_a = e(df_a)
	storedresults save benchmark e()
	
	* 2. Run reghdfe
	reghdfe `lhs' `rhs', absorb(`absvars') keepsingletons verbose(-1)
	notrim
	storedresults compare benchmark e(), tol(1e-12) include(`included_e')
	assert `bench_df_a'==e(df_a)-1

	* 3. Run reghdfe vce(ols)
	reghdfe `lhs' `rhs', absorb(`absvars') keepsingletons verbose(-1) vce(ols)
	notrim
	storedresults compare benchmark e(), tol(1e-12) include(`included_e')
	assert `bench_df_a'==e(df_a)-1

	* 4. Run reghdfe vce(unadjusted)
	reghdfe `lhs' `rhs', absorb(`absvars') keepsingletons verbose(-1) vce(unadjusted)
	notrim
	storedresults compare benchmark e(), tol(1e-12) include(`included_e')
	assert `bench_df_a'==e(df_a)-1

	* Done!

storedresults drop benchmark
exit
