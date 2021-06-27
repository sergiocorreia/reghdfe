noi cscript "reghdfe: ols with if/in" adofile reghdfe

* Dataset
	sysuse auto
	bys turn: gen t = _n
	tsset turn t
	drop if missing(rep)
	
	local included_e ///
		scalar: N rmse tss rss mss r2 r2_a F df_r df_m ll ll_0 ///
		matrix: b V ///
		macros: wexp wtype

* [TEST] Undjusted

	local lhs price
	local rhs weight length
	local absvars turn
	local cond if !foreign in 10/60
	fvunab tmp : `rhs'
	local K : list sizeof tmp
	
	* 1. Run benchmark
	areg `lhs' `rhs' `cond', absorb(`absvars')
	local bench_df_a = e(df_a)
	gen byte bench_sample = e(sample)
	storedresults save benchmark e()
	
	* 2. Run reghdfe
	reghdfe `lhs' `rhs' `cond', absorb(`absvars') keepsingletons verbose(-1)
	storedresults compare benchmark e(), tol(1e-12) include(`included_e')
	gen byte sample = e(sample)
	assert `bench_df_a'==e(df_a)-1
	assert bench_sample==sample

	* Done!

storedresults drop benchmark
exit
