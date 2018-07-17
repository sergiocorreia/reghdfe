noi cscript "reghdfe: bug in st_data" adofile reghdfe

* Dataset
	sysuse auto
	bys turn: gen t = _n
	tsset turn t
	drop if missing(rep)
	
	local included_e ///
		scalar: N rmse tss rss mss r2 r2_a F df_r df_m ll ll_0 ///
		matrix: b V ///
		macros: wexp wtype

* [TEST] 1.x 2.x sometimes selects less variables

	local lhs price
	local rhs 32bn.turn 43.turn 51.turn
	local absvars foreign
	fvunab tmp : `rhs'
	local K : list sizeof tmp
	
	* 1. Run benchmark
	areg `lhs' `rhs', absorb(`absvars')
	matrix list e(b)
	local bench_df_a = e(df_a)
	storedresults save benchmark e()
	
	* 2. Run reghdfe
	reghdfe `lhs' `rhs', absorb(`absvars') keepsingletons verbose(-1)
	matrix list e(b)
	storedresults compare benchmark e(), tol(1e-12) include(`included_e')
	assert `bench_df_a'==e(df_a)-1

	* Done!
	storedresults drop benchmark

* [TEST] Now with a base level

	local lhs price
	local rhs 32.turn 43.turn 51.turn
	local absvars foreign
	fvunab tmp : `rhs'
	local K : list sizeof tmp
	
	* 1. Run benchmark
	areg `lhs' 32.turn 43.turn 51.turn, absorb(`absvars')
	matrix list e(b)
	local bench_df_a = e(df_a)
	storedresults save benchmark e()
	
	* 2. Run reghdfe
	reghdfe `lhs' `rhs', absorb(`absvars') keepsingletons verbose(-1)
	matrix list e(b)
	storedresults compare benchmark e(), tol(1e-12) include(`included_e')
	assert `bench_df_a'==e(df_a)-1

	* Done!
	storedresults drop benchmark


* [TEST] Interactions fail

	local lhs price
	local rhs 0.foreign#31.turn 1.foreign#32.turn
	local absvars trunk
	fvunab tmp : `rhs'
	local K : list sizeof tmp
	
	* 1. Run benchmark
	areg `lhs' `rhs', absorb(`absvars')
	matrix list e(b)
	local bench_df_a = e(df_a)
	storedresults save benchmark e()
	
	* 2. Run reghdfe
	reghdfe `lhs' `rhs', absorb(`absvars') keepsingletons verbose(-1)
	matrix list e(b)
	storedresults compare benchmark e(), tol(1e-12) include(`included_e')
	assert `bench_df_a'==e(df_a)-1

	* Done!
	storedresults drop benchmark

* [Test] rhs collinear
	loc vars price 0.foreign#31.turn 1.foreign#32.turn

	local included_e ///
		scalar: N rmse tss rss mss r2 r2_a df_r df_m ll ll_0 /// F
		matrix: b V ///
		macros: wexp wtype

	* 1. Run benchmark
	areg `vars', a(turn)
	matrix list e(b)
	local bench_df_a = e(df_a)
	assert e(F)==.
	storedresults save benchmark e()

	* 2. Run reghdfe
	reghdfe `vars', a(turn) keepsing v(-1)
	matrix b = e(b)
	matrix list b
	assert e(F)==.
	storedresults compare benchmark e(), tol(1e-12) include(`included_e')
	assert `bench_df_a'==e(df_a)-1

	* Done!
	storedresults drop benchmark

* [Test] rhs collinear and singletons dropped
* after the singletons, both rhs vars are all zero, so they will be omitted

	loc vars price 0.foreign#31.turn 1.foreign#32.turn

	* 0. Get sample
	reghdfe `vars', a(turn) keepsing v(-1)
	gen byte sample = e(sample)

	* 1. Run benchmark
	areg `vars' if sample, a(turn)
	matrix list e(b)
	local bench_df_a = e(df_a)
	assert e(F)==.
	storedresults save benchmark e()

	* 2. Run reghdfe
	reghdfe `vars', a(turn) keepsing v(-1)
	matrix b = e(b)
	matrix list b
	assert e(F)==.
	storedresults compare benchmark e(), tol(1e-12) include(`included_e')
	assert `bench_df_a'==e(df_a)-1

	* Done!
	storedresults drop benchmark

exit
