noi cscript "reghdfe: bug in st_data" adofile reghdfe

* Dataset
	sysuse auto
	bys turn: gen t = _n
	tsset turn t
	drop if missing(rep)
	
	local included_e ///
		scalar: N rmse tss rss mss r2 r2_a F df_r df_m ll ll_0 ///
		matrix: trim_b trim_V ///
		macros: wexp wtype

* [TEST] 1.x 2.x sometimes selects less variables

	local lhs price
	local rhs 32.turn 43.turn 51.turn
	local absvars foreign
	fvunab tmp : `rhs'
	local K : list sizeof tmp
	
	* 1. Run benchmark
	areg `lhs' 32bn.turn 43.turn 51.turn, absorb(`absvars')
	matrix list e(b)
	trim_cons
	local bench_df_a = e(df_a)
	storedresults save benchmark e()
	
	* 2. Run reghdfe
	reghdfe `lhs' `rhs', absorb(`absvars') keepsingletons verbose(-1)
	matrix list e(b)
	notrim
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
	gen double X1 = (foreign==0)*(turn==31)
	gen double X2 = (foreign==1)*(turn==32)
	areg `lhs' X1 X2, absorb(`absvars')
	matrix list e(b)
	trim_cons
	local bench_df_a = e(df_a)
	storedresults save benchmark e()
	
	* 2. Run reghdfe
	reghdfe `lhs' `rhs', absorb(`absvars') keepsingletons verbose(-1)
	matrix list e(b)
	notrim
	storedresults compare benchmark e(), tol(1e-12) include(`included_e')
	assert `bench_df_a'==e(df_a)-1

	* Done!
	storedresults drop benchmark

* [Test] rhs collinear
	reghdfe price 0.foreign#31.turn 1.foreign#32.turn, a(turn) keepsing
	matrix b = e(b)
	mata: assert(st_matrix("b")==J(1, 2, 0))
	matrix V = e(V)
	mata: assert(st_matrix("V")==J(2, 2, 0))

* [Test] rhs collinear and singletons dropped
* after the singletons, both rhs vars are all zero, so they will be omitted
	reghdfe price 0.foreign#31.turn 1.foreign#32.turn, a(turn)
	matrix b = e(b)
	mata: assert(st_matrix("b")==J(0, 0, 0))
	matrix V = e(V)
	mata: assert(st_matrix("V")==J(0, 0, 0))

exit
