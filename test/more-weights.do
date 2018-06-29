cscript "reghdfe with weights (corner cases)" adofile reghdfe

* Dataset
	sysuse auto
	gen n = int(uniform()*10000000+3)
	// drop if missing(rep)
	
	local included_e ///
		scalar: N rmse tss rss mss r2 r2_a F df_r df_m ll ll_0 ///
		matrix: b V ///
		macros: wexp wtype


* [TEST] Freq. weight with huge N (to test numeric stability)

	local lhs price
	local rhs weight length gear disp
	local absvars turn
	local wvar n
	local wtype fw
	fvunab tmp : `rhs'
	local K : list sizeof tmp


	* 1. Run benchmark
	areg `lhs' `rhs' [`wtype'=`wvar'], absorb(`absvars')
	di e(df_a)
	storedresults save benchmark e()


	reghdfe `lhs' `rhs' [`wtype'=`wvar'], absorb(`absvars') keepsingletons resid verbose(-1)
	storedresults compare benchmark e(), tol(1e-11) include(`included_e')
	storedresults drop benchmark

* [TEST] [aw] should fail with negative weights
	gen float m = 1
	replace m = -2 in 1/2
	rcof "reghdfe price weight disp [aw=m], a(turn#foreign)" == 402 // negative weights
	drop m


* [TEST] [fw] should fail with non-integer weights
	gen float m = 1 + uniform()
	rcof "reghdfe price weight disp [fw=m], a(turn#foreign)" == 401
	drop m

exit
