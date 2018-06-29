noi cscript "reghdfe postestimation: predict after pweight" adofile reghdfe

* Dataset
	sysuse auto
	bys turn: gen t = _n
	tsset turn t
	// drop if missing(rep)

	if (c(version)>=14) loc mss mss
	
	local included_e ///
		scalar: N rmse tss rss `mss' r2 r2_a F df_r df_m ll ll_0 ///
		matrix: b V ///
		macros: wexp wtype

* [TEST] predict after reghdfe

	local lhs price
	local rhs weight length gear disp
	local absvars turn
	local wvar turn
	local wtype pw
	fvunab tmp : `rhs'
	local K : list sizeof tmp

	* 1. Run benchmark
	areg `lhs' `rhs' [`wtype'=`wvar'], absorb(`absvars')
	di e(df_a)
	local bench_df_a = e(df_a)
	storedresults save benchmark e()
	predict double xb, xb
	predict double d, d
	predict double xbd, xbd
	predict double resid, resid
	predict double dr, dr
	predict double stdp, stdp

	* 2. Run reghdfe and compare
	reghdfe `lhs' `rhs' [`wtype'=`wvar'], absorb(`absvars') keepsingletons resid verbose(-1)
	assert `bench_df_a'==e(df_a)-1
	predict double xb_test, xb
	predict double d_test, d
	predict double xbd_test, xbd
	predict double resid_test, resid
	predict double dr_test, dr
	predict double stdp_test, stdp
	su d d_test xb xb_test xbd xbd_test resid resid_test dr dr_test stdp stdp_test, sep(2)
	storedresults compare benchmark e(), tol(1e-10) include(`included_e')

	_vassert xb xb_test, tol(1e-10)
	_vassert d d_test, tol(1e-10)
	_vassert xbd xbd_test, tol(1e-10)
	_vassert resid resid_test, tol(1e-10)
	_vassert dr dr_test, tol(1e-10)

storedresults drop benchmark
exit
