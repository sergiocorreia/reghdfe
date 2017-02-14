noi cscript "reghdfe postestimation: predict after pweight" adofile reghdfe

* Dataset
	sysuse auto
	bys turn: gen t = _n
	tsset turn t
	// drop if missing(rep)
	
	local included_e ///
		scalar: N rmse tss rss mss r2 r2_a F df_r df_m ll ll_0 ///
		matrix: trim_b trim_V ///
		macros: wexp wtype

* [TEST] predict after reghdfe

	local lhs price
	local rhs weight length gear disp
	local absvars turn
	local wvar turn
	local wtype fw
	fvunab tmp : `rhs'
	local K : list sizeof tmp

	* 1. Run benchmark
	areg `lhs' `rhs' [`wtype'=`wvar'], absorb(`absvars')
	di e(df_a)
	trim_cons
	local bench_df_a = e(df_a)
	storedresults save benchmark e()
	predict double xb, xb
	predict double d, d
	predict double xbd, xbd
	predict double resid, resid
	predict double dr, dr
	predict double stdp, stdp

	* AREG and REGHDFE disagree because AREG includes _cons in XB instead of D
	replace xb = xb - _b[_cons]
	replace d = d + _b[_cons]
	replace dr = dr + _b[_cons]
	replace stdp = stdp
	su resid, mean

	* 2. Run reghdfe and compare
	
	reghdfe `lhs' `rhs' [`wtype'=`wvar'], absorb(`absvars') keepsingletons resid verbose(-1)
	notrim
	assert `bench_df_a'==e(df_a)-1
	predict double xb_test, xb
	predict double d_test, d
	predict double xbd_test, xbd
	predict double resid_test, resid
	predict double dr_test, dr
	predict double stdp_test, stdp
	su d d_test xb xb_test xbd xbd_test resid resid_test dr dr_test stdp stdp_test, sep(2)
	storedresults compare benchmark e(), tol(1e-11) include(`included_e')

	_vassert xb xb_test, tol(1e-12)
	_vassert d d_test, tol(1e-12)
	_vassert xbd xbd_test, tol(1e-12)
	_vassert resid resid_test, tol(1e-12)
	_vassert dr dr_test, tol(1e-12)

storedresults drop benchmark
exit
