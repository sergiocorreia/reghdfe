noi cscript "reghdfe postestimation: pweight with cluster" adofile reghdfe

* Dataset
	sysuse auto
	bys turn: gen t = _n
	tsset turn t
	
	if (c(version)>=14) loc mss mss
	
	local included_e1 ///
		scalar: N rmse tss rss `mss' r2 r2_a df_r df_m ll ll_0 /// F 
		matrix: b /// trim_V
		macros: wexp wtype

	
	local included_e2 ///
		scalar: N tss rss df_r ll ll_0 F /// rmse mss r2 r2_a df_m
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
	areg `lhs' `rhs' [`wtype'=`wvar'], absorb(`absvars') vce(cluster turn)
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
	
	reghdfe `lhs' `rhs' [`wtype'=`wvar'], absorb(`absvars') keepsingletons resid verbose(-1) vce(cluster turn)
	assert e(df_a)==0
	assert `bench_df_a'==e(df_a)+e(df_a_nested)-1
	predict double xb_test, xb
	predict double d_test, d
	predict double xbd_test, xbd
	predict double resid_test, resid
	predict double dr_test, dr
	predict double stdp_test, stdp
	su d d_test xb xb_test xbd xbd_test resid resid_test dr dr_test stdp stdp_test, sep(2)
	storedresults compare benchmark e(), tol(1e-10) include(`included_e1')

	_vassert xb xb_test, tol(1e-10)
	_vassert d d_test, tol(1e-10)
	_vassert xbd xbd_test, tol(1e-10)
	_vassert resid resid_test, tol(1e-10)
	_vassert dr dr_test, tol(1e-10)

storedresults drop benchmark


* Bench
	xtreg price weight length [pw=turn] , fe robust
	trim_cons
	storedresults save benchmark e()

* reghdfe
	reghdfe price weight length [pw=turn] , a(turn) cluster(turn) keepsing v(-1)
	notrim
	storedresults compare benchmark e(), tol(1e-10) include(`included_e2')
exit
