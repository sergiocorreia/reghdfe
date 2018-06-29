noi cscript "reghdfe postestimation: predict with slope FEs (id##c.var)" adofile reghdfe

* Dataset
	sysuse auto
	bys turn: gen t = _n
	tsset turn t
	
	local included_e ///
		scalar: N rmse tss rss mss r2 r2_a df_r ll /// F df_m ll_0 
		matrix: trim_b trim_V ///
		macros: wexp wtype

* [TEST] predict after reghdfe

	local lhs price
	local rhs weight length gear disp
	local absvars turn
	fvunab tmp : `rhs'
	local K : list sizeof tmp

	* 1. Run benchmark
	areg `lhs' `rhs' ibn.turn#c.mpg , absorb(`absvars')
	di e(df_a)
	trim_cons `K'
	storedresults save benchmark e()
	predict double xbd, xbd
	predict double resid, resid

	* 2. Run reghdfe and compare
	reghdfe `lhs' `rhs', absorb(`absvars'##c.mpg) keepsingletons resid verbose(-1)
	trim_cons
	predict double xbd_test, xbd
	predict double resid_test, resid
	su xbd xbd_test resid resid_test, sep(2)
	storedresults compare benchmark e(), tol(1e-12) include(`included_e')

	_vassert xbd xbd_test, tol(1e-12)
	_vassert resid resid_test, tol(1e-11) // Decreased tolerance!

	* Test that we can't predict resid after dropping e(resid)
	drop *_test `e(resid)'
	cap predict resid_test, resid
	assert c(rc)

storedresults drop benchmark
exit
