noi cscript "reghdfe: ols with robust VCE" adofile reghdfe

* Dataset
	sysuse auto
	bys turn: gen t = _n
	tsset turn t
	drop if missing(rep)
	
* [TEST] Misc. options

	reghdfe price weight gear, absorb(turn trunk) vce(robust) /// verbose(-1) 
		level(90) tol(1e-8) accel(none) transf(cim) noprune
	assert r(level)==90
	assert e(drop_singletons)==1
	assert e(num_singletons)==10
	assert e(N)==59

	reghdfe price weight gear, absorb(turn trunk) vce(robust) /// verbose(-1) 
		level(90) tol(1e-8) accel(none) transf(cim) noprune keepsing
	assert r(level)==90
	assert e(drop_singletons)==0
	assert e(num_singletons)==0
	assert e(N)==69

	reghdfe price weight gear, absorb(turn trunk) tol(1e-5)
	loc ic = e(ic)

	reghdfe price weight gear, absorb(turn trunk) tol(1e-10)
	assert e(ic) < .
	assert `ic' < e(ic)

	reghdfe price weight gear, absorb(turn trunk) tol(1e-5) accel(none)
	assert e(ic) < .
	assert `ic' < e(ic)

	reghdfe price weight gear, absorb(turn trunk) tol(1e-5) transf(cim)
	assert e(ic) < .
	assert `ic' < e(ic)

	* Done!

exit
