noi cscript "reghdfe: saving alphas (FEs)" adofile reghdfe

* Dataset
	sysuse auto
	bys turn: gen t = _n
	tsset turn t

* [TEST] Non-missing alphas
	gen double turn2 = turn
	reghdfe turn2, a(TURN=turn) keepsing v(-1)
	assert abs(TURN - turn) < 1e-12
	drop turn2 TURN

* [TEST] One FE
	areg price weight gear, absorb(turn)
	loc c = _b[_cons]
	predict double bench_d if e(sample), d
	replace bench_d = bench_d + `c' if e(sample)

	reghdfe price weight gear, absorb(FE=turn) keepsing v(-1)
	gen double d1 = FE

	reghdfe price weight gear, absorb(turn, save) keepsing v(-1)
	gen double d2 = __hdfe1__

	gen delta1 = reldif(bench_d, d1)
	gen delta2 = reldif(bench_d, d2)

	assert(delta1 < 1e-10)
	assert(delta2 < 1e-10)

	drop bench_* delta* d1 d2 FE __*

* [TEST] One FE with MVs
	areg price weight gear, absorb(rep)
	loc c = _b[_cons]
	predict double bench_d if e(sample), d
	replace bench_d = bench_d + `c' if e(sample)

	reghdfe price weight gear, absorb(FE=rep) keepsing v(-1)
	gen double d1 = FE

	reghdfe price weight gear, absorb(rep, save) keepsing v(-1)
	gen double d2 = __hdfe1__

	gen double delta1 = reldif(bench_d, d1)
	gen double delta2 = reldif(bench_d, d2)

	assert(delta1 < 1e-9)
	assert(delta2 < 1e-9)

	drop bench_* delta* d1 d2 FE __*


* [TEST] Two FEs
	reghdfe price weight length, a(bench1=turn bench2=trunk) old

	reghdfe price weight length, a(fe1=turn fe2=trunk) v(-1)
	su fe1 if e(sample)
	replace fe1 = fe1-r(mean)
	su bench* fe*

	gen double delta1 = reldif(bench1, fe1)
	gen double delta2 = reldif(bench2, fe2)

	assert(delta1 < 1e-8)
	assert(delta2 < 1e-8)

	drop bench* delta* fe*

* [TEST] Slope
	reg price ibn.turn#c.gear, nocons
	predict double bench_slope, xb

	reghdfe price, a(turn#c.gear, save) keepsing // resid
	gen double slope = __hdfe1__Slope1 * gear
	su *slope
	gen double delta = reldif(bench_slope, slope)
	assert(delta < 1e-8)
	drop delta bench_slope slope

* [TEST] Several slopes
	areg price ibn.turn#c.(gear disp), a(turn)
	predict double bench_slopes, xb
	predict double bench_intercept, d
	replace bench_slopes = bench_slopes - _b[_cons]
	predict double resid, resid

	reghdfe price, a(turn##c.(gear disp), save) keepsing v(-1) tol(1e-12)
	su __hdfe*
	gen double slopes = __hdfe1__Slope1 * gear + __hdfe1__Slope2 * disp

	gen double delta = reldif(bench_slopes, slopes)
	assert(delta <  1e-8 | resid==0)
	// if the residuals are =0 , areg chooses the intercepts first, so it might pick different combinations of the alphas

	drop resid bench_slopes bench_intercept slopes __hdfe* delta


* [TEST] Slope but no constant
	reg price ibn.turn#c.(gear disp), nocons
	predict double bench_slopes, xb
	predict double resid, resid

	reghdfe price, a(turn#c.(gear disp), save) keepsing v(-1) tol(1e-12)
	su __hdfe*
	gen double slopes = __hdfe1__Slope1 * gear + __hdfe1__Slope2 * disp

	gen double delta = reldif(bench_slopes, slopes)
	assert(delta <  1e-8)

	drop resid bench_slopes slopes __hdfe* delta

exit
