*noi cscript "reghdfe: saving alphas (FEs)" adofile reghdfe

* ===========================================================================
* Tests estimates of fixed effects
* ===========================================================================
* need to test for each type of preconditioning... also intercept only, slope only, both, etc.

	loc techniques map lsmr // lsqr
	loc preconditioners none diag block_diag


// --------------------------------------------------------------------------
// Load dataset
// --------------------------------------------------------------------------

	sysuse auto, clear
	bys turn: gen t = _n
	tsset turn t

// --------------------------------------------------------------------------
// [TEST] Verify that alphas are not missing
// --------------------------------------------------------------------------

	reghdfe price, a(turn)
	gen byte sample = e(sample)
	keep if sample
	drop sample
	//sort turn
	gen double turn2 = turn


	di as text "{title:1) intercept (one-way-FE)}"
	foreach precond of local preconditioners {
		foreach tech of local techniques {
			di as text "[test alpha non-missing] technique={res}`tech'{txt} precond={res}`precond'"
			qui reghdfe turn2, a(FE=turn) keepsing v(-1) tech(`tech') precond(`precond')
			_assert abs(e(r2) - 1) < 1e-12 , msg("R2 assertion failed")
			gen double delta = turn - FE - _b[_cons]
			_assert abs(delta) < 1e-12 , msg("FE assertion failed")
			drop FE delta
		}
	}

	di as text "{title:1) intercept 2) slope}"
	foreach precond of local preconditioners {
		foreach tech of local techniques {
			di as text "[test alpha non-missing] technique={res}`tech'{txt} precond={res}`precond'"
			qui reghdfe turn2, a(turn turn#c.price, save) keepsing v(-1) tech(`tech') precond(`precond') tol(1e-12)
			gen double d = _b[_cons] + __hdfe1__ + __hdfe2__Slope1 * price
			_assert abs(e(r2) - 1) < 1e-12 , msg("R2 assertion failed")
			gen double delta = turn - d
			_assert abs(delta) < 1e-6 , msg("FE assertion failed") // lower tolerance (LSMR/LSQR seem less precise)
			drop delta d __*
		}
	}

	di as text "{title:1) joint intercept and slope}"
	foreach precond of local preconditioners {
		foreach tech of local techniques {
			di as text "[test alpha non-missing] technique={res}`tech'{txt} precond={res}`precond'"
			qui reghdfe turn2, a(turn##c.price, save) keepsing v(-1) tech(`tech') precond(`precond') tol(1e-12)
			gen double d = _b[_cons] + __hdfe1__ + __hdfe1__Slope1 * price
			_assert abs(e(r2) - 1) < 1e-12 , msg("R2 assertion failed")
			gen double delta = turn - d
			_assert abs(delta) < 1e-6 , msg("FE assertion failed") // lower tolerance (LSMR/LSQR seem less precise)
			drop delta d __*
		}
	}


	di as text "{title:1) joint intercept and slope 2) slope}"
	foreach precond of local preconditioners {
		foreach tech of local techniques {
			di as text "[test alpha non-missing] technique={res}`tech'{txt} precond={res}`precond'"
			qui reghdfe turn2, a(turn##c.price trunk#c.weight, save) keepsing v(-1) tech(`tech') precond(`precond') tol(1e-12)
			gen double d = _b[_cons] + __hdfe1__ + __hdfe1__Slope1 * price + __hdfe2__Slope1 * weight
			_assert abs(e(r2) - 1) < 1e-12 , msg("R2 assertion failed")
			gen double delta = turn - d
			_assert abs(delta) < 1e-6 , msg("FE assertion failed") // lower tolerance (LSMR/LSQR seem less precise)
			drop delta d __*
		}
	}
	drop turn2


	gen w1 = runiform()
	gen y = w1

	di as text "{title:1) joint intercept and slope; y explained by slope}"
	foreach precond of local preconditioners {
		foreach tech of local techniques {
			di as text "[test alpha non-missing] technique={res}`tech'{txt} precond={res}`precond'"
			qui reghdfe y, a(turn##c.w1, save) keepsing v(-1) tech(`tech') precond(`precond') tol(1e-12)
			gen double d = _b[_cons] + __hdfe1__ + __hdfe1__Slope1 * w1
			_assert abs(e(r2) - 1) < 1e-12 , msg("R2 assertion failed")
			gen double delta = y - d
			_assert abs(delta) < 1e-6 , msg("FE assertion failed") // lower tolerance (LSMR/LSQR seem less precise)
			drop delta d __*
		}
	}
	drop y w1


	gen w1 = runiform()
	gen w2 = runiform()
	gen y = w1 + 3*(turn<=40)*w2 - 8*(turn>40)*w2

	di as text "{title:1) joint intercept and two slopes; y explained by slope}"
	foreach precond of local preconditioners {
		foreach tech of local techniques {
			di as text "[test alpha non-missing] technique={res}`tech'{txt} precond={res}`precond'"
			qui reghdfe y, a(turn##c.(w1 w2), save) keepsing v(-1) tech(`tech') precond(`precond') tol(1e-12)
			gen double d = _b[_cons] + __hdfe1__ + __hdfe1__Slope1 * w1 + __hdfe1__Slope2 * w2
			_assert abs(e(r2) - 1) < 1e-12 , msg("R2 assertion failed")
			gen double delta = y - d
			_assert abs(delta) < 1e-6 , msg("FE assertion failed") // lower tolerance (LSMR/LSQR seem less precise)
			drop delta d __*
		}
	}
	drop y w1 w2



// --------------------------------------------------------------------------
// One FE
// --------------------------------------------------------------------------

	di as text "{title:Validate values of one FE}"
	qui areg price weight gear_ratio, absorb(turn)
	predict double D0 if e(sample), d

	foreach precond of local preconditioners {
		foreach tech of local techniques {

			di as text "[test alpha ok] technique={res}`tech'{txt} precond={res}`precond'"
			qui reghdfe price weight gear_ratio, absorb(FE=turn) keepsing v(-1) tech(`tech') precond(`precond')
			gen double D1 = FE

			qui reghdfe price weight gear_ratio, absorb(turn, save) keepsing v(-1) tech(`tech') precond(`precond')
			gen double D2 = __hdfe1__

			gen delta1 = reldif(D0, D1)
			gen delta2 = reldif(D0, D2)

			qui su D*

			assert(delta1 < 1e-10)
			assert(delta2 < 1e-10)

			drop D1 D2 delta1 delta2 __* FE
		}
	}
	drop D0


// --------------------------------------------------------------------------
// One FE, no regressors
// --------------------------------------------------------------------------
* SEE: https://github.com/sergiocorreia/reghdfe/issues/140

	qui areg price, absorb(turn)
	predict double D0 if e(sample), d

	di as text "{title:Validate values of one FE (no regressors)}"
	foreach precond of local preconditioners {
		foreach tech of local techniques {

			di as text "[test alpha ok] technique={res}`tech'{txt} precond={res}`precond'"

			qui reghdfe price, absorb(FE=turn) keepsing v(-1) tech(`tech') precond(`precond')
			gen double D1 = FE

			qui reghdfe price, absorb(turn, save) keepsing v(-1) tech(`tech') precond(`precond')
			gen double D2 = __hdfe1__

			gen delta1 = reldif(D0, D1)
			gen delta2 = reldif(D0, D2)

			assert(delta1 < 1e-10)
			assert(delta2 < 1e-10)
			drop D1 D2 delta1 delta2 __* FE
		}
	}
	
	drop D0


// --------------------------------------------------------------------------
// Several FEs, all should have mean zero
// --------------------------------------------------------------------------

	di as text "{title:Several FEs, check mean zero}"
	gen fw = ceil(runiform()*10)
	foreach precond of local preconditioners {
		foreach tech of local techniques {
			qui reghdfe price weight length gear_ratio, a(foreign turn trunk rep78, save) tol(1e-16) tech(`tech') precond(`precond')
			di as text "[test mean zero] technique={res}`tech'{txt} precond={res}`precond'"
			foreach var of varlist __hdfe* {
				qui su `var'
				assert(abs(r(mean)) < 1e-7) // this is actually good precision, as price has mean of 6,000
			}
			drop __hdfe*

			qui reghdfe price weight length gear_ratio [fw=fw], a(foreign turn trunk rep78, save) tol(1e-16) tech(`tech') precond(`precond')
			di as text "[test mean zero with fw] technique={res}`tech'{txt} precond={res}`precond'"
			foreach var of varlist __hdfe* {
				qui su `var' [fw=fw]
				assert(abs(r(mean)) < 1e-7) // this is actually good precision, as price has mean of 6,000
			}
			drop __hdfe*
		}
	}

// --------------------------------------------------------------------------
// One FE with MVs
// --------------------------------------------------------------------------

	qui areg price weight gear_ratio, absorb(rep78)
	loc c = _b[_cons]
	predict double bench_d if e(sample), d
	*replace bench_d = bench_d + `c' if e(sample)

	di as text "{title:One FE with MVs}"
	foreach precond of local preconditioners {
		foreach tech of local techniques {

			di as text "[test one fe with mvs] technique={res}`tech'{txt} precond={res}`precond'"
			qui reghdfe price weight gear_ratio, absorb(FE=rep78) keepsing v(-1) tech(`tech') precond(`precond')
			gen double d1 = FE

			qui reghdfe price weight gear_ratio, absorb(rep78, save) keepsing v(-1) tech(`tech') precond(`precond')
			gen double d2 = __hdfe1__

			gen double delta1 = reldif(bench_d, d1)
			gen double delta2 = reldif(bench_d, d2)

			assert(delta1 < 1e-9)
			assert(delta2 < 1e-9)

			drop d1 d2 delta1 delta2 __* FE
		}
	}
	
	drop bench_*


// --------------------------------------------------------------------------
// Two FEs
// --------------------------------------------------------------------------

	di as text "{title:Two FEs}"
	foreach precond of local preconditioners {
		foreach tech of local techniques {

			di as text "[test two fes] technique={res}`tech'{txt} precond={res}`precond'"

			qui reghdfe price weight length, tol(1e-10) a(v3_1=turn v3_2=trunk) version(3) nowarn
			qui reghdfe price weight length, tol(1e-10) a(v5_1=turn v5_2=trunk) version(5) nowarn
			qui reghdfe price weight length, tol(1e-10) a(fe1=turn fe2=trunk) v(-1) tech(`tech') precond(`precond')

			gen double delta_v5_1 = reldif(v5_1, fe1)
			gen double delta_v5_2 = reldif(v5_2, fe2)
			*su delta_*
			_assert(delta_v5_1 < 1e-6)
			_assert(delta_v5_2 < 1e-6)

			gen double delta_v3_1 = reldif(v3_1, fe1)
			gen double delta_v3_2 = reldif(v3_2, fe2)
			*su delta_*
			_assert(delta_v3_1 < 1e-6)
			_assert(delta_v3_2 < 1e-6)

			drop v?_? fe? delta_*
		}
	}	


// --------------------------------------------------------------------------
// Single slope
// --------------------------------------------------------------------------

	di as text "{title:Single slope}"
	foreach precond of local preconditioners {
		foreach tech of local techniques {

			di as text "[single slope] technique={res}`tech'{txt} precond={res}`precond'"

			qui reg price ibn.turn#c.gear_ratio, nocons
			predict double bench_slope, xb

			qui reghdfe price, a(turn#c.gear_ratio, save) keepsing tech(`tech') precond(`precond')
			gen double slope = __hdfe1__Slope1 * gear_ratio
			*su *slope
			gen double delta = reldif(bench_slope, slope)
			assert(delta < 1e-8)
			drop delta bench_slope slope
		}
	}



// --------------------------------------------------------------------------
// Several slopes
// --------------------------------------------------------------------------

	di as text "{title:Several slopes}"
	foreach precond of local preconditioners {
		foreach tech of local techniques {

			di as text "[several slopes] technique={res}`tech'{txt} precond={res}`precond'"
			qui areg price ibn.turn#c.(gear_ratio disp), a(turn)
			predict double benchmark, xbd

			qui reghdfe price, a(turn##c.(gear_ratio disp), save) keepsing v(-1) tol(1e-12) tech(`tech') precond(`precond') resid
			predict double xbd1, xbd
			gen double xbd2 = _b[_cons] + __hdfe1__ + __hdfe1__Slope1 * gear_ratio + __hdfe1__Slope2 * disp

			***gen double delta = reldif(bench_slopes, slopes)
			***assert(delta <  1e-8 | resid==0)
			// if the residuals are =0 , areg chooses the intercepts first, so it might pick different combinations of the alphas
			
			gen double delta1 = reldif(xbd1, benchmark)
			gen double delta2 = reldif(xbd2, benchmark)

			assert delta1 < 1e-6
			assert delta2 < 1e-6

			drop benchmark xbd* delta* __*
		}
	}

// --------------------------------------------------------------------------
// Several slopes without constant
// --------------------------------------------------------------------------

	di as text "{title:Several slopes}"
	foreach precond of local preconditioners {
		foreach tech of local techniques {

			di as text "[several slopes] technique={res}`tech'{txt} precond={res}`precond'"
			qui reg price ibn.turn#c.(gear_ratio disp), nocons
			predict double benchmark, xb

			qui reghdfe price, a(turn#c.(gear_ratio disp), save) keepsing v(-1) tol(1e-12) tech(`tech') precond(`precond') resid
			predict double xbd1, xbd
			gen double xbd2 = __hdfe1__Slope1 * gear_ratio + __hdfe1__Slope2 * disp

			***gen double delta = reldif(bench_slopes, slopes)
			***assert(delta <  1e-8 | resid==0)
			// if the residuals are =0 , areg chooses the intercepts first, so it might pick different combinations of the alphas
			
			gen double delta1 = reldif(xbd1, benchmark)
			gen double delta2 = reldif(xbd2, benchmark)

			assert delta1 < 1e-6
			assert delta2 < 1e-6

			drop benchmark xbd* delta* __*
		}
	}


exit
