noi cscript "reghdfe postestimation: predict" adofile reghdfe

* Setup
	discard
	clear all
	set more off
	* cls

	local included_e scalar: N rmse tss rss r2 r2_a F df_r df_m ///
		matrix: trim_b trim_V ///
		macros: wexp wtype

* Convenience: "Trim <size>" will trim e(b) and e(V)
	capture program drop TrimMatrix
	program define TrimMatrix, eclass
	args size
		assert `size'>0
		matrix trim_b = e(b)
		matrix trim_V = e(V)
		matrix trim_b = trim_b[1, 1..`size']
		matrix trim_V = trim_V[1..`size',1..`size']
		ereturn matrix trim_b = trim_b
		ereturn matrix trim_V = trim_V
	end
set trace off	
* Create fake dataset
	sysuse auto
	*gen n = int(uniform()*10+3) // used for weights
	*replace length = 0 if rep==3 // used for DoF adjustment of cont var
	*replace length = 5 if rep==1
	*gen byte one = 1
	bys turn: gen t = _n
	tsset turn t

	
	noi di as text " -  description of test"
	local lhs price
	local rhs weight length gear disp
	local absvars turn
	fvunab tmp : `rhs'
	local K : list sizeof tmp

	local testvars length gear

* [TEST]

	* 1. Run benchmark
	areg `lhs' `rhs', absorb(`absvars')
	di e(df_a)
	TrimMatrix `K'
	local bench_df_a = e(df_a)
	storedresults save bench_e e()
	predict double xb, xb
	predict double d, d
	predict double xbd, xbd
	predict double resid, resid

	* DIFFERENCE BETWEEN METHODS, XB HAS CONSTANT IN AREG BUT NOT IN REGHDFE
	replace xb = xb - _b[_cons]
	replace d = d + _b[_cons]
	su resid, mean

	* 2. Run reghdfe and compare
	
	reghdfe `lhs' `rhs', absorb(FE=`absvars') vce(, suite(default)) keepsingletons
	TrimMatrix `K'
	assert `bench_df_a'==e(df_a)-1
	predict double xb_test, xb
	predict double d_test, d
	predict double xbd_test, xbd
	predict double resid_test, resid
	su d d_test xb xb_test xbd xbd_test resid resid_test, sep(2)
	storedresults compare bench_e e(), tol(1e-10) include(`included_e') // Note the lowered tol

	_vassert xb xb_test, tol(1e-10)
	_vassert d d_test, tol(1e-10)
	_vassert xbd xbd_test, tol(1e-10)
	_vassert resid resid_test, tol(1e-10)
	
	drop *_test FE
	storedresults drop bench_e

	* 3) Check that we can run with score
	qui reghdfe `lhs' `rhs', absorb(FE=`absvars')
	predict score*, score
	predict alt, resid
	su score* alt
	_vassert alt score, tol(1e-10)
	drop FE alt resid

	* 4) Test that the means of resid are zero , with #c. interactions
	// di as error "WARNING: Note the numerical inaccuracy of interacting with c.vars" // this is for step 5)
	drop d xbd xb 
	reghdfe price weight, a(turn trunk##c.gear, save) tol(1e-12) maxiter(1e4) keepsingletons v(1)
	predict double resid, resid
	predict double d, d
	predict double xbd, xbd
	su d xbd
	su resid
	assert abs(r(mean))<1e-8

	* 5) Check that the means match
	areg price weight ibn.trunk##c.gear, a(turn)
	predict double resid_bench, resid
	corr resid resid_bench
	su resid*
	gen double delta = reldif(resid_bench , resid) // / abs(price)
	su delta
	*_vassert resid resid_bench, tol(1e-10)
	assert delta<1e-4 //  BUGBUG need something better than SD+KAC for the resids
	di as error "BUGBUG need something better than SD+KAC for the resids!"
	drop resid resid_bench
	* BUGBUG:

cd "C:/Git/reghdfe/test"
exit
