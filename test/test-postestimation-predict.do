noi cscript "reghdfe postestimation: predict" adofile reghdfe

* Setup
	discard
	clear all
	set more off
	* cls

	local included_e scalar: N rmse tss rss r2 r2_a F df_r df_a df_m ///
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
	TrimMatrix `K'
	storedresults save bench_e e()
	predict double xb, xb
	predict double d, d
	predict double xbd, xbd
	predict double resid, resid

	* 2. Run reghdfe and compare
	
	* 2a) With constant
	reghdfe `lhs' `rhs', absorb(FE=`absvars') vce(, suite(default))
	TrimMatrix `K'
	storedresults compare bench_e e(), tol(1e-10) include(`included_e') // Note the lowered tol
	predict double xb_test, xb
	predict double d_test, d
	predict double xbd_test, xbd
	predict double resid_test, resid
	*su d d_test xb xb_test xbd xbd_test resid resid_test, sep(2)
	
	_vassert xb xb_test, tol(1e-12)
	_vassert d d_test, tol(1e-12)
	_vassert xbd xbd_test, tol(1e-12)
	_vassert resid resid_test, tol(1e-12)
	
	drop *_test FE

	* 2a) Without constant
	reghdfe `lhs' `rhs', absorb(FE=`absvars') vce(, suite(default)) nocons
	TrimMatrix `K'
	storedresults compare bench_e e(), tol(1e-10) include(`included_e') // Note the lowered tol
	predict double xb_test, xb
	predict double d_test, d
	predict double xbd_test, xbd
	predict double resid_test, resid
	*su d d_test xb xb_test xbd xbd_test resid resid_test, sep(2)
	
	_vassert xb xb_test, tol(1e-12)
	_vassert d d_test, tol(1e-12)
	_vassert xbd xbd_test, tol(1e-12)
	_vassert resid resid_test, tol(1e-12)
	
	drop *_test FE

	storedresults drop bench_e

	* 3) Check that we can run with score
	qui reghdfe `lhs' `rhs', absorb(FE=`absvars')
	predict score*, score
	predict alt, resid
	su score* alt
	_vassert alt score, tol(1e-12)
	drop FE alt resid


cd "D:/Github/reghdfe/test"
exit
