noi cscript "reghdfe with only slope effects" adofile reghdfe

* Setup
	discard
	clear all
	set more off
	* cls

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
	drop if missing(rep)
	
	local included_e scalar: N rmse rss r2 r2_a F df_r ll /// tss F_absorb ll_0 df_m 
		matrix: trim_b trim_V ///
		macros: wexp wtype

* [TEST] Undjusted
	local lhs price
	local rhs weight length
	local absvars i.turn#c.gear
	fvunab tmp : `rhs'
	local K : list sizeof tmp
	
	* 1. Run benchmark
	regress `lhs' `rhs' `absvars', nocons
	predict double d, xb
	foreach var of local rhs {
		replace d = d - _b[`var']*`var'
	}
	rename d d_benchmark
	TrimMatrix `K'
	test `rhs'
	estadd scalar F = r(F), replace
	storedresults save benchmark e()

	* 2. Run reghdfe and compare
	reghdfe `lhs' `rhs', absorb(`absvars', save) keepsingletons
	predict d_reghdfe, d
	TrimMatrix `K'
	storedresults compare benchmark e(), tol(1e-12) include(`included_e')

	su d_*
	corr d_*
	tabstat __hdfe, by(turn)
	assert r(rho)>=0.999

	reghdfe `lhs', absorb(`absvars', save) keepsingletons
	predict xb_reghdfe, xb
	assert abs(xb)<0.0001

	storedresults drop benchmark

cd "D:/Github/reghdfe/test"
exit
