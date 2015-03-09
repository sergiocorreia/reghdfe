cd "D:/Github/reghdfe" // /source
cscript "reghdfe with HAC VCE (bw>1)" adofile reghdfe

* Setup
	discard
	pr drop _all
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
	set seed 43535
	bys turn (rep mpg price): gen t = _n
	tsset turn t

* [TEST] Simple example
	
	noi di as text " -  VCE with absvars==clustervars==tsset, and bw(2)"
	local lhs price
	local rhs weight length gear disp
	local absvars turn t
	fvunab tmp : `rhs'
	local K : list sizeof tmp

	* 1. Run benchmark
	qui tab turn, gen(TURN_)
	qui tab t, gen(T_)
	cap ivreg2 `lhs' `rhs' TURN_* T_*, small cluster(`absvars') bw(2) partial(TURN_* T_*)
	assert _rc==0
	TrimMatrix `K'
	storedresults save benchmark e()

	* 1. Run ANOTHER benchmark just for the R2
	cap ivreg2 `lhs' `rhs' TURN_* T_*, small cluster(`absvars') bw(2) // partial(TURN_* T_*)
	assert _rc==0
	TrimMatrix `K'
	storedresults save benchmark2 e()	
	
	* 2. Run reghdfe
	reghdfe `lhs' `rhs', a(`absvars') vce(cluster `absvars', suite(avar)) nocons
	reghdfe `lhs' `rhs', a(`absvars') vce(cluster `absvars', bw(2)) nocons
	TrimMatrix `K'

	* 3. Compare
	storedresults compare benchmark e(), tol(1e-5) include( ///
		scalar: N rss F df_r ///
		matrix: trim_b trim_V ///
		macros: wexp wtype )

	storedresults compare benchmark2 e(), tol(1e-5) include(scalar: mss r2 r2_a)

	storedresults drop benchmark benchmark2
	
	* NOTE: There seems to be a loss of precision somewhere in the code or in the benchmark
	* Is it in avar/reghdfe? In ivreg2?

cd "D:/Github/reghdfe/test"
exit
