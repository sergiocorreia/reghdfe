noi cscript "reghdfe: bug in st_data" adofile reghdfe

cap pr drop TrimConsAndBase
pr TrimConsAndBase, eclass
args start size adjdof
	assert `start'>=1
	matrix trim_b = e(b)
	matrix trim_V = e(V)

	if ("`size'"=="") {
		loc size = colsof(trim_b) - 1
	}
	assert `size'>0

	matrix trim_b = trim_b[1, `start'..`size']
	matrix trim_V = trim_V[`start'..`size', `start'..`size']

	if ("`adjdof'"!="") {
		matrix trim_V = trim_V * `adjdof'
		ereturn scalar F = e(F) / `adjdof'
	}

	ereturn matrix trim_b = trim_b
	ereturn matrix trim_V = trim_V
end


* Dataset
	sysuse auto
	bys turn: gen t = _n
	tsset turn t
	drop if missing(rep)
	
	local included_e ///
		scalar: N rmse tss rss mss r2 r2_a F df_r df_m ll ll_0 ///
		matrix: trim_b trim_V ///
		macros: wexp wtype

* [TEST] 1.x 2.x sometimes selects less variables

	local lhs price
	local rhs 32bn.turn 43.turn 51.turn
	local absvars foreign
	fvunab tmp : `rhs'
	local K : list sizeof tmp
	
	* 1. Run benchmark
	areg `lhs' `rhs', absorb(`absvars')
	matrix list e(b)
	trim_cons
	local bench_df_a = e(df_a)
	storedresults save benchmark e()
	
	* 2. Run reghdfe
	reghdfe `lhs' `rhs', absorb(`absvars') keepsingletons verbose(-1)
	matrix list e(b)
	notrim
	storedresults compare benchmark e(), tol(1e-12) include(`included_e')
	assert `bench_df_a'==e(df_a)-1

	* Done!
	storedresults drop benchmark

* [TEST] Now with a base level

	local lhs price
	local rhs 32.turn 43.turn 51.turn
	local absvars foreign
	fvunab tmp : `rhs'
	local K : list sizeof tmp
	
	* 1. Run benchmark
	areg `lhs' 32.turn 43.turn 51.turn, absorb(`absvars')
	matrix list e(b)
	TrimConsAndBase 2
	local bench_df_a = e(df_a)
	storedresults save benchmark e()
	
	* 2. Run reghdfe
	reghdfe `lhs' `rhs', absorb(`absvars') keepsingletons verbose(-1)
	matrix list e(b)
	notrim
	storedresults compare benchmark e(), tol(1e-12) include(`included_e')
	assert `bench_df_a'==e(df_a)-1

	* Done!
	storedresults drop benchmark


* [TEST] Interactions fail

	local lhs price
	local rhs 0.foreign#31.turn 1.foreign#32.turn
	local absvars trunk
	fvunab tmp : `rhs'
	local K : list sizeof tmp
	
	* 1. Run benchmark
	//gen double X1 = (foreign==0)*(turn==31)
	//gen double X2 = (foreign==1)*(turn==32)
	areg `lhs' `rhs', absorb(`absvars')
	matrix list e(b)
	TrimConsAndBase 4
	local bench_df_a = e(df_a)
	storedresults save benchmark e()
	
	* 2. Run reghdfe
	reghdfe `lhs' `rhs', absorb(`absvars') keepsingletons verbose(-1)
	matrix list e(b)
	notrim
	storedresults compare benchmark e(), tol(1e-12) include(`included_e')
	assert `bench_df_a'==e(df_a)-1

	* Done!
	storedresults drop benchmark

* [Test] rhs collinear
	reghdfe price 0.foreign#31.turn 1.foreign#32.turn, a(turn) keepsing  v(-1)
	matrix b = e(b)
	mata: assert(st_matrix("b")==J(1, 1, 0))
	matrix V = e(V)
	mata: assert(st_matrix("V")==J(1, 1, 0))

* [Test] rhs collinear and singletons dropped
* after the singletons, both rhs vars are all zero, so they will be omitted
	reghdfe price 0.foreign#31.turn 1.foreign#32.turn, a(turn)
	matrix b = e(b)
	mata: assert(st_matrix("b")==J(0, 0, 0))
	matrix V = e(V)
	mata: assert(st_matrix("V")==J(0, 0, 0))

exit
