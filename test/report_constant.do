noi cscript "reghdfe: test that the constant option does not change anything" adofile reghdfe

* Dataset
	sysuse auto
	bys turn: gen t = _n
	tsset turn t
	drop if missing(rep)
	
	local lhs price
	local rhs weight length
	local absvars turn

	local exclude ///
		scalar: report_constant ///
		matrix: b V ///
		macros: indepvars cmdline


* [TEST] Unadjusted
	reghdfe `lhs' `rhs', absorb(`absvars') keepsingletons verbose(-1) constant
	trim_cons
	storedresults save benchmark e()
	reghdfe `lhs' `rhs', absorb(`absvars') keepsingletons verbose(-1) noconstant
	notrim
	storedresults compare benchmark e(), tol(1e-12) exclude(`exclude')

* [TEST] Robust
	reghdfe `lhs' `rhs', vce(robust) absorb(`absvars') keepsingletons verbose(-1) constant
	trim_cons
	storedresults save benchmark e()
	reghdfe `lhs' `rhs', vce(robust) absorb(`absvars') keepsingletons verbose(-1) noconstant
	notrim
	storedresults compare benchmark e(), tol(1e-12) exclude(`exclude')

* [TEST] Cluster
	reghdfe `lhs' `rhs', cluster(trunk) absorb(`absvars') keepsingletons verbose(-1) constant
	trim_cons
	storedresults save benchmark e()
	reghdfe `lhs' `rhs', cluster(trunk) absorb(`absvars') keepsingletons verbose(-1) noconstant
	notrim
	storedresults compare benchmark e(), tol(1e-12) exclude(`exclude')

* [TEST] Cluster 2
	reghdfe `lhs' `rhs', cluster(turn) absorb(`absvars') keepsingletons verbose(-1) constant
	trim_cons
	storedresults save benchmark e()
	reghdfe `lhs' `rhs', cluster(turn) absorb(`absvars') keepsingletons verbose(-1) noconstant
	notrim
	storedresults compare benchmark e(), tol(1e-12) exclude(`exclude')

* [TEST] Drop singletons

	reghdfe `lhs' `rhs', absorb(`absvars') verbose(-1) constant
	trim_cons
	storedresults save benchmark e()
	reghdfe `lhs' `rhs', absorb(`absvars') verbose(-1) noconstant
	notrim
	storedresults compare benchmark e(), tol(1e-12) exclude(`exclude')
	exit




	fvunab tmp : `rhs'
	local K : list sizeof tmp
	
	* 1. Run benchmark
	areg `lhs' `rhs', absorb(`absvars')
	matrix list e(V)
	trim_cons
	local bench_df_a = e(df_a)
	
	* 2. Run reghdfe
	notrim
	assert `bench_df_a'==e(df_a)-1

	* 3. Run reghdfe vce(ols)
	reghdfe `lhs' `rhs', absorb(`absvars') keepsingletons verbose(-1) vce(ols)
	notrim
	storedresults compare benchmark e(), tol(1e-12) include(`included_e')
	assert `bench_df_a'==e(df_a)-1

	* 4. Run reghdfe vce(unadjusted)
	reghdfe `lhs' `rhs', absorb(`absvars') keepsingletons verbose(-1) vce(unadjusted)
	notrim
	storedresults compare benchmark e(), tol(1e-12) include(`included_e')
	assert `bench_df_a'==e(df_a)-1

	* Done!

storedresults drop benchmark
exit
