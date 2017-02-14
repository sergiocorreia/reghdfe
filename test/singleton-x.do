// See:
// http://www.stata.com/statalist/archive/2003-06/msg00646.html

noi cscript "reghdfe: singleton dummy as regressor" adofile reghdfe

* Dataset
	set obs 10
	gen y = rnormal()
	gen x1 = _n==1
	gen x2 = _n==5
	gen i = _n
	gen c = 1
	
	local included_e ///
		scalar: N tss rss r2 r2_a ll ll_0 /// F mss df_r rmse df_m
		matrix: trim_b trim_V ///
		macros: wexp wtype

* areg reports F missing if there are singleton regressors
* areg never reports mss

* areg reports df_r = N - rank + 1, 
* BUT it seems it uses the rank of inv(XX) when computing the small sample adjustment
* and it uses the rank of inv(V) when computing the DoF of the FStat
* The will disagree with clustered SEs when there are singleton regressors

* as a consequence, rmse and df_m are also different

* [TEST] Unadjusted

	local lhs y
	local rhs x*
	local absvars c
	local clustervars i
	fvunab tmp : `rhs'
	local K : list sizeof tmp
	
	* 1. Run benchmark
	areg `lhs' `rhs', absorb(`absvars') vce(cluster `clustervars')
	trim_cons
	local bench_df_a = e(df_a)
	local bench_F = e(F)
	storedresults save benchmark e()
	
	* 2. Run reghdfe
	reghdfe `lhs' `rhs', absorb(`absvars') vce(cluster `clustervars') // verbose(-1) keepsingletons
	notrim
	storedresults compare benchmark e(), tol(1e-12) include(`included_e')
	assert `bench_df_a'==e(df_a)-1
	assert `bench_F'==. & e(F)==.

	* Done!

storedresults drop benchmark
exit
