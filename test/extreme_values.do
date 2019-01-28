noi cscript "reghdfe: ols with extreme regressor values" adofile reghdfe

* Dataset
use extreme_values

* Params
local included_e ///
	scalar: N rmse tss rss mss r2 df_r ///
	matrix: trim_b trim_V ///
	macros: wexp wtype

* [TEST] Unadjusted

	local lhs y
	local rhs x*
	fvunab tmp : `rhs'
	local K : list sizeof tmp
	
	* 1. Run benchmark
	areg `lhs' `rhs' i.j, absorb(i)
	trim_cons 4
	matrix list e(trim_b)
	local bench_df_a = e(df_a)
	storedresults save benchmark e()
	
	* 2. Run reghdfe
	reghdfe `lhs' `rhs', absorb(i j) keepsingletons
	trim_cons 4
	matrix list e(trim_b)
	storedresults compare benchmark e(), tol(1e-12) include(`included_e')


storedresults drop benchmark
exit
