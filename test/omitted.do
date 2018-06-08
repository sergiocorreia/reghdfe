noi cscript "reghdfe: ensure we omit the correct variables" adofile reghdfe

* Setup
	local included_e ///
		scalar: N rmse tss rss r2 r2_a df_r df_m ll ll_0 /// mss F
		matrix: trim_b trim_V ///
		macros: wexp wtype

* [TEST] Simplest case
	sysuse auto, clear
	gen z = 0
	reghdfe price z, noab
	matrix b = e(b)
	assert b[1,1]==0
	matrix list b
	_ms_omit_info b
	assert r(k_omit)==1


* [TEST] Simple case
* https://github.com/sergiocorreia/reghdfe/issues/93

	sysuse auto, clear
	gen lowrep = rep78<4
	loc lhs price
	loc rhs lowrep
	loc absvars rep78

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

* [TEST] https://github.com/sergiocorreia/reghdfe/issues/95
	clear
	set obs 100
	set seed 123
	gen double y = runiform() * 1e8
	gen double z = runiform()
	gen byte c = 1

	loc lhs z
	loc rhs y
	loc absvars c

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


* [TEST] https://github.com/sergiocorreia/reghdfe/issues/94
	clear
input w x y unity
w          x          y      unity
1    -1.10565   -1.223344       1
1   -1.109139   -1.127615       1
0   -.8497003   -1.207231       1
1   -.2734231   -1.294683       1
0    .0063879   -.8343916       1
1    .7669702   -.8245103       1 
end

	gen z = 0
	reghdfe y z, noab

	reghdfe y i.w##i.unity x, a(unity) tol(1e-12)
	matrix list e(b)
	matrix list e(V)

	areg y i.w##i.unity x, a(unity)
	matrix list e(b)
	matrix list e(V)


	loc lhs y
	loc rhs i.w##i.unity x
	loc absvars unity

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

	* Done!

storedresults drop benchmark
clear matrix
exit
