noi cscript "reghdfe with complex factor variables that fvstrip handles" adofile reghdfe

* Setup
	discard
	clear all
	set more off

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
	
	local included_e scalar: N rmse tss rss r2 r2_a F df_r df_m F_absorb ll ll_0 ///
		matrix: trim_b trim_V ///
		macros: wexp wtype

* [TEST] https://github.com/sergiocorreia/reghdfe/issues/28

	sysuse auto, clear
	gen byte turn43 = (turn==43)

	* 1) Benchmark
	reghdfe turn43 price weight, a(trunk)
	TrimMatrix 2
	local bench_df_a = e(df_a)
	storedresults save benchmark e()
	
	* 2) Reghdfe
	reghdfe 43.turn price weight, a(trunk)
	TrimMatrix 2
	storedresults compare benchmark e(), tol(1e-10) include(`included_e')
	assert `bench_df_a'==e(df_a)-1

	* 3) Cleanup
	storedresults drop benchmark





* [TEST] https://github.com/sergiocorreia/reghdfe/issues/22

	sysuse auto, clear
	gen x = foreign

	* 1) Benchmark
	areg price i.x##c.gear, a(turn)
	TrimMatrix 2
	local bench_df_a = e(df_a)
	storedresults save benchmark e()
	
	* 2) Reghdfe
	reghdfe price i.x##c.gear, a(turn) keepsing
	TrimMatrix 2
	storedresults compare benchmark e(), tol(1e-10) include(`included_e')
	assert `bench_df_a'==e(df_a)-1

	* 3) Cleanup
	storedresults drop benchmark




* [TEST] https://github.com/sergiocorreia/reghdfe/issues/35

	sysuse auto, clear
	fvset base 1 foreign

	* 1) Benchmark
	areg price weight, a(foreign)
	TrimMatrix 2
	local bench_df_a = e(df_a)
	storedresults save benchmark e()
	
	* 2) Reghdfe
	reghdfe price weight, a(foreign) keepsing
	TrimMatrix 2
	storedresults compare benchmark e(), tol(1e-10) include(`included_e')
	assert `bench_df_a'==e(df_a)-1

	* 3) Cleanup
	storedresults drop benchmark

* [TEST] https://github.com/sergiocorreia/reghdfe/issues/43
	
	sysuse nlsw88.dta, clear

	* 1) Benchmark
	areg wage i.industry##c.age i.union##c.age, a(race)
	TrimMatrix 31
	local bench_df_a = e(df_a)
	storedresults save benchmark e()
	
	* 2) Reghdfe
	reghdfe wage i.industry##c.age i.union##c.age, a(race) keepsing
	TrimMatrix 31
	storedresults compare benchmark e(), tol(1e-10) include(`included_e')
	assert `bench_df_a'==e(df_a)-1

	* 3) Cleanup
	storedresults drop benchmark



cd "C:/Git/reghdfe/test"
exit
