noi cscript "reghdfe: complex factor variables that fvstrip handles" adofile reghdfe

cap pr drop TrimBoth
pr TrimBoth, eclass
	matrix trim_b = e(b)
	matrix trim_V = e(V)
	loc size = colsof(trim_b) - 1

	mata: idx = 2, 3, 5
	mata: b = st_matrix("trim_b")[idx]
	mata: st_matrix("trim_b", b)
	mata: V = st_matrix("trim_V")[idx, idx]
	mata: st_matrix("trim_V", V)

	ereturn matrix trim_b = trim_b
	ereturn matrix trim_V = trim_V
end

cap pr drop TrimBoth2
pr TrimBoth2, eclass
	matrix trim_b = e(b)
	matrix trim_V = e(V)
	loc size = colsof(trim_b) - 1

	// mata: idx = 2..13 , 15..25 , 27 , 30
	mata: idx = 2..13 , 15..25 , 27 , 29, 30
	mata: b = st_matrix("trim_b")[idx]
	mata: st_matrix("trim_b", b)
	mata: V = st_matrix("trim_V")[idx, idx]
	mata: st_matrix("trim_V", V)

	ereturn matrix trim_b = trim_b
	ereturn matrix trim_V = trim_V
end




local included_e ///
	scalar: N rmse tss rss mss r2 r2_a F df_r df_m ll ll_0 ///
	matrix: trim_b trim_V ///
	macros: wexp wtype



* [TEST] https://github.com/sergiocorreia/reghdfe/issues/91
	clear
	set obs 8
	gen byte x = _n <= 4
	gen byte y = mod(_n+1, 4) < 2
	gen byte z = mod(_n+1, 2)
	expand 2
	gen w = runiform()
	gen byte c = 1
	li, sepby(x)

	gen byte a1 = 1.x
	gen byte a2 = 1.y
	gen byte a3 = 1.z
	gen byte a4 = 1.x#1.y
	gen byte a5 = 1.x#1.z

	* 1) Benchmark
	reghdfe w a*, keepsing a(c)
	storedresults save benchmark e()
	
	* 2) Reghdfe
	reghdfe w i.x##i.(y z), keepsing a(c) v(1)
	storedresults compare benchmark e(), exclude(macro: indepvars cmdline)

	* 3) Cleanup
	storedresults drop benchmark


* [TEST] https://github.com/sergiocorreia/reghdfe/issues/86

	sysuse auto, clear
	gen byte turn43 = (turn==43)

	* 1) Benchmark
	reghdfe price weight, absorb(rep78) vce(cl foreign)
	storedresults save benchmark e()
	
	* 2) Reghdfe
	reghdfe price weight, absorb(rep78) vce(cl i.foreign)
	storedresults compare benchmark e(), exclude(macro: depvar cmdline)

	* 3) Cleanup
	storedresults drop benchmark

****

	sysuse auto, clear
	gen byte turn43 = (turn==43)

	* 1) Benchmark
	reghdfe price weight, absorb(rep78) vce(cl foreign#turn)
	storedresults save benchmark e()
	
	* 2) Reghdfe
	reghdfe price weight, absorb(rep78) vce(cl i.foreign#i.turn)
	storedresults compare benchmark e(), exclude(macro: depvar cmdline)

	* 3) Cleanup
	storedresults drop benchmark


* [TEST] https://github.com/sergiocorreia/reghdfe/issues/88

	* Setup
	sysuse auto, clear
	bys turn: gen t = _n
	xtset turn t

	* 1) Benchmark
	areg price L.(weight gear), a(turn)
	trim_cons
	storedresults save benchmark e()
	local bench_df_a = e(df_a)

	* 2) Compare w/reghdfe
	qui reghdfe price L.weight, noa keepsing v(1)
	qui reghdfe price L.(weight), noa keepsing v(1)
	reghdfe price L.(weight gear), a(turn) keepsing v(1)
	notrim
	matrix list e(trim_V)
	storedresults compare benchmark e(), tol(1e-10) include(`included_e')
	assert `bench_df_a'==e(df_a)-1

	* 3) Cleanup
	storedresults drop benchmark


* [TEST] https://github.com/sergiocorreia/reghdfe/issues/28

	sysuse auto, clear
	gen byte turn43 = (turn==43)

	* 1) Benchmark
	reghdfe turn43 price weight, a(trunk)
	storedresults save benchmark e()
	
	* 2) Reghdfe
	reghdfe 43.turn price weight, a(trunk)
	storedresults compare benchmark e(), exclude(macro: depvar cmdline)

	* 3) Cleanup
	storedresults drop benchmark


* [TEST] https://github.com/sergiocorreia/reghdfe/issues/22

	sysuse auto, clear
	gen x = foreign

	* 1) Benchmark
	areg price i.x##c.gear, a(turn)
	matrix list e(V)
	TrimBoth
	matrix list e(trim_V)
	local bench_df_a = e(df_a)
	storedresults save benchmark e()
	
	* 2) Reghdfe
	reghdfe price i.x##c.gear, a(turn) keepsing v(-1)
	notrim
	matrix list e(trim_V)
	storedresults compare benchmark e(), tol(1e-10) include(`included_e')
	assert `bench_df_a'==e(df_a)-1

	* 3) Cleanup
	storedresults drop benchmark


* [TEST] https://github.com/sergiocorreia/reghdfe/issues/35

	sysuse auto, clear
	fvset base 1 foreign

	* 1) Benchmark
	areg price weight, a(foreign)
	trim_cons
	local bench_df_a = e(df_a)
	storedresults save benchmark e()
	
	* 2) Reghdfe
	reghdfe price weight, a(foreign) keepsing v(-1)
	notrim
	storedresults compare benchmark e(), tol(1e-10) include(`included_e')
	assert `bench_df_a'==e(df_a)-1

	* 3) Cleanup
	storedresults drop benchmark


* [TEST] https://github.com/sergiocorreia/reghdfe/issues/43
	
	sysuse nlsw88.dta, clear

	* 1) Benchmark
	areg wage i.industry##c.age i.union##c.age, a(race)
	matrix list e(b)
	TrimBoth2
	matrix list e(trim_b)
	local bench_df_a = e(df_a)
	storedresults save benchmark e()
	
	* 2) Reghdfe
	reghdfe wage i.industry##c.age i.union##c.age, a(race) keepsing  v(-1)
	notrim
	storedresults compare benchmark e(), tol(1e-9) include(`included_e')
	assert `bench_df_a'==e(df_a)-1

	* 3) Cleanup
	storedresults drop benchmark

* [TEST] Bug in ms_fvstrip reported by Gabriel C-R
	
	sysuse auto
	bys turn: gen t = _n
	xtset turn t
	gen F = 1
	gen L = 1
	reghdfe F.price L.weight, noa



cd "C:/Git/reghdfe/test"
exit
