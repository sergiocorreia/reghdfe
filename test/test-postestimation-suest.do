cd "D:/Github/reghdfe/source"
cscript "reghdfe postestimation: suest" adofile reghdfe

* Setup
	discard
	clear all
	set more off
	* cls

	local included_e scalar: N rmse rss r2 r2_a df_r /// tss df_a df_m F 
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
	local rhs1 weight length      disp trunk
	local rhs2 weight        gear disp trunk
	local absvars turn
	fvunab tmp : `rhs1'
	local K1 : list sizeof tmp
	fvunab tmp : `rhs2'
	local K2 : list sizeof tmp

	local testvars disp trunk

* [TEST]

	* 1. Run benchmark
	qui tab `absvars', gen(ABS_)

	reg `lhs' `rhs1' ABS_*, hascons
	TrimMatrix `K1'
	storedresults save bench_e1 e()
	estimates store bench_e1

	reg `lhs' `rhs2' ABS_*, hascons
	TrimMatrix `K2'
	storedresults save bench_e2 e()
	estimates store bench_e2

	suest bench_e1 bench_e2
	test ([bench_e1_mean]displacement==[bench_e2_mean]displacement) ([bench_e1_mean]trunk==[bench_e2_mean]trunk)

	* 2. Run reghdfe and compare
	reghdfe `lhs' `rhs1', absorb(FE1=`absvars') nocons notes(cmd=regress) // vce(, suite(default))
	TrimMatrix `K1'
	storedresults compare bench_e1 e(), tol(1e-10) include(`included_e')
	estimates store m1

	reghdfe `lhs' `rhs2', absorb(FE2=`absvars') nocons notes(cmd=regress) // vce(, suite(default))
	TrimMatrix `K2'
	storedresults compare bench_e2 e(), tol(1e-10) include(`included_e')
	estimates store m2

	storedresults drop bench_e1 bench_e2
set trace off
	suest m1 m2
	test ([m1_mean]displacement==[m2_mean]displacement) ([m1_mean]trunk==[m2_mean]trunk)
	

	estimates clear

cd "D:/Github/reghdfe/test"
exit
