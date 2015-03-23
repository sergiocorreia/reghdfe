noi cscript "reghdfe with ivreg2 should use nocons" adofile reghdfe

* Setup
	discard
	clear all
	set more off

* Convenience
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

	
* [TEST] One core
set trace off
	
	areg price weight length, absorb(turn)
	TrimMatrix 2
	storedresults save areg e()

	reghdfe price weight length, a(turn)
	TrimMatrix 2
	storedresults compare areg e(), tol(1e-10) include(scalar: df_r df_a F N r2 matrix: trim_b trim_V macros: wexp wtype)

	reghdfe price weight length, a(turn) nocons
	TrimMatrix 2
	storedresults compare areg e(), tol(1e-10) include(scalar: df_r df_a F N r2 matrix: trim_b trim_V macros: wexp wtype)
	storedresults drop areg

* [TEST] IV
set trace off

	reghdfe price weight (length=disp), a(turn) ivsuite(ivreg2)
	TrimMatrix 2
	storedresults save benchmark e()

	reghdfe price weight (length=disp), a(turn) ivsuite(ivregress)
	TrimMatrix 2
	storedresults compare benchmark e(), tol(1e-10) include(scalar: df_r df_a F N r2 matrix: trim_b trim_V macros: wexp wtype)

	reghdfe price weight (length=disp), a(turn) ivsuite(ivregress) nocons
	TrimMatrix 2
	storedresults compare benchmark e(), tol(1e-10) include(scalar: df_r df_a /* F */ N r2 matrix: trim_b trim_V macros: wexp wtype)
	* Sadly the FStat is missing since -ivregress- believes there is no constant

	storedresults drop benchmark
	
	
* Misc
	reghdfe price weight (length=disp), a(turn)
	reghdfe price weight (length=disp), a(turn) vce(cluster rep) ivsuite(ivregress)
	*reghdfe price weight (length=disp), a(turn) vce(cluster rep) ivsuite(ivregress) cores(2)
	
	reghdfe price weight (length=disp), a(turn) vce(cluster rep)
	
cd "D:/Github/reghdfe/test"
exit
