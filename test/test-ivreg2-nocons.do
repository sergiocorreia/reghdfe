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

	
* [TEST] One core (what was the point of this?)
	
	areg price weight length, absorb(turn)
	TrimMatrix 2
	local areg_df_a = e(df_a)
	storedresults save areg e()

	reghdfe price weight length, a(turn) keepsingletons
	TrimMatrix 2
	storedresults compare areg e(), tol(1e-10) include(scalar: df_r F N r2 matrix: trim_b trim_V macros: wexp wtype)
	assert `areg_df_a'==e(df_a)-1

	storedresults drop areg

* [TEST] IV

	reghdfe price weight (length=disp), a(turn) ivsuite(ivreg2) keepsingletons
	TrimMatrix 2
	storedresults save benchmark e()

	reghdfe price weight (length=disp), a(turn) ivsuite(ivregress) keepsingletons
	TrimMatrix 2
	storedresults compare benchmark e(), tol(1e-10) include(scalar: df_r df_a F N r2 matrix: trim_b trim_V macros: wexp wtype)

	storedresults drop benchmark
	
* Misc
	reghdfe price weight (length=disp), a(turn) keepsingletons
	reghdfe price weight (length=disp), a(turn) vce(cluster rep) ivsuite(ivregress) keepsingletons
	*reghdfe price weight (length=disp), a(turn) vce(cluster rep) ivsuite(ivregress) cores(2)
	
	reghdfe price weight (length=disp), a(turn) vce(cluster rep) keepsingletons
	
cd "C:/Git/reghdfe/test"
exit
