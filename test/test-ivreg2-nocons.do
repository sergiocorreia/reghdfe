cd "D:/Github/reghdfe/source"
cscript "reghdfe with ivreg2 should use nocons" adofile reghdfe

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
	reghdfe price weight length, a(turn)
	reghdfe price weight (length=disp), a(turn)
	reghdfe price weight (length=disp), a(turn) vce(cluster rep) ivsuite(ivregress)
	*reghdfe price weight (length=disp), a(turn) vce(cluster rep) ivsuite(ivregress) cores(2)
	
	reghdfe price weight (length=disp), a(turn) vce(cluster rep)
	
	
exit
