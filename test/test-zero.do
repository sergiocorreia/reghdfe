noi cscript "reghdfe with rhs equal to zero" adofile reghdfe

* Setup
	discard
	clear all
	set more off
	* cls

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
	drop if missing(rep)
	gen byte zero = 0
	
	local included_e scalar: N rmse tss rss r2 r2_a F df_r df_m F_absorb ll ll_0 ///
		matrix: trim_b trim_V ///
		macros: wexp wtype

* [TEST]
	keep price zero turn trunk
	reghdfe price zero, a(turn trunk) v(5) maxiter(20) pool(5)

cd "C:/Git/reghdfe/test"
exit
