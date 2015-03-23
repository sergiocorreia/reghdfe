cscript "reghdfe with singletons" adofile reghdfe

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

	local excluded 	macro: cmdline corr1 corr2 corr3 corr4 prettynames ///
					scalar: N df_a K1 K2 M2 r2 r2_a tss mss _cons mobility F_absorb ///
					matrix: b V
	* Note that R2-within shouldn't change; also we don't compare -b- instead we compare trim_b

* [TEST]
	local lhs price
	local rhs weight length
	local absvars turn trunk foreign
	fvunab tmp : `rhs'
	local K : list sizeof tmp
	
* Unadjusted
	reghdfe `lhs' `rhs', absorb(`absvars') nocons
	TrimMatrix `K'
	storedresults save benchmark e()

	reghdfe `lhs' `rhs', absorb(`absvars') dropsingle
	TrimMatrix `K'
	storedresults compare benchmark e(), tol(1e-10) exclude(`excluded')
	storedresults drop benchmark

* Unadj with constant
	reghdfe `lhs' `rhs', absorb(`absvars')
	TrimMatrix `K'
	storedresults save benchmark e()

	reghdfe `lhs' `rhs', absorb(`absvars') dropsingle
	TrimMatrix `K'
	storedresults compare benchmark e(), tol(1e-10) exclude(`excluded')
	storedresults drop benchmark

* Unadjusted w/cont
	* This DOESNT MATCH (!!) because when dropping singletons we get a better DoF (closer to the truth). All in all, this makes the variances *lower*
	local absvars turn trunk##c.gear foreign
	
	reghdfe `lhs' `rhs', absorb(`absvars') nocons
	tempname trim_diag_V
	TrimMatrix `K'
	matrix `trim_diag_V' = vecdiag(e(trim_V))
	storedresults save benchmark e()

	reghdfe `lhs' `rhs', absorb(`absvars') dropsingle
	TrimMatrix `K'
	storedresults compare benchmark e(), tol(1e-10) include(matrix: trim_b) // exclude(`excluded')
	tempname delta zeros
	matrix `delta' = `trim_diag_V' - vecdiag(e(trim_V))
	mata: assert(all(st_matrix("`delta'") :> 0))
	storedresults drop benchmark

* Many more absvars
	local absvars turn trunk rep foreign

	reghdfe `lhs' `rhs', absorb(`absvars') nocons
	TrimMatrix `K'
	storedresults save benchmark e()

	reghdfe `lhs' `rhs', absorb(`absvars') dropsingle
	TrimMatrix `K'
	storedresults compare benchmark e(), tol(1e-10) exclude(`excluded')
	storedresults drop benchmark

* Robust  =====> It seems there IS a difference here..
* I kinda expected it for cluster-robust, not just for robust itself.

	**local absvars turn trunk // foreign
	**
	**reghdfe `lhs' `rhs', absorb(`absvars') nocons vce(robust)
	**TrimMatrix `K'
	**storedresults save benchmark e()
	**
	**reghdfe `lhs' `rhs', absorb(`absvars') dropsingle vce(robust)
	**TrimMatrix `K'
	**storedresults compare benchmark e(), tol(1e-10) exclude(`excluded')
	**storedresults drop benchmark

cd "D:/Github/reghdfe/test"
exit
