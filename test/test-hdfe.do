cd "D:/Github/reghdfe" // /source
cscript "reghdfe with clusters" adofile reghdfe

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

* [TEST] Cluster
	local lhs price
	local rhs weight length
	local absvars turn
	local clustervar rep78
	fvunab tmp : `rhs'
	local K : list sizeof tmp

	* Custom adjustment: in this simple case we can compare _cons
	local K = `K' + 1
	drop if missing(rep)

	* 1. Run -reghdfe- as benchmark
	di as result "reghdfe `lhs' `rhs', absorb(`absvars') vce(cluster `clustervar')"
	reghdfe `lhs' `rhs', absorb(`absvars') vce(cluster `clustervar') // dof(none)
	matrix list e(V)
	TrimMatrix `K'

	* 2. Run -hdfe-
	preserve
	hdfe `lhs' `rhs', abs(`absvars') cluster(`clustervar')
	return list
	regress `lhs' `rhs' , vce(cluster `clustervar') // this will have different DoF so of course different VCE
	restore

	hdfe `lhs' `rhs', abs(`absvars') cluster(`clustervar') gen(resid_)
	reg resid_*
	
cd "D:/Github/reghdfe/test"
exit
