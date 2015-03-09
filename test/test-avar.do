cd "D:/Github/reghdfe" // /source
cscript "reghdfe with advanced robustness options from avar" adofile reghdfe

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
	version 13
	set seed 3453 // this seed causes a collinearity problem
	**expand 100
	gen u = uniform()
	bys turn (u): gen t1 = _n
	drop u
	version 13
	set seed 34535 // this seed DOESNT causes a collinearity problem
	gen u = uniform()
	bys turn (u): gen t2 = _n
	drop u

forv i=1/2 {
	tsset turn t`i'
	cap drop ABS_*
	do "D:/Github/reghdfe/test/test-avar-inner" `i'
}
cd "D:/Github/reghdfe/test"
exit
