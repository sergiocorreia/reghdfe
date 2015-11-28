cscript "reghdfe with iv-stages" adofile reghdfe

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
	local endogvar gear
	local iv head displace
	local absvars turn
	local clustervar trunk


	* Custom adjustment: in this simple case we can compare _cons
	local K = `K' + 1
	drop if missing(rep)

	* 1. Run benchmarks
	qui tab `absvars', gen(ABS_)
	
	* IV
	local cmd ivreg2 `lhs' `rhs' (`endogvar'=`iv') ABS_* , small partial(ABS_*) cluster(`clustervar') 
	di as input "`cmd'"
	`cmd'
	fvunab tmp : `rhs' `endogvar'
	local K_iv : list sizeof tmp
	TrimMatrix `K_iv'
	estadd local stage = "iv"
	storedresults save b_iv e()

	* OLS
	areg `lhs' `endogvar' `rhs' , abs(`absvars') cluster(`clustervar')
	fvunab tmp : `rhs'
	local K_ols : list sizeof tmp
	TrimMatrix `K_ols'
	estadd local stage = "ols"
	storedresults save b_ols e()

	* Reduced
	areg `lhs' `iv' `rhs', abs(`absvars') cluster(`clustervar')
	fvunab tmp : `rhs' `iv'
	local K_reduced : list sizeof tmp
	TrimMatrix `K_reduced'
	estadd local stage = "reduced"
	storedresults save b_reduced e()

	* Acid
	areg `lhs' `endogvar' `iv' `rhs', abs(`absvars') cluster(`clustervar')
	fvunab tmp : `rhs' `endogvar' `iv'
	local K_acid : list sizeof tmp
	TrimMatrix `K_acid'
	estadd local stage = "acid"
	storedresults save b_acid e()

	* First stage
	areg `endogvar' `iv' `rhs', abs(`absvars') cluster(`clustervar')
	fvunab tmp : `rhs' `iv'
	local K_first1 : list sizeof tmp
	TrimMatrix `K_first1'
	estadd local stage = "first1"
	storedresults save b_first1 e()

	* 2. Run reghdfe
	reghdfe `lhs' `rhs' (`endogvar'=`iv'), absorb(`absvars') ///
		stages(ols first acid reduced) ///
		v(3) vce(cluster `clustervar')  keepsingletons // will be saved as reghdfe_`stage'
		
	estimates store reghdfe_iv

	* Compare
	foreach stage in iv ols first1 acid reduced  {
		estimates restore reghdfe_`stage'
		TrimMatrix `K_`stage''
		storedresults compare b_`stage' e(), tol(1e-6) include( ///
			scalar: N rss df_r `extra' /// rmse
			matrix: trim_b trim_V ///
			macros: wexp wtype depvar)
		local extra r2 // Activate after comparing IV
		
	}
	
	storedresults drop b_iv b_ols b_acid b_reduced b_first1
	estimates drop reghdfe_* // accepts wildchars? or just -estimates clear-?

	* Extra
	reghdfe `lhs' `rhs' (`endogvar'=`iv'), absorb(`absvars') stages(ols first acid reduced) ffirst
	reghdfe `lhs' `rhs' (`endogvar'=`iv'), absorb(`absvars') stages(ols first acid reduced, nosave)

cd "C:/Git/reghdfe/test"
exit
