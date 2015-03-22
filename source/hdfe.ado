*! hdfe 1.0 12mar2015
*! Sergio Correia (sergio.correia@duke.edu)
// -------------------------------------------------------------------------------------------------
// Partial-out a list of variables with respect to any number of fixed effects
// -------------------------------------------------------------------------------------------------
/// Does NOT accept factor/time series, and variables MUST BE FULLY spelled out (i.e. you need to use unab beforehand!)

cap pr drop hdfe
program define hdfe, rclass
	syntax varlist [fweight aweight pweight/] , Absorb(string) [CLUSTERvars(string) Verbose(integer 0) TOLerance(real 1e-7) MAXITerations(integer 10000)] [GENerate(name)]
	if ("`weight'"!="") {
		local weightvar `exp'
		conf var `weightvar' // just allow simple weights
		local weighttype `weight'
	}
	
* Assert that base program exists
	qui which reghdfe_absorb

* Preserve if asked to
	if ("`generate'"!="") {
		tempvar uid
		gen double `uid' = _n
		preserve
	}
	
* Clear previous errors
	cap reghdfe_absorb, step(stop)

* Time/panel variables
	cap conf var `_dta[_TStvar]'
	if (!_rc) local timevar `_dta[_TStvar]'
	cap conf var `_dta[_TSpanel]'
	if (!_rc) local panelvar `_dta[_TSpanel]'

* Set Verbosity
	mata: VERBOSE = `verbose' // Pick a number between 0 (quiet) and 4 (lots of debugging info)

* Parse: absorb, clusters, and weights
	reghdfe_absorb, step(start) absorb(`absorb') clustervars(`clustervars') weight(`weighttype') weightvar(`weightvar')
	local absorb_keepvars = r(keepvars)
	local N_hdfe = r(N_hdfe)
	
* Keep relevant observations
	marksample touse, novar
	markout `touse' `varlist' `absorb_keepvars'
	qui keep if `touse'
	
* Keep relevant variables
	keep `varlist' `clustervars' `weightvar' `panelvar' `timevar' `absorb_keepvars' `uid'
	
* Construct Mata objects and auxiliary variables
	reghdfe_absorb, step(precompute) keep(`varlist' `clustervars' `weightvar' `panelvar' `timevar' `uid') tsvars(`panelvar' `timevar')

* Compute e(df_a)
	reghdfe_absorb, step(estimatedof) dofadjustments(pairwise clusters continuous)
	* return list // what matters is r(kk) which will be e(df_a)
	local kk = r(kk)
	
* Demean variables wrt to the fixed effects
	reghdfe_absorb, step(demean) varlist(`varlist') tol(`tolerance') maxiterations(`maxiterations') // Other maximize/parallel options
	
	return scalar df_a = `kk'
	return scalar N_hdfe = `N_hdfe'
	forv g=1/`N_hdfe' {
		reghdfe_absorb, fe2local(`g') // copies Mata structure into locals
		* Will inject the following with c_local:
		* ivars cvars target varname varlabel is_interaction is_cont_interaction is_bivariate is_mock levels
		return local hdfe`g' = "`varlabel'"
		return scalar df_a`g' = `levels'
	}
	
* Clean up Mata objects
	reghdfe_absorb, step(stop)

	if ("`generate'"!="") {
		foreach var of local varlist {
			rename `var' `generate'`var'
		}

		tempfile output
		sort `uid'
		qui save "`output'"
		restore
		SafeMerge, uid(`uid') file("`output'")
	}
end

* [SafeMerge: ADAPTED FROM THE ONE IN ESTIMATE.ADO]
* The idea of this program is to keep the sort order when doing the merges
cap pr drop SafeMerge
program define SafeMerge, eclass sortpreserve
syntax, uid(varname numeric) file(string)
	* Merging gives us e(sample) and the FEs / AvgEs
	merge 1:1 `uid' using "`file'", assert(master match) nolabel nonotes noreport nogen
end
