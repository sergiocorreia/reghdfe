cap pr drop partialfe
program define partialfe, rclass
	syntax varlist [fweight aweight pweight/] , Absorb(string) [CLUSTERvars(string) Verbose(integer 0) TOLerance(real 1e-7) MAXITerations(integer 10000)]
	if ("`weight'"!="") {
		local weightvar `exp'
		conf var `weightvar' // just allow simple weights
		local weighttype `weight'
	}
	
* Assert that program exists
	qui which reghdfe_absorb	
	
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
	keep if `touse'
	
* Keep relevant variables
	keep `varlist' `clustervars' `weightvar' `panelvar' `timevar' `absorb_keepvars'
	
* Construct Mata objects and auxiliary variables
	reghdfe_absorb, step(precompute) keep(`varlist' `clustervars' `weightvar' `panelvar' `timevar') tsvars(`panelvar' `timevar')

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
end
