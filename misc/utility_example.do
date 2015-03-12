**** Preamble ****

which reghdfe_absorb
sysuse auto, clear
gen n = 1
cap cls

* Relevant variables
	local absvars 		trunk rep
	local clustervars	turn
	local depvar 		price
	local indepvars 	weight length
	local endogvars		gear
	local instruments	head displace
	
* Weights, if needed
	local weighttype 	// fweight
	local weightvar 	// n
	
* Time/panel variables
	cap conf var `_dta[_TStvar]'
	if (!_rc) local timevar `_dta[_TStvar]'
	cap conf var `_dta[_TSpanel]'
	if (!_rc) local panelvar `_dta[_TSpanel]'

* Clear previous errors
	cap reghdfe_absorb, step(stop)
	
* Benchmark regressions
	qui tab trunk, gen(ABS1_)
	qui tab rep, gen(ABS2_)
	if ("`weightvar'"!="") local weightexp "[`weighttype'=`weightvar']"
	
	* Benchmark 1: IVREG WITH PARTIAL
	ivreg2 `depvar' `indepvars' ABS1_* ABS2_* (`endogvars'=`instruments') `weightexp', cluster(`clustervars') small nocons partial(ABS1_* ABS2_*)
	drop ABS1_* ABS2_*
	
	* Benchmark 2: REGHDFE
	reghdfe `depvar' `indepvars' (`endogvars'=`instruments') `weightexp', vce(cluster `clustervars') absorb(`absvars')


**** Start ****

* 1) Set Verbosity
	mata: VERBOSE = 0 // Pick a number between 0 (quiet) and 4 (lots of debugging info)

* 2) Parse: absorb, clusters, and weights
	reghdfe_absorb, step(start) absorb(`absvars') clustervars(`clustervars') weight(`weighttype') weightvar(`weightvar')
	local absorb_keepvars = r(keepvars)
	local N_hdfe = r(N_hdfe)
	
* x) Interim steps:
	preserve

	* Expand factor variables, if needed
	* fvrevar ...
	local main_keepvars `absvars' `clustervars' `depvar' `indepvars' `endogvars' `instruments' `weightvar' `panelvar' `timevar' // add created factor/time variables to this list
	
	* Keep relevant observations
	marksample touse, novar
	markout `touse' `main_keepvars' `absorb_keepvars'
	keep if `touse'
	
	* Keep relevant variables
	keep `main_keepvars' `absorb_keepvars'
	
* 3) Construct Mata objects and auxiliary variables
	reghdfe_absorb, step(precompute) keep(`main_keepvars') depvar("`depvar'") tsvars(`panelvar' `timevar')

* 4) Compute e(df_a)
	reghdfe_absorb, step(estimatedof) dofadjustments(pairwise clusters continuous)
	return list // what matters is r(kk) which will be e(df_a)
	local kk = r(kk)
	
* 5) Demean variables wrt to the fixed effects
	reghdfe_absorb, step(demean) varlist(`depvar' `indepvars' `endogvars' `instruments') tol(1e-10) maxiter(1000) // Other maximize/parallel options
	
* 6) Now we can run the regression
	if ("`weightvar'"!="") local weightexp "[`weighttype'=`weightvar']"
	ivreg2 `depvar' `indepvars' (`endogvars'=`instruments') `weightexp', cluster(`clustervars') small sdofminus(`=`kk'+1') nocons
	
* 7) Report FEs:
	forv g=1/`N_hdfe' {
		reghdfe_absorb, fe2local(`g') // copies Mata structure into locals
		* Will inject the following with c_local:
		* ivars cvars target varname varlabel is_interaction is_cont_interaction is_bivariate is_mock levels
		di as text" - Fixed effect `g': `varlabel'; `levels' categories"
	}

* 8) Clean up Mata objects
	reghdfe_absorb, step(stop)
	
	restore	
	
exit
