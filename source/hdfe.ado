*! hdfe VERSION_NUMBER
*! Sergio Correia (sergio.correia@duke.edu)
* (built from multiple source files using build.py)

include _mata/reghdfe.mata

cap pr drop hdfe
program define hdfe, rclass
	local version `=clip(`c(version)', 11.2, 13.1)' // 11.2 minimum, 13+ preferred
	qui version `version'

* Intercept version calls
	cap syntax, version
	local rc = _rc
	 if (`rc'==0) {
		Version
		exit
	}

* Intercept multiprocessor/parallel calls
	cap syntax, instance [*]
	local rc = _rc
	 if (`rc'==0) {
		ParallelInstance, `options'
		exit
	}

* Parse
	syntax varlist [if] [in] [fweight aweight pweight/] , Absorb(string) ///
		[PARTIAL(varlist numeric)] ///
		[CORES(integer 1)] ///
		[DROPSIngletons] ///
		[SAMPLE(name)] ///
		[GENerate(name)] [CLEAR] ///
		[CLUSTERVARs(string) Verbose(integer 0) TOLerance(real 1e-7) MAXITerations(integer 10000)] ///
		[noACCELerate*]

	Assert ("`generate'"!="") + ("`clear'"!="") == 1 , msg("hdfe error: you need to specify one and only one of the following options: clear generate(...)")

	if ("`accelerate'"!="") local options `options' accelerate(0) // Deal with noACCELerate


	* Check that intersection(partial,varlist) = Null
	local intersection : list varlist & partial
	Assert "`intersection'"=="", msg("variables in varlist cannot appear in partial()")

	if ("`sample'"!="") conf new var `sample'

	if ("`weight'"!="") {
		local weightvar `exp'
		conf var `weightvar' // just allow simple weights
		local weighttype `weight'
		local weightequal =
	}

* Preserve if asked to
	if ("`generate'"!="") {

		* The new var must not exist!
		foreach var of varlist `varlist' {
			conf new var `generate'`var', exact
		}

		tempvar uid
		gen double `uid' = _n
		preserve
	}

* Clear previous errors
	Stop

* From now on, we will pollute the Mata workspace, so wrap this in case of error
cap noi {

* Time/panel variables
	cap conf var `_dta[_TStvar]'
	if (!_rc) local timevar `_dta[_TStvar]'
	cap conf var `_dta[_TSpanel]'
	if (!_rc) local panelvar `_dta[_TSpanel]'

* Set Verbosity
	mata: VERBOSE = `verbose' // Pick a number between 0 (quiet) and 4 (lots of debugging info)

* Parse: absorb, clusters, and weights
	Start, absorb(`absorb') clustervars(`clustervars') weight(`weighttype') weightvar(`weightvar')
	local absorb_keepvars = r(keepvars)
	local N_hdfe = r(N_hdfe)

* Check if we can save FEs
	forval g = 1/`N_hdfe' {
		mata: fe2local(`g')
		local targets "`targets'`target'"
	}
	if ("`targets'"!="") {
		Assert ("`partial'"==""), msg("hdfe error: partial() not allowed when saving fixed effects")
		local numvars : word count `varlist'
		Assert `numvars'==1 , msg("hdfe error: to save the fixed effects, you need to demean only one variable")
		local opt_savefe "save_fe(1)"
	}

* Keep relevant observations
	marksample touse, novar
	markout `touse' `varlist' `partial' `absorb_keepvars'
	qui keep if `touse'
	
* Keep relevant variables
	keep `varlist' `partial' `clustervars' `weightvar' `panelvar' `timevar' `uid' `absorb_keepvars'

* Drop singletons
	if ("`dropsingletons'"!="") DropSingletons, num_absvars(`N_hdfe')
	
* Construct Mata objects and auxiliary variables
	Precompute, ///
		keep(`varlist' `partial' `clustervars' `weightvar' `panelvar' `timevar' `uid') ///
		tsvars(`panelvar' `timevar')
	
* Compute e(df_a)
	EstimateDoF, dofadjustments(pairwise clusters continuous)
	* return list // what matters is r(kk) which will be e(df_a)
	local kk = r(kk)
	forval g = 1/`N_hdfe' {
		local df_a`g' = r(K`g') - r(M`g')
	}
	
* We don't need the FE variables (they are in mata objects now)
	*drop __FE*__

* Demean variables wrt to the fixed effects
	local opt varlist(`varlist' `partial') tol(`tolerance') maxiterations(`maxiterations') `options' `opt_savefe'
	if (`cores'>1) {
		DemeanParallel, `opt' self(hdfe) cores(`cores')
	}
	else {
		Demean, `opt'
	}

	if ("`opt_savefe'"!="") {
		Save, original_depvar(`varlist')
		local saved_fe = r(keepvars)
	}
	
	return scalar df_a = `kk'
	return scalar N_hdfe = `N_hdfe'
	forv g=1/`N_hdfe' {
		mata: fe2local(`g') // copies Mata structure into locals
		* Will inject the following with c_local:
		* ivars cvars target varname varlabel is_interaction is_cont_interaction is_bivariate is_mock levels
		return local hdfe`g' = "`varlabel'"
		return scalar df_a`g' = `df_a`g'' // `levels'
	}
* Clean up Mata objects
	Stop

* Deal with partial() option (Alternative: do as ivreg-partial: precompute inv x'x )
	if ("`partial'"!="") {
		tempvar resid
		_rmcoll `partial', forcedrop
		local partial = r(varlist)
		foreach var of local varlist {
			_regress `var' `partial' `weightexp' [`weighttype'`weightequal'`weightvar'], nohead notable
			_predict double `resid', resid
			qui replace `var' = `resid' // preserve labels
			drop `resid'
		}
		local numpartial : word count `partial'
		return scalar df_partial = `numpartial'
	}

	if ("`generate'"!="") {
		keep `varlist' `uid' `saved_fe'
		foreach var of local varlist {
			rename `var' `generate'`var'
		}

		tempfile output
		sort `uid'
		qui save "`output'"
		restore
		SafeMerge, uid(`uid') file("`output'") sample(`sample')
	}

}
if (_rc) {
	local rc = _rc
	Stop
	exit `rc'
}
end

* [SafeMerge: ADAPTED FROM THE ONE IN ESTIMATE.ADO]
* The idea of this program is to keep the sort order when doing the merges
cap pr drop SafeMerge
program define SafeMerge, eclass sortpreserve
syntax, uid(varname numeric) file(string) [sample(string)]
	* Merging gives us e(sample) and the FEs / AvgEs
	if ("`sample'"!="") {
		tempvar smpl
		merge 1:1 `uid' using "`file'", assert(master match) nolabel nonotes noreport gen(`smpl')
		gen byte `sample' = (`smpl'==3)
		drop `smpl' // redundant
	}
	else {
		merge 1:1 `uid' using "`file'", assert(master match) nolabel nonotes noreport nogen
	}
end

include "_common/Assert.ado"
include "_common/Debug.ado"
include "_common/Version.ado"

include "_hdfe/ConnectedGroups.ado"
include "_hdfe/GenerateID.ado"
include "_hdfe/AverageOthers.ado"
include "_hdfe/EstimateDoF.ado"

include "_hdfe/Start.ado"
include "_hdfe/ParseOneAbsvar.ado"
include "_hdfe/Precompute.ado"
include "_hdfe/Demean.ado"
include "_hdfe/DemeanParallel.ado"
include "_hdfe/ParallelInstance.ado"
include "_hdfe/Save.ado"
include "_hdfe/Stop.ado"
include "_hdfe/CheckCorrectOrder.ado"
include "_hdfe/DropSingletons.ado"
