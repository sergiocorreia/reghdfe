*! hdfe 1.1 22mar2015
*! Sergio Correia (sergio.correia@duke.edu)
* (built from multiple source files using build.py)

include _mata/reghdfe.mata

cap pr drop hdfe
program define hdfe, rclass
	local version `clip(`c(version)', 11.2, 13.1)' // 11.2 minimum, 13+ preferred
	qui version `version'

	syntax varlist [if] [in] [fweight aweight pweight/] , Absorb(string) ///
		[CORES(integer 1)]
		[CLUSTERvars(string) Verbose(integer 0) TOLerance(real 1e-7) MAXITerations(integer 10000)] [GENerate(name)]
	if ("`weight'"!="") {
		local weightvar `exp'
		conf var `weightvar' // just allow simple weights
		local weighttype `weight'
	}

* Intercept multiprocessor/parallel calls
	cap syntax, instance [*]
	local rc = _rc
	 if (`rc'==0) {
		ParallelInstance, `options'
		exit
	}
	
* Preserve if asked to
	if ("`generate'"!="") {
		tempvar uid
		gen double `uid' = _n
		preserve
	}
	
* Clear previous errors
	cap Stop

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
	
* Keep relevant observations
	marksample touse, novar
	markout `touse' `varlist' `absorb_keepvars'
	qui keep if `touse'
	
* Keep relevant variables
	keep `varlist' `clustervars' `weightvar' `panelvar' `timevar' `absorb_keepvars' `uid'
	
* Construct Mata objects and auxiliary variables
	Precompute, keep(`varlist' `clustervars' `weightvar' `panelvar' `timevar' `uid') tsvars(`panelvar' `timevar')

* Compute e(df_a)
	EstimateDoF, dofadjustments(pairwise clusters continuous)
	* return list // what matters is r(kk) which will be e(df_a)
	local kk = r(kk)
	
* Demean variables wrt to the fixed effects
	local opt varlist(`varlist') tol(`tolerance') maxiterations(`maxiterations') // Other maximize/parallel options
	if (`cores'>1) {
		DemeanParallel, `opt' self(hdfe) cores(`cores')
	}
	else {
		Demean, `opt'	
	}
	
	return scalar df_a = `kk'
	return scalar N_hdfe = `N_hdfe'
	forv g=1/`N_hdfe' {
		mata: fe2local(`g') // copies Mata structure into locals
		* Will inject the following with c_local:
		* ivars cvars target varname varlabel is_interaction is_cont_interaction is_bivariate is_mock levels
		return local hdfe`g' = "`varlabel'"
		return scalar df_a`g' = `levels'
	}
	
* Clean up Mata objects
	Stop

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

include "_common/Assert.ado"
include "_common/Debug.ado"

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
