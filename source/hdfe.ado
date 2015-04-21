
// Mata code is first, then main hdfe.ado, then auxiliary .ado files
clear mata
include "mata/map.mata"

capture program drop hdfe
program define hdfe, rclass

* Set Stata version
	version `=clip(`c(version)', 11.2, 13.1)' // 11.2 minimum, 13+ preferred

* Intercept version calls
	cap syntax, version
	if !c(rc) {
		Version
		exit
	}

* Parse
	syntax varlist(numeric) [if] [in] [fw aw pw/] , ///
	/// Main Options ///
		Absorb(string) ///
		[ ///
		PARTIAL(varlist numeric) /// Additional regressors besides those in absorb()
		SAMPLE(name) ///
		Generate(name) CLEAR /// Replace dataset, or just add new variables
		GROUPVAR(name) /// Variable that will contain the first connected group between FEs
		CLUSTERVARs(varlist numeric fv max=10) /// Used to estimate the DoF
	/// Optimization /// Defaults are handled within Mata
		GROUPsize(string) /// Process variables in groups of #
		TRANSFORM(string) ///
		ACCELeration(string) ///
		Verbose(string) ///
		TOLerance(string) ///
		MAXITerations(string) ///
		KEEPSINGLETONS(string) /// Only use this option for debugging
		SUBCMD(string) /// Regression package
		] [*] // Remaining options 

* Time/panel variables
	local timevar `_dta[_TStvar]'
	local panelvar `_dta[_TSpanel]'

* Validation
	local clustervars : subinstr local clustervars "i." "", all // Remove i. prefixes
	if ("`options'"!="") di as error "unused options: `options'"
	if ("`sample'"!="") confirm new variable `sample'
	Assert ("`generate'"!="") + ("`clear'"!="") == 1 , ///
		msg("hdfe error: you need to specify one and only one of the following options: clear generate(...)")
	Assert "`: list varlist & partial'"=="", ///
		msg("variables in varlist cannot appear in partial()")
	if ("`weight'"!="") {
		local weightvar `exp'
		local weighttype `weight'
		confirm var `weightvar', exact // just allow simple weights
	}
	if ("`group'"!="") confirm new var `group'

* From now on, we will pollute the Mata workspace, so wrap this in case of error
	cap noi {

* Create Mata structure
	ParseAbsvars `absorb' // Stores results in r()
	// return list
	mata: HDFE_S = map_init() // Reads results from r()
	// return list
	local save_fe = r(save_fe)

	if ("`weightvar'"!="") mata: map_init_weights(HDFE_S, "`weightvar'", "`weighttype'")
	* String options
	local optlist transform acceleration clustervars panelvar timevar
	foreach opt of local optlist {
		if ("``opt''"!="") mata: map_init_`opt'(HDFE_S, "``opt''")
	}
	* Numeric option
s	local optlist groupsize verbose tolerance maxiterations keepsingletons
	foreach opt of local optlist {
		if ("``opt''"!="") mata: map_init_`opt'(HDFE_S, ``opt'')
	}

* (Optional) Preserve
	if ("`generate'"!="" | `save_fe') {
		tempvar uid
		local uid_type = cond(c(N)>c(maxlong), "double", "long")
		gen `uid_type' `uid' = _n // Useful for later merges
		la var `uid' "[UID]"
		preserve
	}

	if ("`generate'"!="") {
		foreach var of varlist `varlist' {
			confirm new var `generate'`var'
			local newvars `newvars' `generate'`var'
		}
	}

* Precompute Mata objects
	mata: map_init_keepvars(HDFE_S, "`varlist' `partial' `uid'") // Non-essential vars will be deleted
	mata: map_precompute(HDFE_S)

* Compute e(df_a)
	mata: map_estimate_dof(HDFE_S, "pairwise clusters continuous", "`group'")
	//return list

* (Optional) Drop IDs, unless i) they are also clusters or ii) we want to save fe
	// TODO

* (Optional) Need to backup dataset if we want to save FEs
	if (`save_fe') {
		tempfile untransformed
		qui save "`untransformed'"
	}

* Within Transformation
	mata: map_solve(HDFE_S, "`varlist'", "`newvars'", "`partial'")

* Run regression
	if ("`subcmd'"!="") {
		// TODO
		// `subcmd' ..  `varlist' `options' ...
	}
	else if (`save_fe') { // Need to regress before predicting
		regress `varlist', noheader notable // qui  // BUGBUG: _regress?
	}

* (Optional) Save FEs
	if (`save_fe') {
		tempvar resid
		predict double `resid', resid
		keep `uid' `resid'
		tempfile transformed
		qui save "`transformed'"

		qui use "`untransformed'"
		erase "`untransformed'"

		merge 1:1 `uid' using "`transformed'", assert(match) nolabel nonotes noreport nogen
		erase "`transformed'"
		tempvar resid_d
		predict double `resid_d', resid
		if ("`weightvar'"!="") local tmp_weight "[fw=`weightvar']" // summarize doesn't work with pweight
		su `resid_d' `tmp_weight', mean
		qui replace `resid_d' = `resid_d' - r(mean)
		tempvar d
		gen double `d' = `resid_d' - `resid'
		//clonevar dd = `d'
		mata: map_solve(HDFE_S, "`d'", "", "`partial'", 1) // Save FE (should fail if partial is set)
		//regress dd __hdfe*, nocons
	}

* (Optional) Tempsave, restore and merge with
	//if ("`esample'"!="" | `save_fe' | "`generate'"!="")
	//maso menos xq las variables transformadas ya las chanque (estaban en transformed!)
	//esta parte es medio rara
	// ...
	// need to add vars if i) i want e(sample), ii) i want to merge the FEs, iii) i want to merge

	if ("`generate'"!="" | `save_fe') restore

* Cleanup after an error
	} // cap noi
	if c(rc) {
		local rc = c(rc)
		cap mata: mata drop HDFE_S // overwrites c(rc)
		exit `rc'
	}
end
