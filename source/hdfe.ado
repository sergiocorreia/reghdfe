
// Mata code is first, then main reghdfe.ado, then auxiliary .ado files
clear mata
include "mata/map.mata"

capture program drop hdfe
program define hdfe, eclass

* Set Stata version
	version `=clip(c(version), 11.2, 14.1)'

* Intercept version
	cap syntax, version
	if !c(rc) {
		Version
		exit
	}

* Parse options specific to hdfe.ado: partial sample generate clear clustervars
	syntax anything(everything) [fw aw pw/], [*] ///
		[CACHE(string) VCE(string) CLUSTER(string) FAST] /// Disabled options
		Absorb(string) [ ///
		/// PARTIAL(varlist numeric) /// Additional regressors besides those in absorb()
		CLUSTERVars(string) /// Used to estimate the DoF
		Generate(name) SAMPLE(name) ///
		CLEAR KEEPIDs KEEPVars(varlist) ///
		 ///
		]

	Assert "`cache'"=="", msg("invalid option {cmd:cache()}")
	Assert "`vce'"=="", msg("invalid option {cmd:vce()}")
	Assert "`fast'"=="", msg("invalid option {cmd:fast()}")
	Assert "`cluster'"=="", msg("invalid option {cmd:cluster()}, perhaps you meant {opt clusterv:ars()}?")
	Assert ("`generate'"!="") + ("`clear'"!="") == 1 , ///
		msg("hdfe error: you need to specify one and only one of the following options: clear generate(...)")

	if ("`clear'"!="") Assert "`sample'"=="", msg("option {cmd:sample()} not compatible with {cmd:clear}")
	if ("`generate'"!="") Assert "`keepids'"=="", msg("option {cmd:keepids} not compatible with {cmd:generate()}")
	if ("`generate'"!="") Assert "`keepvars'"=="", msg("option {cmd:keepvars()} not compatible with {cmd:generate()}")

	cap drop __ID*__
	cap drop __CL*__
	
	if ("`keepvars'"!="") {
		qui ds `keepvars'
		local keepvars `r(varlist)'
	}

	if ("`weight'"!="") local weight_part [`weight'=`exp']
	local 0 `partial' `anything' `weight_part' , absorb(`absorb') `options' ///
		cluster(`clustervars') cache(save, keepvars(`keepvars'))

* INITIAL CLEANUP
	ereturn clear // Clear previous results and drops e(sample)

* From now on, we will pollute the Mata workspace, so wrap this in case of error
	cap noi {

* PARSE - inject opts with c_local, create Mata structure HDFE_S (use verbose>2 for details)
	Parse `0'
	assert `savecache'
	Assert !`will_save_fe', msg("savecache disallows saving FEs")

	if ("`clear'"=="") {
		preserve
		tempvar uid
		GenUID `uid'
	}

* PROBLEM:
	* I can translate L(1/2).x into __L__x __L2__x
	* But how can I translate i.x if I don't have the original anymore?

* SOLUTION
	* The cache option of ExpandFactorVariables (called from Compact.ado)

* COMPACT - Expand time and factor variables, and drop unused variables and obs.
	Compact, basevars(`basevars') depvar(`depvar' `indepvars') uid(`uid') timevar(`timevar') panelvar(`panelvar') weightvar(`weightvar') weighttype(`weighttype') ///
		absorb_keepvars(`absorb_keepvars') clustervars(`clustervars') ///
		if(`if') in(`in') verbose(`verbose') vceextra(`vceextra') savecache(1) more_keepvars(`keepvars')
	// Injects locals: depvar indepvars endogvars instruments expandedvars cachevars

* PRECOMPUTE MATA OBJECTS (means, counts, etc.)
	mata: map_init_keepvars(HDFE_S, "`expandedvars' `uid' `cachevars' `keepvars'") 	// Non-essential vars will be deleted (e.g. interactions of a clustervar)
	mata: map_precompute(HDFE_S)
	global updated_clustervars = "`r(updated_clustervars)'"

* Store UID in Mata so we can then attach the variables without an expensive merge
	if ("`clear'"=="") mata: store_uid(HDFE_S, "`uid'")
	
* COMPUTE DOF
	if (`timeit') Tic, n(62)
	mata: map_estimate_dof(HDFE_S, "`dofadjustments'", "`groupvar'") // requires the IDs
	if (`timeit') Toc, n(62) msg(estimate dof)
	assert e(df_a)<. // estimate_dof() only sets e(df_a); map_ereturn_dof() is for setting everything aferwards
	local kk = e(df_a) // we need this for the regression step

* MAP_SOLVE() - WITHIN TRANFORMATION (note: overwrites variables)
	qui ds `expandedvars'
	local NUM_VARS : word count `r(varlist)'
	Debug, msg("(computing residuals for `NUM_VARS' variables)")
	
* Build newvar string (trick: if generate is empty, then newvars==expandedvars)
	foreach var of local expandedvars {
		local newvars `newvars' `generate'`var'
	}

	local restore_dta = ("`clear'"=="")
	mata: map_solve(HDFE_S, "`expandedvars'", "`newvars'", `restore_dta')

	foreach var of local expandedvars {
		local label : char `generate'`var'[name]
		if ("`label'"=="") local label `var'
		la var `generate'`var' "Residuals: `label'"
	}

	if ("`groupvar'"!="") mata: groupvar2dta(HDFE_S, `restore_dta')

	if ("`sample'"!="") {
		mata: esample2dta(HDFE_S, "`sample'")
		qui replace `sample' = 0 if `sample'==.
		la var `sample' "[HDFE Sample]"
		mata: drop_uid(HDFE_S)
	}

* Absorbed-specific returns
	* e(N_hdfe) e(N_hdfe_extended) e(mobility)==M e(df_a)==K-M
	* e(M#) e(K#) e(M#_exact) e(M#_nested) -> for #=1/e(N_hdfe_extended)
	mata: map_ereturn_dof(HDFE_S)
	local N_hdfe = e(N_hdfe)
	ereturn local cmd = "hdfe"
	ereturn local extended_absvars "`extended_absvars'"
	ereturn local absvars "`original_absvars'"

* Cleanup
	} // cap noi
	local rc = c(rc)
	cap mata: mata drop HDFE_S // overwrites c(rc)
	*cap mata: mata drop varlist_cache
	cap global updated_clustervars

	local keys absorb N_hdfe original_absvars extended_absvars vce vceoption vcetype ///
		vcesuite vceextra num_clusters clustervars bw kernel dkraay kiefer twicerobust ///
		reghdfe_cache // cache_obs
	foreach key of local keys {
		char _dta[`key']
	}

	if ("`keepids'"=="" | `rc') cap drop __ID*__
	if ("`keepids'"=="" | `rc') cap drop __CL*__

	if (`rc') exit `rc'
end

// -------------------------------------------------------------------------------------------------
include "common/Assert.ado"
include "common/Debug.ado"
include "common/Version.ado"
include "common/Tic.ado"
include "common/Toc.ado"
	include "internal/Parse.ado"
		include "internal/ParseCache.ado"
		include "internal/ParseIV.ado"
		include "internal/ParseStages.ado"
		include "internal/ParseVCE.ado"
		include "internal/ParseAbsvars.ado"
		include "internal/ParseDOF.ado"
		include "internal/ParseImplicit.ado"
	include "internal/GenUID.ado"
	include "internal/Compact.ado"
		include "internal/ExpandFactorVariables.ado"
		include "internal/GenerateID.ado"
// -------------------------------------------------------------------------------------------------
