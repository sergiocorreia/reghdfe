capture program drop InnerSaveCache
program define InnerSaveCache, eclass
* (note: based on Inner.ado)

* INITIAL CLEANUP
	ereturn clear // Clear previous results and drops e(sample)

* PARSE - inject opts with c_local, create Mata structure HDFE_S (use verbose>2 for details)
	Parse `0'
	assert `savecache'
	Assert !`will_save_fe', msg("savecache disallows saving FEs")

* PROBLEM:
	* I can translate L(1/2).x into __L__x __L2__x
	* But how can I translate i.x if I don't have the original anymore?

* SOLUTION
	* The cache option of ExpandFactorVariables (called from Compact.ado)

* COMPACT - Expand time and factor variables, and drop unused variables and obs.
	Compact, basevars(`basevars') depvar(`depvar' `indepvars') uid(`uid') timevar(`timevar') panelvar(`panelvar') weightvar(`weightvar') ///
		absorb_keepvars(`absorb_keepvars') clustervars(`clustervars') ///
		if(`if') in(`in') verbose(`verbose') vceextra(`vceextra') savecache(1) more_keepvars(`keepvars')
	// Injects locals: depvar indepvars endogvars instruments expandedvars cachevars

* PRECOMPUTE MATA OBJECTS (means, counts, etc.)
	mata: map_init_keepvars(HDFE_S, "`expandedvars' `uid' `cachevars' `by' `keepvars'") 	// Non-essential vars will be deleted (e.g. interactions of a clustervar)
	mata: map_precompute(HDFE_S)
	global updated_clustervars = "`r(updated_clustervars)'"
	
* STORE BY LEVELS
	if ("`by'"!="") {
		qui levelsof `by'
		ereturn local levels `r(levels)'
	}

* PREPARE - Compute untransformed tss *OF ALL THE VARIABLES*
	mata: tss_cache = asarray_create()
	mata: asarray_notfound(tss_cache, .)
	local tmpweightexp = subinstr("`weightexp'", "[pweight=", "[aweight=", 1)
	if ("`by'"=="") {
		foreach var of local expandedvars {
			qui su `var' `tmpweightexp' // BUGBUG: Is this correct?!
			local tss = r(Var)*(r(N)-1)
			mata: asarray(tss_cache, "`var'", `tss')
		}		
	}
	else {
		// ...
	}
	*NOTE: r2c is too slow and thus won't be saved
	*ALTERNATIVE: Allow a varlist of the form (depvars) (indepvars) and only compute for those

* COMPUTE e(stats) - Summary statistics for the all the regression variables
	if ("`stats'"!="" & "`by'"=="") {
		Stats `expandedvars', weightexp(`weightexp') stats(`stats') statsmatrix(reghdfe_statsmatrix)
	}

* MAP_SOLVE() - WITHIN TRANFORMATION (note: overwrites variables)
	qui ds `expandedvars'
	local NUM_VARS : word count `r(varlist)'
	Debug, msg("(computing residuals for `NUM_VARS' variables)")
	mata: map_solve(HDFE_S, "`expandedvars'")

* This was in -parse- but we are dropping observations through the code
	char _dta[cache_obs] `c(N)'

end
