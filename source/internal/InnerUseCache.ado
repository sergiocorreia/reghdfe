capture program drop InnerUseCache
program define InnerUseCache, eclass

* INITIAL CLEANUP
	ereturn clear // Clear previous results and drops e(sample)

* PARSE - inject opts with c_local, create Mata structure HDFE_S (use verbose>2 for details)
	Parse `0'
	assert `usecache'

	local original_depvar "`depvar'"

* Match "L.price" --> __L__price
* Expand factor and time-series variables
* (based on part of Compact.ado)
	local expandedvars
	local sets depvar indepvars endogvars instruments // depvar MUST be first
	Debug, level(4) newline
	Debug, level(4) msg("{title:Expanding factor and time-series variables:}")
	foreach set of local sets {
		local varlist ``set''
		local `set' // empty
		if ("`varlist'"=="") continue
		fvunab factors : `varlist', name("error parsing `set'")
		foreach factor of local factors {
			mata: st_local("var", asarray(varlist_cache, "`factor'"))
			Assert "`var'"!="", msg("couldn't find the match of {res}`factor'{error} in the cache (details: set=`set'; factors=`factors')")
			local `set' ``set'' `var'
		}
		local expandedvars `expandedvars' ``set''
	}

* Replace vceoption with the correct cluster names (e.g. if it's a FE or a new variable)
	if (`num_clusters'>0) {
		assert "$updated_clustervars"!=""
		local vceoption : subinstr local vceoption "<CLUSTERVARS>" "$updated_clustervars"
	}

* LEVEL check
	if ("`by'"!="") {
		cap cou if `by'==`level'
		Assert r(N)>0 , msg("reghdfe by() error: there are no cases where `by'==`level'")
		local cond if `by'==`level'
	}

* PREPARE - Compute untransformed tss, R2 of eqn w/out FEs
	if ("`by'"=="") {
		mata: st_local("tss", strofreal(asarray(tss_cache, "`depvar'")))
		Assert `tss'<., msg("tss of depvar `depvar' not found in cache")
		foreach var of local endogvars {
			mata: st_local("tss_`var'", strofreal(asarray(tss_cache, "`var'")))
		}
		local r2c = . // BUGBUG!!!
	}

* STAGES SETUP - Deal with different stages
	assert "`stages'"!=""
	if ("`stages'"!="none") {
		Debug, level(1) msg(_n " {title:Stages to run}: " as result "`stages'" _n)
		* Need to backup some locals
		local backuplist residuals groupvar fast will_save_fe depvar indepvars endogvars instruments original_depvar tss
		foreach loc of local backuplist {
			local backup_`loc' ``loc''
		}

		local num_stages : word count `stages'
		local last_stage : word `num_stages' of `stages'
		assert "`last_stage'"=="iv"
	}

* STAGES LOOPS
foreach stage of local stages {
Assert inlist("`stage'", "none", "iv", "first", "ols", "reduced", "acid")
local lhs_endogvars = cond("`stage'"=="first", "`backup_endogvars'", "<none>")
local i_endogvar = cond("`stage'"=="first", "0", "")
foreach lhs_endogvar of local lhs_endogvars {

	if ("`stage'"!="none") {
		* Start with backup values
		foreach loc of local backuplist {
			local `loc' `backup_`loc''
		}

		if ("`stage'"=="ols") {
			local indepvars `indepvars' `endogvars'
		}
		else if ("`stage'"=="reduced") {
			local indepvars `indepvars' `instruments'
		}
		else if ("`stage'"=="acid") {
			local indepvars `indepvars' `endogvars' `instruments'
		}
		else if ("`stage'"=="first") {
			local ++ i_endogvar
			local tss = `tss_`lhs_endogvar''
			assert `tss'<.
			local depvar `lhs_endogvar'
			local indepvars `indepvars' `instruments'
			local original_depvar : char `depvar'[name]
			if ("`original_depvar'"=="") local original_depvar `depvar' 
		}

		if ("`stage'"!="iv") {
			local fast 1
			local will_save_fe 0
			local vcesuite avar
			local endogvars
			local instruments
			local groupvar
			local residuals
		}
	}

* COMPUTE DOF
* NOTE: could move this to before backup untransformed variables for the stages=none case!!!
	mata: map_estimate_dof(HDFE_S, "`dofadjustments'", "`groupvar'", "`cond'") // requires the IDs
	assert e(df_a)<. // estimate_dof() only sets e(df_a); ereturn_dof() is for setting everything aferwards
	local kk = e(df_a) // we need this for the regression step

* REGRESS - Call appropiate wrapper (regress, avar, mwc for ols; ivreg2, ivregress for iv)
	ereturn clear
	if ("`stage'"=="none") Debug, level(2) msg("(running regresion: `model'.`ivsuite')")
	local wrapper "Wrapper_`subcmd'" // regress ivreg2 ivregress
	if ("`subcmd'"=="regress" & "`vcesuite'"=="avar") local wrapper "Wrapper_avar"
	if ("`subcmd'"=="regress" & "`vcesuite'"=="mwc") local wrapper "Wrapper_mwc"
	if (!inlist("`stage'","none", "iv")) local wrapper "Wrapper_avar" // Compatible with ivreg2
	local opt_list
	local opts /// cond // BUGUBG: Add by() (cond) options
		depvar indepvars endogvars instruments ///
		vceoption vcetype vcesuite ///
		kk suboptions showraw vceunadjusted first weightexp ///
		estimator twicerobust // Whether to run or not two-step gmm
	foreach opt of local opts {
		local opt_list `opt_list' `opt'(``opt'')
	}
	Debug, level(3) msg(_n "call to wrapper:" _n as result "`wrapper', `opt_list'")
	`wrapper', `opt_list'

* COMPUTE AND STORE RESIDS (based on SaveFE.ado)
	local drop_resid_vector
	if ("`residuals'"!="") {
		local drop_resid_vector drop_resid_vector(0)
		local subpredict = e(predict)
		local score = cond("`model'"=="ols", "score", "resid")
		if e(df_m)>0 {
			`subpredict' double `residuals', `score' // equation: y = xb + d + e, we recovered "e"
		}
		else {
			gen double `residuals' = `depvar'
		}
		// No need to store in Mata
	}

* FIX VARNAMES - Replace tempnames in the coefs table (run AFTER regress and BEFORE restore)
	* (e.g. __00001 -> L.somevar)
	tempname b
	matrix `b' = e(b)
	local backup_colnames : colnames `b'
	FixVarnames `backup_colnames'
	local newnames "`r(newnames)'"
	matrix colnames `b' = `newnames'
	ereturn local depvar = "`original_depvar'" // Run after SaveFE

* POST ERETURN - Add e(...) (besides e(sample) and those added by the wrappers)	
	local opt_list
	local opts dofadjustments subpredict model stage stages subcmd cmdline vceoption equation_d original_absvars extended_absvars vcetype vcesuite tss r2c savefirst diopts weightvar gmm2s cue liml dkraay by level num_clusters clustervars timevar
	foreach opt of local opts {
		local opt_list `opt_list' `opt'(``opt'')
	}
	Post, `opt_list' coefnames(`b')

* REPLAY - Show the regression table	
	if ("`stage'"!="none") Debug, level(0) msg(_n "{title:Stage: `stage'}" _n)
	if ("`lhs_endogvar'"!="<none>") Debug, level(0) msg("{title:Endogvar: `original_depvar'}")
	Replay

* STAGES - END
	if (!inlist("`stage'","none", "iv")) {
		local estimate_name reghdfe_`stage'`i_endogvar'
		local stored_estimates `stored_estimates' `estimate_name'
		local cmd estimates store `estimate_name', nocopy
		Debug, level(2) msg(" - Storing estimate: `cmd'")
		`cmd'
	}
	else if ("`stage'"=="iv") {
		* On the last stage, save list of all stored estimates
		assert "`stored_estimates'"!=""
		ereturn `hidden' local stored_estimates = "`stored_estimates'"
	}
} // lhs_endogvar
} // stage

* ATTACH - Add e(stats) and e(notes)
	if ("`by'"=="") {
		cap conf matrix reghdfe_statsmatrix
		if (!c(rc)) local statsmatrix reghdfe_statsmatrix
		Attach, notes(`notes') statsmatrix(`statsmatrix') summarize_quietly(`summarize_quietly') // Attach only once, not per stage
	}
end
