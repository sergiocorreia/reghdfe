capture program drop Inner
program define Inner, eclass

* INITIAL CLEANUP
	ereturn clear // Clear previous results and drops e(sample)

* PARSE - inject opts with c_local, create Mata structure HDFE_S (use verbose>2 for details)
	Parse `0'
	assert !`savecache'
	assert !`usecache'

* PRESERVE (optional)
	preserve
	Debug, level(2) newline
	Debug, level(2) msg("(dataset preserved)")

* MEMORY REPORT - Store dataset size
	qui de, simple
	local old_mem = string(r(width) * r(N)  / 2^20, "%6.2f")
	local raw_n = c(N)
	local raw_k = c(k)

* CREATE UID - allows attaching e(sample) and the FE estimates into the restored dataset
	if (!`fast') {
		tempvar uid
		GenUID `uid'
	}

* COMPACT - Expand time and factor variables, and drop unused variables and obs.
	local original_depvar "`depvar'"
	Compact, basevars(`basevars') depvar(`depvar') indepvars(`indepvars') endogvars(`endogvars') instruments(`instruments') uid(`uid') timevar(`timevar') panelvar(`panelvar') weightvar(`weightvar') absorb_keepvars(`absorb_keepvars') clustervars(`clustervars') if(`if') in(`in') verbose(`verbose') vceextra(`vceextra')
	// Injects locals: depvar indepvars endogvars instruments expandedvars

* PRECOMPUTE MATA OBJECTS (means, counts, etc.)
	mata: map_init_keepvars(HDFE_S, "`expandedvars' `uid'") 	// Non-essential vars will be deleted (e.g. interactions of a clustervar)
	mata: map_precompute(HDFE_S)
	
	* Replace vceoption with the correct cluster names (e.g. if it's a FE or a new variable)
	if (`num_clusters'>0) {
		assert "`r(updated_clustervars)'"!=""
		local vceoption : subinstr local vceoption "<CLUSTERVARS>" "`r(updated_clustervars)'"
	}

* MEMORY REPORT
	Debug, level(2) msg("(dataset compacted: observations " as result "`raw_n' -> `c(N)'" as text " ; variables " as result "`raw_k' -> `c(k)'" as text ")")
	qui de, simple
	local new_mem = string(r(width) * r(N) / 2^20, "%6.2f")
	Debug, level(2) msg("(dataset compacted, c(memory): " as result "`old_mem'" as text "M -> " as result "`new_mem'" as text "M)")
	if (`verbose'>3) {
		di as text "(memory usage including mata:)"
		memory
		di as text ""
	}

* PREPARE - Compute untransformed tss, R2 of eqn w/out FEs
	Prepare, weightexp(`weightexp') depvar(`depvar') stages(`stages') model(`model') expandedvars(`expandedvars') vcetype(`vcetype') endogvars(`endogvars')
	* Injects tss, tss_`endogvar' (with stages), and r2c

* STORE UID - Used to add variables to original dataset: e(sample), mobility group, and FE estimates
	if (!`fast') mata: store_uid(HDFE_S, "`uid'")
	if (`fast') Debug, msg("(option {opt fast} specified; will not save e(sample))")

* BACKUP UNTRANSFORMED VARIABLES - If we are saving the FEs, we need to backup the untransformed variables
	if (`will_save_fe') {
		tempfile untransformed
		qui save "`untransformed'"
	}

* COMPUTE e(stats) - Summary statistics for the all the regression variables
	if ("`stats'"!="") {
		local tabstat_weight : subinstr local weightexp "[pweight" "[aweight"
		qui tabstat `expandedvars' `tabstat_weight' , stat(`stats') col(stat) save
		tempname statsmatrix
		matrix `statsmatrix' = r(StatTotal)
	}

* MAP_SOLVE() - WITHIN TRANFORMATION (note: overwrites variables)
	qui ds `expandedvars'
	local NUM_VARS : word count `r(varlist)'
	Debug, msg("(computing residuals for `NUM_VARS' variables)")
	mata: map_solve(HDFE_S, "`expandedvars'")

* STAGES SETUP - Deal with different stages
	assert "`stages'"!=""
	if ("`stages'"!="none") {
		Debug, level(1) msg(_n " {title:Stages to run}: " as result "`stages'" _n)
		* Need to backup some locals
		local backuplist groupvar fast will_save_fe depvar indepvars endogvars instruments original_depvar tss
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
		}
	}

* COMPUTE DOF
* NOTE: could move this to before backup untransformed variables for the stages=none case!!!
	mata: map_estimate_dof(HDFE_S, "`dofadjustments'", "`groupvar'") // requires the IDs
	assert e(df_a)<. // estimate_dof() only sets e(df_a); ereturn_dof() is for setting everything aferwards
	local kk = e(df_a) // we need this for the regression step
* DROP FE IDs - Except if they are also a clustervar or we are saving their respecting alphas
	mata: drop_ids(HDFE_S)

* REGRESS - Call appropiate wrapper (regress, avar, mwc for ols; ivreg2, ivregress for iv)
	ereturn clear
	if ("`stage'"=="none") Debug, level(2) msg("(running regresion: `model'.`ivsuite')")
	local wrapper "Wrapper_`subcmd'" // regress ivreg2 ivregress
	if ("`subcmd'"=="regress" & "`vcesuite'"=="avar") local wrapper "Wrapper_avar"
	if ("`subcmd'"=="regress" & "`vcesuite'"=="mwc") local wrapper "Wrapper_mwc"
	if (!inlist("`stage'","none", "iv")) local wrapper "Wrapper_avar" // Compatible with ivreg2
	Debug, level(3) msg(_n "call to wrapper:" _n as result "`wrapper', `options'")
	local opt_list
	local opts ///
		depvar indepvars endogvars instruments ///
		vceoption vcetype vcesuite ///
		kk suboptions showraw vceunadjusted first weightexp ///
		estimator twicerobust // Whether to run or not two-step gmm
	foreach opt of local opts {
		local opt_list `opt_list' `opt'(``opt'')
	}
	`wrapper', `opt_list'

* SAVE FE
	if (`will_save_fe') {
		local subpredict = e(predict) // used to recover the FEs
		SaveFE, model(`model') depvar(`depvar') untransformed(`untransformed') weightexp(`weightexp') subpredict(`subpredict')
	}

* FIX VARNAMES - Replace tempnames in the coefs table (run AFTER regress and BEFORE restore)
	* (e.g. __00001 -> L.somevar)
	tempname b
	matrix `b' = e(b)
	local backup_colnames : colnames `b'
	FixVarnames `backup_colnames'
	local newnames "`r(newnames)'"
	matrix colnames `b' = `newnames'
	ereturn repost b=`b', rename
	ereturn local depvar = "`original_depvar'" // Run after SaveFE

* (optional) Restore
	if inlist("`stage'","none", "iv") {
		restore
		Debug, level(2) newline
		Debug, level(2) msg("(dataset restored)")
		// TODO: Format alphas
	}

* (optional) Save mobility groups
	if ("`groupvar'"!="") mata: groupvar2dta(HDFE_S)

* (optional) Save alphas (fixed effect estimates)
	if (`will_save_fe') mata: alphas2dta(HDFE_S)

* (optional) Add e(sample)
	if (!`fast') {
		tempvar sample
		mata: esample2dta(HDFE_S, "`sample'")
		qui replace `sample' = 0 if `sample'==.
		la var `sample' "[HDFE Sample]"
		ereturn repost , esample(`sample')
		mata: drop_uid(HDFE_S)
	}

* POST ERETURN - Add e(...) (besides e(sample) and those added by the wrappers)	
	local opt_list
	local opts dofadjustments subpredict model stage stages subcmd cmdline vceoption equation_d original_absvars extended_absvars vcetype vcesuite tss r2c savefirst diopts weightvar gmm2s cue liml dkraay by level num_clusters clustervars timevar
	foreach opt of local opts {
		local opt_list `opt_list' `opt'(``opt'')
	}
	Post, `opt_list'

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
	Attach, notes(`notes') statsmatrix(`statsmatrix') summarize_quietly(`summarize_quietly') // Attach only once, not per stage

* CLEANUP
	mata: mata drop HDFE_S // cleanup
end
