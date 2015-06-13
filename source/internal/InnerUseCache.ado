capture program drop InnerUseCache
program define InnerUseCache, eclass

* INITIAL CLEANUP
	ereturn clear // Clear previous results and drops e(sample)

* PARSE - inject opts with c_local, create Mata structure HDFE_S (use verbose>2 for details)
	Parse `0'
	assert `usecache'
	if (`timeit') Tic, n(50)

	foreach cat in depvar indepvars endogvars instruments {
		local original_`cat' "``cat''"
	}

* Match "L.price" --> __L__price
* Expand factor and time-series variables
* (based on part of Compact.ado)
	if (`timeit') Tic, n(52)
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
	if (`timeit') Toc, n(52) msg(fix names)

* Replace vceoption with the correct cluster names (e.g. if it's a FE or a new variable)
	if (`num_clusters'>0) {
		assert "$updated_clustervars"!=""
		local vceoption : subinstr local vceoption "<CLUSTERVARS>" "$updated_clustervars"
	}

* PREPARE - Compute untransformed tss, R2 of eqn w/out FEs
	if (`timeit') Tic, n(54)
	mata: st_local("tss", asarray(tss_cache, "`depvar'"))
	Assert `tss'<., msg("tss of depvar `depvar' not found in cache")
	foreach var of local endogvars {
		mata: st_local("tss_`var'", asarray(tss_cache, "`var'"))
	}
	local r2c = . // BUGBUG!!!
	if (`timeit') Toc, n(54) msg(use cached tss)

 * COMPUTE DOF - Already precomputed in InnerSaveCache.ado
	if (`timeit') Tic, n(62)
	mata: map_ereturn_dof(HDFE_S) // this gives us e(df_a)==`kk', which we need
	assert e(df_a)<.
	local kk = e(df_a) // we need this for the regression step
	if (`timeit') Toc, n(62) msg(load dof estimates)

* STAGES SETUP - Deal with different stages
	assert "`stages'"!=""
	if ("`stages'"!="none") {
		Debug, level(1) msg(_n "{title:Stages to run}: " as result "`stages'")
		* Need to backup some locals
		local backuplist residuals groupvar fast will_save_fe depvar indepvars endogvars instruments original_depvar tss suboptions
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
			local endogvars
			local instruments
			local groupvar
			local residuals
			local suboptions `stage_suboptions'
		}
	}

* REGRESS - Call appropiate wrapper (regress, avar, mwc for ols; ivreg2, ivregress for iv)
	ereturn clear
	if ("`stage'"=="none") Debug, level(2) msg("(running regresion: `model'.`ivsuite')")
	local wrapper "Wrapper_`subcmd'" // regress ivreg2 ivregress
	if ("`subcmd'"=="regress" & "`vcesuite'"=="avar") local wrapper "Wrapper_avar"
	if ("`subcmd'"=="regress" & "`vcesuite'"=="mwc") local wrapper "Wrapper_mwc"
	if (!inlist("`stage'","none", "iv")) {
		if ("`vcesuite'"=="default") local wrapper Wrapper_regress
		if ("`vcesuite'"!="default") local wrapper Wrapper_`vcesuite'
	}
	local opt_list
	local opts /// cond // BUGUBG: Add by() (cond) options
		depvar indepvars endogvars instruments ///
		vceoption vcetype ///
		kk suboptions ffirst weightexp ///
		estimator twicerobust /// Whether to run or not two-step gmm
		num_clusters clustervars // Used to fix e() of ivreg2 first stages
	foreach opt of local opts {
		local opt_list `opt_list' `opt'(``opt'')
	}
	Debug, level(3) msg(_n "call to wrapper:" _n as result "`wrapper', `opt_list'")
	if (`timeit') Tic, n(66)
	`wrapper', `opt_list'
	if (`timeit') Toc, n(66) msg(regression)

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

* (optional) Save mobility groups (note: group vector will stay on HDFE_S)
	if ("`groupvar'"!="") mata: groupvar2dta(HDFE_S, 0)

* FIX VARNAMES - Replace tempnames in the coefs table (run AFTER regress)
	* (e.g. __00001 -> L.somevar)
	if (`timeit') Tic, n(68)
	tempname b
	matrix `b' = e(b)
	local backup_colnames : colnames `b'
	FixVarnames `backup_colnames'
	local newnames "`r(newnames)'"
	matrix colnames `b' = `newnames'
	ereturn local depvar = "`original_depvar'" // Run after SaveFE
	if (`timeit') Toc, n(68) msg(fix varnames)

* POST ERETURN - Add e(...) (besides e(sample) and those added by the wrappers)	
	local opt_list
	local opts dofadjustments subpredict model stage stages subcmd cmdline vceoption equation_d original_absvars extended_absvars vcetype vcesuite tss r2c savestages diopts weightvar estimator dkraay by level num_clusters clustervars timevar backup_original_depvar original_indepvars original_endogvars original_instruments
	foreach opt of local opts {
		local opt_list `opt_list' `opt'(``opt'')
	}
	if (`timeit') Tic, n(69)
	Post, `opt_list' coefnames(`b')
	if (`timeit') Toc, n(69) msg(Post)

* REPLAY - Show the regression table	
	Replay

* STAGES - END
	if (`timeit') Tic, n(70)
	if (!inlist("`stage'","none", "iv") & `savestages') {
		local estimate_name reghdfe_`stage'`i_endogvar'
		local stored_estimates `stored_estimates' `estimate_name'
		local cmd estimates store `estimate_name', nocopy
		Debug, level(2) msg(" - Storing estimate: `cmd'")
		`cmd'
	}
	else if ("`stage'"=="iv") {
		* On the last stage, save list of all stored estimates
		if ("`stored_estimates'"!="") ereturn `hidden' local stored_estimates = "`stored_estimates'"
	}
	if (`timeit') Toc, n(70) msg(store estimates if needed)

} // lhs_endogvar
} // stage

* ATTACH - Add e(stats) and e(notes)
	if ("`stats'"!="") {
		if (`timeit') Tic, n(71)
		tempname statsmatrix
		Stats `expandedvars', weightexp(`weightexp') stats(`stats') statsmatrix(`statsmatrix') usecache
		// stats() will be ignored
		if (`timeit') Tic, n(71) msg(Stats.ado)
	}
	if (`timeit') Tic, n(72)
	Attach, notes(`notes') statsmatrix(`statsmatrix') summarize_quietly(`summarize_quietly') // Attach only once, not per stage
	if (`timeit') Toc, n(72) msg(Attach.ado)

	if (`timeit') Toc, n(50) msg([TOTAL])
end
