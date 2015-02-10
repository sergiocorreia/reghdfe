*! version 1.3.1 09feb2015
*! By Sergio Correia (sergio.correia@duke.edu)
* (built from multiple source files using build.py)
program define reghdfe
	local version `=clip(`c(version)', 11.2, 13.1)' // 11.2 minimum, 13+ preferred
	qui version `version'

	if replay() {
		if (`"`e(cmd)'"'!="reghdfe") error 301
		* Undocumented feature for debugging
		cap syntax, ALTernative
		if ("`alternative'"!="") {
			AlternativeCMD
			exit
		}
		else {	
			Replay `0'
		}
	}
	else {
		* Estimate, and then clean up Mata in case of failure
		mata: st_global("reghdfe_pwd",pwd())
		cap noi Estimate `0'
		if (_rc) {
			local rc = _rc
			reghdfe_absorb, step(stop)
			exit `rc'
		}
	}
end


// -------------------------------------------------------------
// Transform data and run the regression
// -------------------------------------------------------------
program define Estimate, eclass
/* Notation
	__FE1__   Fixed effect categories
	__Z1__    Fixed effect coefficients (estimates)
	__W1__    AvgE transformed variables (avg of depvar by category) */

**** PART I - PREPARE DATASET FOR REGRESSION ****

* 1) Parse main options
	reghdfe_absorb, step(stop) // clean Mata leftovers before running -Parse-
	Parse `0' // save all arguments into locals (verbose>=3 shows them)
	local sets depvar indepvars endogvars instruments // depvar MUST be first

* 2) Parse identifiers (absorb variables, avge, clustervar)
	reghdfe_absorb, step(start) absorb(`absorb') over(`over') avge(`avge') clustervar1(`clustervar1') weight(`weight') weightvar(`weightvar')
	// In this step, it doesn't matter if the weight is FW or AW
	local N_hdfe = r(N_hdfe)
	local N_avge = r(N_avge)
	local absorb_keepvars "`r(keepvars)'" // Vars used in hdfe,avge,cluster
	local RAW_N = c(N)
	local RAW_K = c(k)
	qui de, simple
	local old_mem = string(r(width) * r(N)  / 2^20, "%6.2f") // In MBs

* 3) Preserve
if ("`usecache'"!="") {
	local uid __uid__
}
else {
	tempvar uid
	local uid_type = cond(`RAW_N'>c(maxlong), "double", "long")
	gen `uid_type' `uid' = _n // Useful for later merges
	la var `uid' "[UID]" // So I can recognize it in -describe-
}

	if (`savingcache') {
		cap drop __uid__
		rename `uid' __uid__
		local uid __uid__
		local handshake = int(uniform()*1e8)
		char __uid__[handshake] `handshake'
		char __uid__[tolerance] `tolerance'
		char __uid__[maxiterations] `maxiterations'
	}

	preserve
	Debug, msg("(dataset preserved)") level(2)

* 4) Drop unused variables
	local exp "= `weightvar'"
	marksample touse, novar // Uses -if- , -in- ; -weight-? and -exp- ; can't drop any var until this
	keep `uid' `touse' `timevar' `panelvar' `absorb_keepvars' `basevars' `over' `weightvar'

* 5) Expand factor and time-series variables
	local expandedvars
	foreach set of local sets {
		local varlist ``set''
		if ("`varlist'"=="") continue
		local original_`set' `varlist'
		* the -if- prevents creating dummies for categories that have been excluded
		ExpandFactorVariables `varlist' if `touse', setname(`set')
		local `set' "`r(varlist)'"
		local expandedvars `expandedvars' ``set''
	} 

* 6) Drop unused basevars and tsset vars (no longer needed)
	keep `uid' `touse' `absorb_keepvars' `expandedvars' `over' `weightvar'

* 7) Drop all observations with missing values (before creating the FE ids!)
	markout `touse' `expandedvars'
	markout `touse' `expandedvars' `absorb_keepvars'
	qui keep if `touse'
	Assert c(N)>0, rc(2000)
	drop `touse'
	if ("`over'"!="" & `savingcache') qui levelsof `over', local(levels_over)

* 8) Fill Mata structures, create FE identifiers, avge vars and cluster if needed
	reghdfe_absorb, step(precompute) keep(`uid' `expandedvars') depvar("`depvar'") `excludeself'
	
	* Cluster name to show in the regression table
	*local original_clustervar1 `clustervar1'
	mata: st_local("ivars_clustervar1", ivars_clustervar1)
	local original_clustervar1 : subinstr local ivars_clustervar1 " " "#", all
	
	local clustervar1 "`r(clustervar1)'"


	Debug, level(2) msg("(dataset compacted: observations " as result "`RAW_N' -> `c(N)'" as text " ; variables " as result "`RAW_K' -> `c(k)'" as text ")")
	local avgevars = cond("`avge'"=="", "", "__W*__")
	local vars `expandedvars' `avgevars'

	* qui compress `expandedvars' // will recast to -double- later on
	qui de, simple
	local new_mem = string(r(width) * r(N) / 2^20, "%6.2f")
	Debug, level(2) msg("(dataset compacted, c(memory): " as result "`old_mem'" as text "M -> " as result "`new_mem'" as text "M)")

* ?) Check that weights have acceptable values
if ("`weightvar'"!="") {
	local require_integer = ("`weight'"=="fweight")
	local num_type = cond(`require_integer', "integers", "reals")

	local basenote "weight -`weightvar'- can only contain strictly positive `num_type', but"
	qui cou if `weightvar'<0
	Assert (`r(N)'==0), msg("`basenote' `r(N)' negative values were found!")
	qui cou if `weightvar'==0
	Assert (`r(N)'==0), msg("`basenote' `r(N)' zero values were found!")
	qui cou if `weightvar'>=.
	Assert (`r(N)'==0), msg("`basenote' `r(N)' missing values were found!")
	if (`require_integer') {
		qui cou if mod(`weightvar',1)
		Assert (`r(N)'==0), msg("`basenote' `r(N)' non-integer values were found!")
	}
}

* 9) Save the statistics we need before transforming the variables
if (`savingcache') {
	cap drop __FE*__
	cap drop __clustervar*__
}
else {
	* Compute TSS of untransformed depvar
	local tmpweightexp = subinstr("`weightexp'", "[pweight=", "[aweight=", 1)
	qui su `depvar' `tmpweightexp' // BUGBUG: Is this correct?!
	local tss = r(Var)*(r(N)-1)
	assert `tss'<.

* 10) Calculate the degrees of freedom lost due to the FEs
	if ("`group'"!="") {
		tempfile groupdta
		local opt group(`group') groupdta(`groupdta') uid(`uid')
	}
	EstimateDoF, dofmethod(`dofmethod') clustervar1(`clustervar1') `opt'
	// For simplicity, -EstimateDoF- will save its results in locals HERE
	// Locals saved are: Mg, Kg, Mg_exact for all g=1..N_hdfe
	// and M (sum of Mg), kk (sum of Kg), fe_nested_in_cluster

	* Sanity checks
	Assert `kk'<.
	Assert inlist(`fe_nested_in_cluster',0,1)
	forv g=1/`N_hdfe' {
		assert inlist(`M`g'_exact',0,1)
		// 1 or 0 whether M`g' was calculated exactly or not
		assert `M`g''<. & `K`g''<.
		assert `M`g''>=0 & `K`g''>=0
	}

* 11) Drop IDs for the absorbed FEs (except if its the clustervar)
* Useful b/c regr. w/cluster takes a lot of memory
	if ( `fe_nested_in_cluster' | ("`clustervar1'"!="" & "`dofmethod'"=="naive") ) {
		rename `clustervar1' __clustervar1__
		local clustervar1 __clustervar1__
	}
	if ("`vcetype'"=="cluster") {
		conf numeric var `clustervar1'
		local vceoption : subinstr local vceoption "<CLUSTERVAR1>" "`clustervar1'"
	}
	cap drop __FE*__ // FE IDs no longer needed (except if clustervar)
	* cap is required in case there is only one FE which is the clustervar1
}

* 12) Save untransformed data.
*	This allows us to:
*	i) nested ftests for the FEs,
*	ii) to recover the FEs, compute their correlations with xb, check that FE==1

	* We can avoid this if i) nested=check=0 ii) targets={} iii) fast=1
	mata: st_local("any_target_avge", strofreal(any(avge_target :!= "")) ) // saving avge?
	local any_target_hdfe 0 // saving hdfe?
	forv g=1/`N_hdfe' {
		reghdfe_absorb, fe2local(`g')
		if (!`is_bivariate' | `is_mock') local hdfe_cvar`g' `cvars'
		// If it's the intercept part of the bivariate absorbed effect, don't add the cvar!
		local hdfe_target`g' `target'
		if ("`target'"!="") local any_target_hdfe 1
	}

	if (`fast') {
		if (`nested' | `check' | `any_target_hdfe' | `any_target_avge' | "`group'"!="" | `cores'>1) {
			Debug, msg(as text "(option {it:fast} not compatible with other options; disabled)") level(0) // {opt ..} is too boldy
			local fast 0
		}
		else {
			Debug, msg("(option {opt fast} specified; will not save e(sample) or compute correlations)")
		}
	}

	if (!`fast') {
		sort `uid'
		tempfile original_vars
		qui save "`original_vars'"
		if (`cores'>1) local parallel_opt `" filename("`original_vars'") uid(`uid') cores(`cores') "'
		Debug, msg("(untransformed dataset saved)") level(2)
	}

* 13) (optional) Compute R2/RSS to run nested Ftests on the FEs
	* a) Compute R2 of regression without FE, to build the joint FTest for all the FEs
	* b) Also, compute RSS of regressions with less FEs so we can run nested FTests on the FEs
	if ("`model'"=="ols" & !`savingcache') {
		qui _regress `vars' `weightexp', noheader notable
		local r2c = e(r2)

		if (`nested') {
			local rss0 = e(rss)
			local subZs
			forv g=1/`=`N_hdfe'-1' {
				Debug, msg("(computing nested model w/`g' FEs)")
				reghdfe_absorb, step(demean) varlist(`vars') `maximize_options' num_fe(`g') `parallel_opt'
				qui _regress `vars' `weightexp', noheader notable
				local rss`g' = e(rss)
				qui use "`original_vars'", clear // Back to untransformed dataset
			}
		}
	}

	* Get normalized string of the absvars (i.e. turn -> i.turn)
	local original_absvars
	forv g=1/`N_hdfe' {
		reghdfe_absorb, fe2local(`g')
		local original_absvars `original_absvars'  `varlabel'
	}

* 14) Compute residuals for all variables including the AvgEs (overwrites vars!)
	qui ds `vars'
	local NUM_VARS : word count `r(varlist)'
	Debug, msg("(computing residuals for `NUM_VARS' variables)")
	Debug, msg(" - tolerance = `tolerance'")
	Debug, msg(" - max. iter = `maxiterations'")
	if ("`usecache'"=="") {
		reghdfe_absorb, step(demean) varlist(`vars') `maximize_options' `parallel_opt'
	}
	else {
		Debug, msg("(using cache data)")
		drop `vars'
		local handshake_master : char __uid__[handshake]
		char __uid__[handshake]
		// An error in the merge most likely means different # of obs due to missing values in a group but not in other
		// try with if !missing(__uid__) // TODO: Auto-add this by default?
		// TODO: Make this fool-proof when using -over-
		if ("`over'"!="") local using using // This is dangerous
		sort __uid__ // The user may have changed the sort order of the master data
		qui merge 1:1 __uid__ using "`usecache'", keepusing(`vars') assert(match master `using') keep(master match) nolabel sorted
		qui cou if _merge!=3
		if (r(N)>0) {
			Debug, level(0) msg(as error "Warning: the cache has `r(N)' less observations than the master data" _n as text ///
				" - This is possibly because, when created, it included variables that were missing in cases where the current ones are not." _n ///
				" - It may or may not be an error depending on your objective.")
		}
		qui drop if _merge!=3
		drop _merge

		local handshake_using : char __uid__[handshake]
		local tolerance_using : char __uid__[tolerance]
		local maxiterations_using : char __uid__[maxiterations]
		Assert (`handshake_master'==`handshake_using'), msg("using dataset does not have the same __uid__")
		Assert abs(`tolerance'-`tolerance_using')<epsdouble(), msg("using dataset not computed with the same tolerance (`tolerance_using')")
		Assert (`maxiterations'==`maxiterations_using'), msg("using dataset not computed with the same maxiterations (`maxiterations_using')")

		local absvar_master `original_absvars'
		local absvar_using : char __uid__[absvars_key]
		Assert ("`absvar_master'"=="`absvar_using'"), msg("using dataset not created with the same absvars")
		char __uid__[absvars_key]
	}

if (`savingcache') {
	Debug, msg("(saving cache and exiting)")
	char __uid__[absvars_key] `original_absvars'
	sort __uid__
	save "`savecache'", replace
	if ("`levels_over'"!="") ereturn local levels_over = "`levels_over'"
	exit
}

**** PART II - REGRESSION ****

* 1) Cleanup
	ereturn clear

* 2) Debugging - This will allow to run an equivalent -regress- command with -reghdfe, alt-
	
	* If we didn't save (set a target) for all avge, we can't replicate the command
	local alternative_ok 1
	if (`N_avge'>0) {
		mata: st_local("alternative_ok", strofreal(all(avge_target :!= "")) )
		if (`alternative_ok') mata: st_local("avge_targets", strtrim(invtokens(avge_target)) )
	}
	if (`: word count `ivars_clustervar1''>1) {
		local alternative_ok 0
	}

* 3) Regression
	Debug, level(2) msg("(running regresion: `model'.`ivsuite')")
	local avge = cond(`N_avge'>0, "__W*__", "")
	local options
	local option_list ///
		depvar indepvars endogvars instruments avgevars ///
		original_depvar original_indepvars original_endogvars ///
		original_instruments original_absvars avge_targets ///
		estimator vceoption vcetype kk suboptions showraw dofminus first  weightexp
	foreach opt of local option_list {
		if ("``opt''"!="") local options `options' `opt'(``opt'')
	}
	Debug, level(3) msg(_n "call to wrapper:" _n as result "Wrapper_`subcmd', `options'")
	Wrapper_`subcmd', `options' 
	local subpredict = e(predict) // used to recover the FEs

	if ("`weightvar'"!="") {
		qui su `weightvar', mean
		local sumweights = r(sum)
	}

**** PART III - RECOVER FEs AND SAVE RESULTS ****

if (`fast') {
	* Copy pasted from below
	Debug, level(3) msg("(avoiding -use- of temporary dataset")
	tempname b
	matrix `b' = e(b)
	local backup_colnames : colnames `b'
	FixVarnames `backup_colnames'
	local newnames "`r(newnames)'"
	local prettynames "`r(prettynames)'"
	matrix colnames `b' = `newnames'

	clear // can comment out after debugging
}
else {

* 1) Restore untransformed dataset
	qui use "`original_vars'", clear

* 2) Recover the FEs

	* Predict will get (e+d) from the equation y=xb+d+e
	tempvar resid_d
	`subpredict' double `resid_d', resid // Auto-selects the program based on the estimation method
	Debug, level(2) msg("(loaded untransformed variables, predicted residuals)")

	* Absorb the residuals to obtain the FEs (i.e. run a regression on just the resids)
	Debug, level(2) tic(31)
	reghdfe_absorb, step(demean) varlist(`resid_d') `maximize_options' save_fe(1)
	Debug, level(2) toc(31) msg("mata:make_residual on final model took")
	drop `resid_d'

* 3) Compute corr(FE,xb) (do before rescaling by cvar or deleting)
	if ("`model'"=="ols") {
		tempvar xb
		_predict double `xb', xb // -predict- overwrites sreturn, use _predict if needed
		forv g=1/`N_hdfe' { 
			qui corr `xb' __Z`g'__
			local corr`g' = r(rho)
		}
		drop `xb'
	}

* 4) Replace tempnames in the coefs table
	* (e.g. __00001 -> L.somevar)
	* (this needs to be AFTER predict but before deleting FEs and AvgEs)
	tempname b
	matrix `b' = e(b)
	local backup_colnames : colnames `b'
	FixVarnames `backup_colnames'
	local newnames "`r(newnames)'"
	local prettynames "`r(prettynames)'"
	matrix colnames `b' = `newnames'

* 5) Save FEs w/proper name, format
	reghdfe_absorb, step(save) original_depvar(`original_depvar')
	local keepvars `r(keepvars)'
	if ("`keepvars'"!="") format `fe_format' `keepvars'

* 6) Save AvgEs
	forv g=1/`N_avge' {
		local var __W`g'__
		local target : char `var'[target]
		if ("`target'"!="") {
			rename `var' `target'
			local avge_target`g' `target' // Used by -predict-
			local keepvars `keepvars' `target'
		}
	}

	if ("`keepvars'"!="") format `fe_format' `keepvars' // The format of depvar, saved by -Parse-

* 7) Save dataset with FEs and e(sample)
	keep `uid' `keepvars'
	tempfile output
	qui save "`output'"
} // fast

* 8) Restore original dataset and merge
	restore // Restore user-provided dataset
	if (!`fast') {
		// `saved_group' was created by EstimateDoF.ado
		if (!`saved_group')  local groupdta
		SafeMerge, uid(`uid') file("`output'") groupdta("`groupdta'")
		*cap tsset, noquery // we changed -sortby- when we merged (even if we didn't really resort)
	}

**** PART IV - ERETURN OUTPUT ****

	if (`c(version)'>=12) local hidden hidden // ereturn hidden requires v12+

* Ereturns common to all commands
	ereturn local cmd = "reghdfe"
	ereturn local subcmd = "`subcmd'"
	ereturn local cmdline `"`cmdline'"'
	ereturn scalar alternative_ok = `alternative_ok'
	if ("`e(model)'"!="" & "`e(model)'"!="`model'") di as error "`e(model) was <`e(model)'>" // ?
	ereturn local model = "`model'"
	ereturn local dofmethod = "`dofmethod'"
	ereturn local title = "HDFE " + e(title)
	ereturn local subtitle =  "Absorbing `N_hdfe' HDFE " + plural(`N_hdfe', "indicator")
	ereturn local predict = "reghdfe_p"
	ereturn local estat_cmd = "reghdfe_estat"
	ereturn local footnote = "reghdfe_footnote"
	ereturn local absvars = "`original_absvars'"
	ereturn `hidden' local diopts = "`diopts'"

	if ("`e(clustvar)'"!="") ereturn local clustvar "`original_clustervar1'"

	* Besides each cmd's naming style (e.g. exogr, exexog, etc.) keep one common one
	foreach cat in depvar indepvars endogvars instruments {
		local vars ``cat''
		if ("`vars'"=="") continue
		ereturn local `cat' "`original_`cat''"
	}
	ereturn local avgevars "`avge'" // bugbug?

	ereturn `hidden' local subpredict = "`subpredict'"
	ereturn `hidden' local prettynames "`prettynames'"
	forv g=1/`N_avge' {
		ereturn `hidden' local avge_target`g' "`avge_target`g''" // Used by -predict-
	}

	* Stata uses e(vcetype) for the SE column headers
	* In the default option, leave it empty.
	* In the cluster and robust options, set it as "Robust"
	ereturn local vcetype = proper("`vcetype'") //
	if (e(vcetype)=="Cluster") ereturn local vcetype = "Robust"
	if (e(vcetype)=="Unadjusted") ereturn local vcetype
	ereturn local vce = "`vcetype'" // Not a mistake! Saves unadjusted,robust,cluster,etc. into e(vce)
	Assert inlist("`e(vcetype)'", "", "Robust", "Jackknife", "Bootstrap")
	Assert inlist(e(vce), "unadjusted", "robust", "cluster", "jacknife", "bootstrap") // Redundant, already checked in Parse.ado
	
	* Clear results that are wrong
	ereturn local ll
	ereturn local ll_0

	ereturn scalar N_hdfe = `N_hdfe'
	if ("`N_avge'"!="") ereturn scalar N_avge = `N_avge'

* Absorbed-specific returns
	ereturn scalar mobility = `M'
	ereturn scalar df_a = `kk'
	if ("`vcetype'"=="cluster") ereturn scalar fe_nested_in_cluster = `fe_nested_in_cluster'
	forv g=1/`N_hdfe' {
		ereturn scalar M`g' = `M`g''
		ereturn scalar K`g' = `K`g''
		ereturn `hidden' scalar M`g'_exact = `M`g'_exact' // 1 or 0 whether M`g' was calculated exactly or not
		ereturn `hidden' local corr`g' = "`corr`g''" //  cond("`corr`g''"=="", ., "`corr`g''")
		ereturn `hidden' local hdfe_target`g' = "`hdfe_target`g''"
		ereturn `hidden' local hdfe_cvar`g' = "`hdfe_cvar`g''"
	}

	Assert e(df_r)<. , msg("e(df_r) is missing")
	ereturn scalar tss = `tss'
	ereturn scalar mss = e(tss) - e(rss)
	ereturn scalar r2 = 1 - e(rss) / `tss'

	ereturn scalar r2_a = 1 - (e(rss)/e(df_r)) / (`tss' / (e(N)-1) ) // After fixing e(df_r)
	ereturn scalar rmse = sqrt( e(rss) / e(df_r) )

	if ("`weightvar'"!="") ereturn scalar sumweights = `sumweights'

	if ("`model'"=="ols" & "`vcetype'"=="unadjusted") {
		ereturn scalar F_absorb = (e(r2)-`r2c') / (1-e(r2)) * e(df_r) / `kk'
		if (`nested') {
			local rss`N_hdfe' = e(rss)
			local temp_dof = e(N) - 1 - e(df_m) // What if there are absorbed collinear with the other RHS vars?
			local j 0
			ereturn `hidden' scalar rss0 = `rss0'
			forv g=1/`N_hdfe' {
				local temp_dof = `temp_dof' - e(K`g') + e(M`g')
				*di in red "g=`g' RSS=`rss`g'' and was `rss`j''.  dof=`temp_dof'"
				ereturn `hidden' scalar rss`g' = `rss`g''
				ereturn `hidden' scalar df_a`g' = e(K`g') - e(M`g')
				ereturn scalar F_absorb`g' = (`rss`j''-`rss`g'') / `rss`g'' * `temp_dof' / e(df_a`g')
				ereturn `hidden' scalar df_r`g' = `temp_dof'
				local j `g'
			}   
		}
	}

	if (e(N_clust)<.) ereturn scalar df_r = e(N_clust) - 1
	// There is a big assumption hre, that the number of other parameters does not increase asymptotically
	// BUGBUG: We should allow the option to indicate what parameters do increase asympt.
	// BUGBUG; xtreg does this: est scalar df_r = min(`df_r':=N-1-K, `df_cl') why was that?

	if ("`savefirst'"!="") ereturn `hidden' scalar savefirst = `savefirst'

* Show table and clean up
	ereturn repost b=`b', rename // why here???
	Replay
	reghdfe_absorb, step(stop)

end

* The idea of this program is to keep the sort order when doing the merges
program define SafeMerge, eclass sortpreserve
syntax, uid(varname numeric) file(string) [groupdta(string)]
	* Merging gives us e(sample) and the FEs / AvgEs
	tempvar merge
	merge 1:1 `uid' using "`file'", assert(master match) nolabel nonotes noreport gen(`merge')
	
	* Add e(sample) from _merge
	tempvar sample
	gen byte `sample' = (`merge'==3)
	la var `sample' "[HDFE Sample]"
	ereturn repost , esample(`sample')
	drop `merge'

	* Add mobility group
	if ("`groupdta'"!="") merge 1:1 `uid' using "`groupdta'", assert(master match) nogen nolabel nonotes noreport sorted
end


	
// -------------------------------------------------------------
// Parsing and basic sanity checks for REGHDFE.ado
// -------------------------------------------------------------
// depvar: dependent variable
// indepvars: included exogenous regressors
// endogvars: included endogenous regressors
// instruments: excluded exogenous regressors
program define Parse

* Remove extra spacing from cmdline (just for aesthetics, run before syntax)
	cap syntax anything(name=indepvars) [if] [in] [fweight aweight pweight/] , SAVEcache(string) [*]
	local savingcache = (`=_rc'==0)

if (`savingcache') {

	* Disable these options
	local fast
	local nested

	syntax anything(name=indepvars) [if] [in] [fweight aweight pweight/] , ///
		Absorb(string) SAVEcache(string) ///
		[Verbose(integer 0) CHECK TOLerance(real 1e-7) MAXITerations(integer 1000) noACCELerate ///
		bad_loop_threshold(integer 1) stuck_threshold(real 5e-3) pause_length(integer 20) ///
		accel_freq(integer 3) accel_start(integer 6) /// Advanced optimization options
		CORES(integer 1) OVER(varname numeric)]

	cap conf file "`savecache'.dta"
	if (`=_rc'!=0) {
		cap conf new file "`savecache'.dta"
		Assert (`=_rc'==0), msg("reghdfe will not be able to save `savecache'.dta")
	}

}
else {
	mata: st_local("cmdline", stritrim(`"reghdfe `0'"') )
	ereturn clear // Clear previous results and drops e(sample)
	syntax anything(id="varlist" name=0 equalok) [if] [in] ///
		[fweight aweight pweight/] , ///
		Absorb(string) ///
		[GROUP(name) VCE(string) Verbose(integer 0) CHECK NESTED FAST] ///
		[avge(string) EXCLUDESELF] ///
		[TOLerance(real 1e-7) MAXITerations(integer 1000) noACCELerate] /// See reghdfe_absorb.Annihilate
		[noTRACK] /// Not used here but in -Track-
		[IVsuite(string) ESTimator(string) SAVEFIRST FIRST SHOWRAW dofminus(string)] ///
		[dofmethod(string)] ///
		[SMALL noCONstant Hascons TSSCONS] /// ignored options
		[gmm2s liml kiefer cue] ///
		[SUBOPTions(string)] /// Options to be passed to the estimation command (e.g . to regress)
		[bad_loop_threshold(integer 1) stuck_threshold(real 5e-3) pause_length(integer 20) ///
		accel_freq(integer 3) accel_start(integer 6)] /// Advanced optimization options
		[CORES(integer 1)] [USEcache(string)] [OVER(varname numeric)] ///
		[*] // For display options
}

* Weight
* We'll have -weight- (fweight|aweight|pweight), -weightvar-, -exp-, and -weightexp-
	if ("`weight'"!="") {
		local weightvar `exp'
		conf var `weightvar' // just allow simple weights
		local weightexp [`weight'=`weightvar']
		local backupweight `weight'
	}

* Cache options
	if ("`usecache'"!="") {
		conf file "`usecache'.dta"
		conf var __uid__
		Assert ("`avge'"==""), msg("option -avge- not allowed with -usecache-")
		Assert ("`avge'"==""), msg("option -nested- not allowed with -usecache-")
	}

* Save locals that will be overwritten by later calls to -syntax-
	local ifopt `if'
	local inopt `in'

* Coef Table Options
if (!`savingcache') {
	_get_diopts diopts options, `options'
	Assert `"`options'"'=="", msg(`"invalid options: `options'"')
	if ("`constant'`hascons'`tsscons'"!="") di in ye "(option `constant'`hascons'`tsscons' ignored)"
}

* Over
	if ("`over'"!="") {
		unab over : `over', max(1)
		Assert ("`usecache'"!="" | "`savecache'"!=""), msg("-over- needs to be used together with either -usecache- or -savecache-")
	}

* Verbose
	assert inlist(`verbose', 0, 1, 2, 3, 4) // 3 and 4 are developer options
	mata: VERBOSE = `verbose' // Ugly hack to avoid using a -global-

* Model settings
if (!`savingcache') {

	// local model = cond(strpos(`"`0'"', " ("), "iv", "ols") // Fails with long strs in stata 12<
	local model ols
	foreach _ of local 0 {
		if (substr(`"`_'"', 1, 1)=="(") {
			local model iv
			continue, break
		}
	}
	

	if ("`model'"=="iv") {

		// get part before parentheses
		local wrongparens 1
		while (`wrongparens') {
			gettoken tmp 0 : 0 ,p("(")
			local left `left'`tmp'
			* Avoid matching the parens of e.g. L(-1/2) and L.(var1 var2)
			* Using Mata to avoid regexm() and trim() space limitations
			mata: st_local("tmp1", subinstr("`0'", " ", "") ) // wrong parens if ( and then a number
			mata: st_local("tmp2", substr(strtrim("`left'"), -1) ) // wrong parens if dot
			local wrongparens = regexm("`tmp1'", "^\([0-9-]") | ("`tmp2'"==".")
			if (`wrongparens') {
				gettoken tmp 0 : 0 ,p(")")
				local left `left'`tmp'
			}
		}

		// get part in parentheses
		gettoken right 0 : 0 ,bind match(parens)
		Assert trim(`"`0'"')=="" , msg("error: remaining argument: `0'")

		// now parse part in parentheses
		gettoken endogvars instruments : right ,p("=")
		gettoken equalsign instruments : instruments ,p("=")

		Assert "`endogvars'"!="", msg("iv: endogvars required")
		local 0 `endogvars'
		syntax varlist(fv ts numeric)

		Assert "`instruments'"!="", msg("iv: instruments required")
		local 0 `instruments'
		syntax varlist(fv ts numeric)
		
		local 0 `left' // So OLS part can handle it
		Assert "`endogvars'`instruments'"!=""
		
		if ("`ivsuite'"=="") local ivsuite ivreg2
		Assert inlist("`ivsuite'","ivreg2","ivregress") , msg("error: wrong IV routine (`ivsuite'), valid options are -ivreg2- and -ivregress-")
		cap findfile `ivsuite'.ado
		Assert !_rc , msg("error: -`ivsuite'- not installed, please run ssc install `ivsuite' or change the option -ivsuite-")
	}

* OLS varlist
	syntax varlist(fv ts numeric)
	gettoken depvar indepvars : 0
	_fv_check_depvar `depvar'

* Extract format of depvar so we can format FEs like this
	fvrevar `depvar', list
	local fe_format : format `r(varlist)' // The format of the FEs and AvgEs that will be saved

* Variables shouldn't be repeated
* This is not perfect (e.g. doesn't deal with "x1-x10") but still helpful
	local allvars `depvar' `indepvars' `endogvars' `instruments'
	local dupvars : list dups allvars
	Assert "`dupvars'"=="", msg("error: there are repeated variables: <`dupvars'>")

	Debug, msg(_n " {title:REGHDFE} Verbose level = `verbose'")
	*Debug, msg("{hline 64}")

* VCE
	gettoken vcetype vcerest : vce, parse(" ,")
	if ("`vcetype'"=="") {
		local vcetype unadjusted
		if ("`ivsuite'"=="ivregress" & "`estimator'"=="gmm") local vcetype robust // Default for ivregress gmm
		if ("`backupweight'"=="pweight") local vcetype robust // pweight requires robust
	}

	* Abbreviations:
	if (substr("`vcetype'",1,2)=="un") local vcetype unadjusted
	if (substr("`vcetype'",1,1)=="r") local vcetype robust
	if (substr("`vcetype'",1,2)=="cl") local vcetype cluster
	if (substr("`vcetype'",1,4)=="boot") local vcetype bootstrap
	if (substr("`vcetype'",1,4)=="jack") local vcetype jackknife
	if ("`vcetype'"=="conventional") local vcetype unadjusted // Conventional is the name given in e.g. xtreg
	assert inlist("`vcetype'", "unadjusted", "robust", "cluster", "bootstrap", "jackknife")

	* Cluster vars
	local num_clusters 0
	if ("`vcetype'"=="cluster") {
		local num_clusters = `: word count `vcerest''
		assert inlist(`num_clusters', 1, 2)
		gettoken clustervar1 clustervar2 : vcerest
		cap unab clustervar1 : `clustervar1'
		cap unab clustervar2 : `clustervar2'
		local vcerest
	}

	//local vceoption = trim("`vcetype' `clustervar1' `clustervar2' `vcerest'")
	local temp_clustervar1 = cond("`clustervar1'"=="","","<CLUSTERVAR1>")
	local vceoption = trim("`vcetype' `temp_clustervar1'`vcerest'")
	local vceoption vce(`vceoption')
	
	if ("`model'"=="ols") {
		local ok = inlist("`vcetype'", "unadjusted", "robust", "cluster") & (`num_clusters'<=1)
	}
	else if ("`ivsuite'"=="ivregress") {
		local ok = inlist("`vcetype'", "unadjusted", "robust", "cluster", "bootstrap", "jackknife", "hac") & (`num_clusters'<=1)
	}
	else if ("`ivsuite'"=="ivreg2") {
		local ok = inlist("`vcetype'", "unadjusted", "robust", "cluster", "bw") & (`num_clusters'<=2)
	}
	Assert `ok', msg("invalid vce type or number of clusters")
}

* Optimization
	if (`maxiterations'==0) local maxiterations 1e7
	Assert (`maxiterations'>0)
	local accelerate = cond("`accelerate'"!="", 0, 1) // 1=Yes
	local check = cond("`check'"!="", 1, 0) // 1=Yes
	local fast = cond("`fast'"!="", 1, 0) // 1=Yes
	local tolerance = strofreal(`tolerance', "%9.1e") // Purely esthetic
	Assert `cores'<=32 & `cores'>0 , msg("At most 32 cores supported")
	if (`cores'>1) {
		cap findfile parallel.ado
		Assert !_rc , msg("error: -parallel- not installed, please run {cmd:ssc install parallel}")
	}
	local opt_list tolerance maxiterations check accelerate ///
		bad_loop_threshold stuck_threshold pause_length accel_freq accel_start
	foreach opt of local opt_list {
		if ("``opt''"!="") local maximize_options `maximize_options' `opt'(``opt'')
	}

* IV options
if (!`savingcache') {
	if ("`model'"=="iv" & "`ivsuite'"=="ivregress") {
		if ("`estimator'"=="") local estimator 2sls
		Assert inlist("`estimator'","2sls","gmm"), msg("liml estimator not allowed")
		if ("`estimator'"!="2sls") di as error "WARNING! GMM not fully implemented, robust option gives wrong output"
	}
	else if ("`model'"=="iv" & "`ivsuite'"=="ivreg2") {
		if ("`estimator'"=="") local estimator 2sls
    	Assert inlist("`estimator'","2sls","gmm2s") , msg("Estimator is `estimator', instead of empty (2sls/ols) or gmm2s")
    	* di in ye "WARNING: IV estimates not fully tested; methods like GMM2S will not work correctly"
		if ("`showraw'"!="") local showraw 1

		Assert inlist("`dofminus'","","large","small"), msg("option -dofminus- is either -small- or -large-")
		local dofminus = cond("`dofminus'"=="large", "dofminus", "sdofminus") // Default uses sdofminus
	}
	else {
		local estimator
	}
	if ("`small'"!="") di in ye "(note: reghdfe will always use the option -small-, no need to specify it)"

	Assert ("`gmm2s'`liml'`kiefer'`cue'"==""), msg("Please use estimator() option instead of gmm2s/liml/etc")
	
	if ("`model'"=="iv") {
		local savefirst = ("`savefirst'"!="")
		local first = ("`first'"!="")
		if (`savefirst') Assert `first', msg("Option -savefirst- requires -first-")
	}

* DoF Estimation
* For more than 2 FEs, no exact soln
* Can estimate M3 and so on using
* i) bounds using FE1 and FE2, ii) assuming M3=1 (fast), iii) bootstrap
	if ("`dofmethod'"=="") local dofmethod bounds
	assert inlist("`dofmethod'", "naive", "simple", "bounds", "bootstrap")

* Mobility groups
	if ("`group'"!="") conf new var `group'
}

* tsset variables, if any
	cap conf var `_dta[_TStvar]'
	if (!_rc) local timevar `_dta[_TStvar]'
	cap conf var `_dta[_TSpanel]'
	if (!_rc) local panelvar `_dta[_TSpanel]'

* Varnames underlying tsvars and fvvars (e.g. i.foo L(1/3).bar -> foo bar)
	foreach vars in depvar indepvars endogvars instruments {
		if ("``vars''"!="") {
			fvrevar ``vars'' , list
			local basevars `basevars' `r(varlist)'
		}
	}

if (!`savingcache') {
* Nested
	local nested = cond("`nested'"!="", 1, 0) // 1=Yes
	if (`nested' & !("`model'"=="ols" & "`vcetype'"=="unadjusted") ) {
		Debug, level(0) msg("(option nested ignored, only works with OLS and conventional/unadjusted VCE)") color("error")
	}

* How can we do the same regression from a standard stata command?
* (useful for benchmarking and testing correctness of results)
	local subcmd = cond("`model'"=="ols" ,"regress", "`ivsuite'")

* _fv_check_depvar overwrites the local -weight-
	local weight `backupweight'
	Assert inlist( ("`weight'"!="") + ("`weightvar'"!="") + ("`weightexp'"!="") , 0 , 3 ) , msg("not all 3 weight locals are set")

* Return values
	local names cmdline diopts model ///
		ivsuite estimator showraw dofminus ///
		depvar indepvars endogvars instruments savefirst first ///
		vceoption vcetype num_clusters clustervar1 clustervar2 vcerest ///
		if in group dofmethod check fast nested fe_format ///
		tolerance maxiterations accelerate maximize_options ///
		subcmd suboptions ///
		absorb avge excludeself ///
		timevar panelvar basevars ///
		weight weightvar exp weightexp /// type of weight (fw,aw,pw), weight var., and full expr. ([fw=n])
		cores savingcache usecache over
}

if (`savingcache') {
	local names maximize_options cores if in timevar panelvar indepvars basevars ///
		absorb savecache savingcache fast nested check over ///
		weight weightvar exp weightexp /// type of weight (fw,aw), weight var., and full expr. ([fw=n])
		tolerance maxiterations // Here just used for -verbose- and cache handshake purposes
}

	local if `ifopt'
	local in `inopt'

	Debug, level(3) newline
	Debug, level(3) msg("Parsed options:")
	foreach name of local names {
		if (`"``name''"'!="") Debug, level(3) msg("  `name' = " as result `"``name''"')
		c_local `name' `"``name''"' // Inject values into caller (reghdfe.ado)
	}
	// Debug, level(3) newline
end

	
//------------------------------------------------------------------------------
// Expand Factor Variables, interactions, and time-series vars
//------------------------------------------------------------------------------
// This basically wraps -fvrevar-, adds labels, and drops omitted/base
program define ExpandFactorVariables, rclass
syntax varlist(min=1 numeric fv ts) [if] [,setname(string)] [CACHE]

	local expanded_msg `"" - variable expansion for `setname': " as result "`varlist'" as text " ->""'

	* It's (usually) a waste to add base and omitted categories
	* EG: if we use i.foreign#i.rep78 , several categories will be redundant, seen as e.g. "0b.foreign" in -char list-
	* We'll also exclude base categories that don't have the "bn" option (to have no base)

	* Loop for each var and then expand them into i.var -> 1.var.. and loop
	* Why two loops? B/c I want to save each var expansion to allow for a cache

	if ("`cache'"!="") mata: varlist_cache = asarray_create()

	local newvarlist
	* I can't do a simple foreach!
	* Because a factor expression could be io(3 4).rep78
	* and foreach would split the parens in half
	while (1) {
	gettoken fvvar varlist : varlist, bind
	if ("`fvvar'"=="") continue, break

		fvrevar `fvvar' `if' // , stub(__V__) // stub doesn't work in 11.2
		local contents

		foreach var of varlist `r(varlist)' {
			
			* Get readable varname
			local fvchar : char `var'[fvrevar]
			local tschar : char `var'[tsrevar]
			local name `fvchar'`tschar'
			local color input
			if ("`name'"=="") {
				local name `var'
				local color result
			}
			char `var'[name] `name'
			la var `var' "[Tempvar] `name'"

			* See if the factor can be dropped safely
			if (substr("`var'", 1, 2)=="__") {
				local color result
				local parts : subinstr local fvchar "#" " ", all
				foreach part of local parts {
					* "^[0-9]+b\." -> "b.*\."
					if regexm("`part'", "b.*\.") | regexm("`part'", "o.*\.") {
						local color error	
						drop `var'
						continue, break
					}
				}


				* Need to rename it, or else it gets dropped since its a tempvar
				if ("`color'"!="error") {
					local newvarbase : subinstr local name "." "__", all // pray that no variable has three _
					local newvarbase : subinstr local newvarbase "#" "_X_", all // idem
					local newvarbase : permname __`newvarbase', length(30)
					local i 0
					while (1) {
						local newvar "`newvarbase'`++i'"
						Assert `i'<1000, msg("Couldn't create tempvar for `var' (`name')")
						cap conf new var `newvar', exact
						if _rc==0 {
							continue, break
						}
					}
					rename `var' `newvar'
					local var `newvar'
				}
			}

			* Save contents of the expansion for optional -cache-			
			if ("`color'"!="error") {
				local contents `contents' `var'
			}
			
			* Set debug message
			local expanded_msg `"`expanded_msg' as `color' " `name'" as text " (`var')""'
		}

		if ("`cache'"!="") mata: asarray(varlist_cache, "`fvvar'", "`contents'")
		Assert "`contents'"!="", msg("error: variable -`fvvar'- in varlist -`varlist'- in category -`setname'- is  empty after factor expansion")
		local newvarlist `newvarlist' `contents'
	}

	* Yellow=Already existed, White=Created, Red=NotCreated (omitted or base)
	Debug, level(3) msg(`expanded_msg')
	return local varlist "`newvarlist'"
end

	
// -------------------------------------------------------------
// Calculate DoF lost from the FEs
// -------------------------------------------------------------
program define EstimateDoF
syntax, dofmethod(string) [clustervar1(string) group(name) uid(varname) groupdta(string)]
	
	Assert inlist("`dofmethod'", "bounds", "simple", "naive", "bootstrap")
	Assert "`dofmethod'"!="bootstrap" , msg("DoF Bootstrap: not yet implemented!")
	Debug, level(1) msg("(calculating degrees of freedom lost by the FEs)")
	Debug, level(2) msg(" - dofmethod: `dofmethod'")
	mata: st_local("G", strofreal(G))

* Start conservatively, assuming M`g'=1 (or =0 if interacted with cont. var)
	local g_list // List of FEs where M is still unknown
	forv g=1/`G' {
		reghdfe_absorb, fe2local(`g')
		// ivars cvars target varname varlabel is_interaction is_cont_interaction is_bivariate is_mock levels
		local M`g' = ("`cvars'"=="") | (`is_bivariate' & !`is_mock')
		* We know M1=1 and Mg=0 if g is has a continuous interaction. The others we don't
		if (`M`g''==1 & `levels'>1) local g_list `g_list' `g'
		// If the FE is just an intercept, it's redundant (useful with -over-)
	}


* (ADDENDUM) Look for nested within cluster for the categorical components of a FE with cont. interaction
	if ("`clustervar1'"!="" & "`dofmethod'"!="naive") {
		forv g=1/`G' {
			reghdfe_absorb, fe2local(`g')
			if (`M`g''==0) {
				local gg = `g' - `is_mock'
				cap _xtreg_chk_cl2 `clustervar1' __FE`gg'__
				assert inlist(_rc, 0, 498)
				if (!_rc) {
					local prettyname : char __FE`gg'__[name]
					Debug, msg("(Absorbed variable " as result "`prettyname'" ///
							as text " is nested within cluster " ///
							as result "`clustervar1'" as text ", adjusting DoF)")
					local M`g' `levels'
				}
			}
		}
	}

* Is the cluster var one of the FEs, or contains one of the FEs?
* Since the number of clusters is the effective "number of observations", 
* we shouldn't penalize the DoF for estimating means within a cluster ("within an obs")
	local g_cluster 0
	if ("`clustervar1'"!="" & "`dofmethod'"!="naive") {
	
* 1) See if it's exactly the same variable
		local regex = regexm("`clustervar1'","^__FE([0-9]+)__$")
		if (`regex') {
			local g_cluster = `=regexs(1)'
			local prettyname : char __FE`g_cluster'__[name]
			Debug, msg("(cluster variable " as result "`prettyname'" ///
				as text " is also an absorbed variable, adjusting DoF)")
		}
* 2) If that failed, see if one of the panels is nested within cluster
* i.e., if whenever two obs are in the same group (in the FE), 
* they are also in the same cluster
		else {
			foreach g of local g_list {
				reghdfe_absorb, fe2local(`g')
				cap _xtreg_chk_cl2 `clustervar1' __FE`g'__
				assert inlist(_rc, 0, 498)

				if (!_rc) {
					local g_cluster = `g'
					local prettyname : char __FE`g'__[name]
					Debug, msg("(Absorbed variable " as result "`prettyname'" ///
							as text " is nested within cluster " ///
							as result "`clustervar1'" as text ", adjusting DoF)")
					continue, break
				}
			}
		}

		if (`g_cluster'>0) {
			reghdfe_absorb, fe2local(`g_cluster')
			local M`g_cluster' = `levels'
			local g_list : list g_list - g_cluster
		}
	} // clustervar1

* Compute connected groups for the remaining FEs (except those with cont interactions)

	if ("`group'"!="") local group_option ", gen(`group')"
	local length : list sizeof g_list
	tokenize `g_list' // Saves g1 g2 .. into `1' `2' etc

	if (`length'==0) {
		// pass
	}
	else if (`length'==1) {
		local M`1' = 1
		local g_list // Only one FE besides cont. interaction and clustervar
	}
	else if (`length'==2) {
		Debug, level(1) msg(" - exact DoF computation using connected groups")
		ConnectedGroups __FE`1'__ __FE`2'__ `group_option'
		local M`1' = 1
		local M`2' = r(groups)
		local g_list
		local label "`: char __FE`1'__[name]' and `: char __FE`2'__[name]'"
	}
	else {
		ConnectedGroups __FE`1'__ __FE`2'__ `group_option'
		local M`1' = 1
		local M`2' = r(groups)
		local g_list : list g_list - 1
		local g_list : list g_list - 2

		* Get a conservative bound for M3 M4, etc:
		* Calculate the mobility group wrt M1 and M2 and use the max of both
		* Of course, this excludes effects in FE3 that are collinear with FE1 and FE2 but not FE1 or FE2 separately (and also excludes searching for FE4 vs FE3, etc)

		if ("`dofmethod'"=="bounds") {
			local msg_done = 0
			* Iterate over the remaining elements
			foreach g of local g_list {
				assert `1'!=`g' & `2'!=`g'
				ConnectedGroups __FE`1'__ __FE`g'__
				local candidate1 = r(groups)
				ConnectedGroups __FE`2'__ __FE`g'__
				local candidate2 = r(groups)
				local M`g' = max(`candidate1', `candidate2')

				reghdfe_absorb, fe2local(`g')
				Assert (`M`g''<=`levels'), msg("Mg should be at most Kg: `M`g''<=`levels'")
				
				* If Mg==Kg , the bound is exact and the FE is redundant!
				if (`M`g''==`levels') local g_list : list g_list - g
			}
			if ("`g_list'"!="") di in ye "(note: conservative DoF estimates; they do not account for all possible collinearities between FEs; dofmethod=`dofmethod')"
		}
		else if inlist("`dofmethod'","simple","naive") {
			di in ye "(note: DoF estimates ignore most possible collinearities between FEs, as dofmethod=`dofmethod')"
		}
		else {
			error 999
		}
	}

	local SumM 0
	local SumK 0
	forv g=1/`G' {
		reghdfe_absorb, fe2local(`g')
		// ivars cvars target varname varlabel is_interaction is_cont_interaction is_bivariate is_mock levels
		assert !missing(`M`g'') & !missing(`levels')
		local SumM = `SumM' + `M`g''
		local SumK = `SumK' + `levels'
		local is_exact = !`: list g in g_list'

		c_local M`g' `M`g''
		c_local K`g' `levels'
		c_local M`g'_exact `is_exact'
		Debug, level(2) msg(" - parameters of FE`g': K=`levels' M=`M`g'' is_exact=`is_exact'")
	}
	local NetSumK = `SumK' - `SumM'
	Debug, level(2) msg(" - DoF loss due to FEs: Sum(Kg)=`SumK', M:Sum(Mg)=`SumM' --> KK:=SumK-SumM=`NetSumK'")

* Save mobility group if needed
	c_local saved_group = 0
	if ("`group'"!="" & `length'>=2) {
		tempfile backup
		qui save "`backup'"
		
		keep `uid' `group'
		sort `uid'
		la var `group' "Mobility group between `label'"
		qui save "`groupdta'" // A tempfile from the caller program
		Debug, level(2) msg(" - mobility group saved")
		qui use "`backup'", clear
		cap erase "`backup'"
		c_local saved_group = 1
	}

	c_local M `SumM'
	c_local kk `NetSumK'
	local fe_nested_in_cluster = (`g_cluster'>0)
	c_local fe_nested_in_cluster `fe_nested_in_cluster'
end

		
// -------------------------------------------------------------
// Faster alternative to -makegps-, but with some limitations
// -------------------------------------------------------------
* To avoid backuping the data, use option -clear-
* For simplicity, disallow -if- and -in- options
program ConnectedGroups, rclass
syntax varlist(min=2 max=2) [, GENerate(name) CLEAR]

* To avoid backuping the data, use option -clear-
* For simplicity, disallow -if- and -in- options

    if ("`generate'"!="") conf new var `generate'
    gettoken id1 id2 : varlist
    Debug, level(2) msg(" - computing connected groups between `id1' and`id2'")
    tempvar group copy

    tempfile backup
    if ("`clear'"=="") qui save "`backup'"
    keep `varlist'
    qui bys `varlist': keep if _n==1


    clonevar `group' = `id1'
    clonevar `copy' = `group'
    capture error 100 // We want an error
    while _rc {
        qui bys `id2' (`group'): replace `group' = `group'[1]
        qui bys `id1' (`group'): replace `group' = `group'[1]
        capture assert `copy'==`group'
        qui replace `copy' = `group'
    }

    assert !missing(`group')
    qui bys `group': replace `group' = (_n==1)
    qui replace `group' = sum(`group')
    
    su `group', mean
    local num_groups = r(max)
    
    if ("`generate'"!="") rename `group' `generate'
    
    if ("`clear'"=="") {
        if ("`generate'"!="") {
            tempfile groups
            qui compress
            la var `generate' "Mobility group for (`varlist')"
            qui save "`groups'"
            qui use "`backup'", clear
            qui merge m:1 `id1' `id2' using "`groups'" , assert(match) nogen
        }
        else {
            qui use "`backup'", clear
        }
    }
    
    return scalar groups=`num_groups'
end

	
//------------------------------------------------------------------------------
// Name tempvars into e.g. L.x i1.y i2.y AvgE:z , etc.
//------------------------------------------------------------------------------
program define FixVarnames, rclass
local vars `0'

	foreach var of local vars {
		local newname
		local pretyname

		* -var- can be <o.__W1__>
		if ("`var'"=="_cons") {
			local newname `var'
			local prettyname `var'
		}
		else {
			fvrevar `var', list
			local basevar "`r(varlist)'"
			local label : var label `basevar'
			local is_avge = regexm("`basevar'", "^__W[0-9]+__$")
			local is_temp = substr("`basevar'",1,2)=="__"
			local is_omitted = strpos("`var'", "o.")
			local prefix = cond(`is_omitted'>0, "o.", "")
			local name : char `basevar'[name]

			if (`is_avge') {
				local avge_str : char `basevar'[avge_equation]
				local name : char `basevar'[name]
				local prettyname `avge_str':`prefix'`name'

				local newname : char `basevar'[target]
				if ("`newname'"=="") local newname `var'
			}
			else if (`is_temp' & "`name'"!="") {
				local newname `prefix'`name' // BUGBUG
				local prettyname `newname'
			}
			else {
				local newname `var'
				local prettyname `newname'
			}
		}
		
		* di in red " var=<`var'> --> new=<`newname'> pretty=<`prettyname'>"
		Assert ("`newname'"!="" & "`prettyname'"!=""), ///
			msg("var=<`var'> --> new=<`newname'> pretty=<`prettyname'>")
		local newnames `newnames' `newname'
		local prettynames `prettynames' `prettyname'
	}

	local A : word count `vars'
	local B : word count `newnames'
	local C : word count `prettynames'
	Assert `A'==`B', msg("`A' vars but `B' newnames")
	Assert `A'==`C', msg("`A' vars but `C' newnames")
	
	***di as error "newnames=`newnames'"
	***di as error "prettynames=`prettynames'"

	return local newnames "`newnames'"
	return local prettynames "`prettynames'"
end

	
/* Notes:
- For -cluster- I need to run two regressions as in xtreg; the first one is to get the df_m

- El FTEST es distinto entre areg/reg y xtreg porque xtreg hace un ajuste extra
para pasar de areg a xtreg, multiplicar el F por Q^2
Donde Q =  (e(N) - e(rank)) / (e(N) - e(rank) - e(df_a))
Es decir, en vez de dividir entre N-K-KK, me basta con dividir entre N-K
Asi que me bastaria usar -test- despues de correr la regresion y deberia salir igual que el FTEST ajustado del areg!!!
(tambien igual al del xreg pero eso es mas limitante , aunq igual probar para 1 HDFE creando t=_n a nivel del ID1)
*/
*/
program define Wrapper_regress, eclass
	syntax , depvar(varname) [indepvars(varlist) avgevars(varlist)] ///
		original_absvars(string) original_depvar(string) [original_indepvars(string) avge_targets(string)] ///
		vceoption(string asis) kk(integer) vcetype(string) [weightexp(string)] ///
		[SUBOPTions(string)] [*] // [*] are ignored!

	if ("`options'"!="") Debug, level(3) msg("(ignored options: `options')")

	local vceoption = regexr("`vceoption'", "vce\( *unadjusted *\)", "vce(ols)")
	mata: st_local("vars", strtrim(stritrim( "`depvar' `indepvars' `avgevars'" )) ) // Just for esthetic purposes

* Run regression just to compute true DoF
	local subcmd _regress `vars' `weightexp', noheader notable `suboptions'
	Debug, level(3) msg("Subcommand: " in ye "`subcmd'")
	qui `subcmd'
	local N = e(N)
	local K = e(df_m) // Should also be equal to e(rank)+1
	*** scalar `sse' = e(rss)
	local WrongDoF = `N' - 1 - `K'
	local CorrectDoF = `WrongDoF' - `kk' // kk = Absorbed DoF
	Assert !missing(`CorrectDoF')
	
* Now run intended regression and fix VCV
	qui regress `vars' `weightexp', `vceoption' noheader notable `suboptions'
	* Fix DoF
	tempname V
	matrix `V' = e(V) * (`WrongDoF' / `CorrectDoF')

	* Avoid corner case error when all the RHS vars are collinear with the FEs
	if (`K'>0) {
		cap ereturn repost V=`V' // Else the fix would create MVs and we can't post
		Assert inlist(_rc,0,504), msg("error `=_rc' when adjusting the VCV")
	}
	else {
		ereturn scalar rank = 1 // Set e(rank)==1 when e(df_m)=0 , due to the constant
		* (will not be completely correct if model is already demeaned?)
	}
	
	*** if ("`vcetype'"!="cluster") ereturn scalar rank = e(rank) + `kk'

* ereturns specific to this command
	ereturn scalar df_r = max(`CorrectDoF', 0)
	mata: st_local("original_vars", strtrim(stritrim( "`original_depvar' `original_indepvars' `avge_targets' `original_absvars'" )) )
	ereturn local alternative_cmd regress `original_vars', `vceoption' `options'

	if ("`vcetype'"!="cluster") { // ("`vcetype'"=="unadjusted")
		ereturn scalar F = e(F) * `CorrectDoF' / `WrongDoF'
		if missing(e(F)) di as error "WARNING! Missing FStat"
	}


	local run_test = ("`vcetype'"=="cluster") | ( e(df_m)+1!=e(rank) )
	if (`run_test') {
		if ("`vcetype'"!="cluster") {
			Debug, level(0) msg("Note: equality df_m+1==rank failed (is there a collinear variable in the RHS?), running -test- to get correct values")
		}
		return clear
		if (`K'>0) {
			qui test `indepvars' `avge' // Wald test
			ereturn scalar F = r(F)
			ereturn scalar df_m = r(df)
			ereturn scalar rank = r(df)+1 // Add constant
		}
		else {
			ereturn scalar F = 0
			ereturn scalar df_m = 0
			ereturn scalar rank = 1
		}
	}

* Fstat
	* _U: Unrestricted, _R: Restricted
	* FStat = (RSS_R - RSS_U) / RSS * (N-K) / q
	*       = (R2_U - R2_R) / (1 - R2_U) * DoF_U / (DoF_R - DoF_U)
	Assert e(df_m)+1==e(rank) , rc(0) /// rc(322)
		msg("Error: expected e(df_m)+1==e(rank), got (`=`e(df_m)'+1'!=`e(rank)')")
end

* Cluster notes (see stata PDFs):
* We don't really have "N" indep observations but "M" (aka `N_clust') superobservations,
* and we are replacing (N-K) DoF with (M-1) (used when computing the T and F tests)
		
* For the VCV matrix, the multiplier (small sample adjustement) is q := (N-1)/(N-K) * M / (M-1)
* Notice that if every obs is its own cluster, M=N and q = N/(N-K) (the usual multiplier for -ols- and -robust-)
		
* Also, if one of the absorbed FEs is nested within the cluster variable, then we don't need to include that variable in K
* (this is the adjustment that xtreg makes that areg doesn't)
program define Wrapper_ivregress, eclass
	syntax , depvar(varname) endogvars(varlist) instruments(varlist) ///
		[indepvars(varlist) avgevars(varlist)] ///
		original_depvar(string) original_endogvars(string) original_instruments(string) ///
		[original_indepvars(string) avge_targets(string)] ///
		estimator(string) vceoption(string asis) KK(integer) [weightexp(string)] ///
		[SUBOPTions(string)] [*] // [*] are ignored!

	mata: st_local("vars", strtrim(stritrim( "`depvar' `indepvars' `avgevars' (`endogvars'=`instruments')" )) )
	if ("`estimator'"=="gmm") local vceoption = "`vceoption' " + subinstr("`vceoption'", "vce(", "wmatrix(", .)
	
	* Note: the call to -ivregress- could be optimized.
	* EG: -ivregress- calls ereturn post .. ESAMPLE(..) but we overwrite the esample and its SLOW
	* But it's a 1700 line program so let's not worry about it
	*profiler on

	local subcmd ivregress `estimator' `vars' `weightexp', `vceoption' small `suboptions'
	Debug, level(3) msg("Subcommand: " in ye "`subcmd'")
	local noise qui
	if (strpos("`options'", "first")>0) local noise noi
	`noise' `subcmd'
	
	*profiler off
	*profiler report
	
	* Fix DoF if needed
	local N = e(N)
	local K = e(df_m)
	local WrongDoF = `N' - 1 - `K'
	local CorrectDoF = `WrongDoF' - `kk'
	Assert !missing(`CorrectDoF')
	if ("`estimator'"!="gmm" | 1) {
		tempname V
		matrix `V' = e(V) * (`WrongDoF' / `CorrectDoF')
		ereturn repost V=`V'
	}
	ereturn scalar df_r = `CorrectDoF'

	* ereturns specific to this command
	mata: st_local("original_vars", strtrim(stritrim( "`original_depvar' `original_indepvars' `avge_targets' `original_absvars' (`original_endogvars'=`original_instruments')" )) )
	ereturn local alternative_cmd ivregress `estimator' `original_vars', small `vceoption' `options'
	if ("`estimator'"!="gmm") ereturn scalar F = e(F) * `CorrectDoF' / `WrongDoF'

end
program define Wrapper_ivreg2, eclass
	syntax , depvar(varname) endogvars(varlist) instruments(varlist) ///
		[indepvars(varlist) avgevars(varlist)] ///
		original_depvar(string) original_endogvars(string) original_instruments(string) ///
		[original_indepvars(string) avge_targets(string)] ///
		[original_absvars(string) avge_targets] ///
		estimator(string) vceoption(string asis) KK(integer) ///
		[SHOWRAW(integer 0) dofminus(string)] first(integer) [weightexp(string)] ///
		[SUBOPTions(string)] [*] // [*] are ignored!
	if ("`options'"!="") Debug, level(3) msg("(ignored options: `options')")
	if (`c(version)'>=12) local hidden hidden
	
	* Disable some options
	local 0 , `suboptions'
	syntax , [SAVEFPrefix(name)] [*] // Will ignore SAVEFPREFIX
	local suboptions `options'

	if ("`estimator'"=="2sls") local estimator

	if strpos("`vceoption'","unadj") {
		local vceoption
	}
	else if strpos("`vceoption'","robust")>0 {
		local vceoption "robust" 
	}
	else {
		local vceoption = substr(trim("`vceoption'"), 5, strlen("`vceoption'")-5)
		gettoken vcefirst vceoption : vceoption
		local vceoption = trim("`vceoption'")
		local vceoption `vcefirst'(`vceoption')
	}
	
	mata: st_local("vars", strtrim(stritrim( "`depvar' `indepvars' `avgevars' (`endogvars'=`instruments')" )) )
	
	if (`first') {
		local firstoption "first savefirst"
	}
	Assert inlist("`dofminus'","dofminus","sdofminus")

	local subcmd ivreg2 `vars' `weightexp', `estimator' `vceoption' `firstoption' small `dofminus'(`kk') `suboptions'
	Debug, level(3) msg(_n "call to subcommand: " _n as result "`subcmd'")
	local noise = cond(`showraw', "noi", "qui")
	`noise' `subcmd'
	if ("`noise'"=="noi") di in red "{hline 64}" _n "{hline 64}"

	if !missing(e(ecollin)) {
		di as error "endogenous covariate <`e(ecollin)'> was perfectly predicted by the instruments!"
		error 2000
	}

	if (`first') {
		ereturn `hidden' local first_prefix = "_ivreg2_"
		ereturn `hidden' local ivreg2_firsteqs = e(firsteqs)
		ereturn local firsteqs
	}

	foreach cat in exexog insts instd {
		FixVarnames `e(`cat')'
		ereturn local `cat' = "`r(newnames)'"
	}

	if (`first') {
		* May be a problem if we ran out of space for storing estimates
		local ivreg2_firsteqs "`e(ivreg2_firsteqs)'"
		tempname hold
		estimates store `hold' , nocopy
		foreach fs_eqn in `ivreg2_firsteqs' {
			qui estimates restore `fs_eqn'
			FixVarnames `e(depvar)'
			ereturn local depvar = r(prettynames)
			FixVarnames `e(inexog)'
			ereturn local inexog = r(prettynames)

			tempname b
			matrix `b' = e(b)
			local backup_colnames : colnames `b'
			FixVarnames `backup_colnames'
			matrix colnames `b' = `r(prettynames)' // newnames? prettynames?
			ereturn repost b=`b', rename

			estimates store `fs_eqn', nocopy
		}
		qui estimates restore `hold'
		estimates drop `hold'
	}

	* ereturns specific to this command
	mata: st_local("original_vars", strtrim(stritrim( "`original_depvar' `original_indepvars' `avge_targets' `original_absvars' (`original_endogvars'=`original_instruments')" )) )
	ereturn local alternative_cmd ivreg2 `original_vars', small `vceoption' `options' `estimator'
	***if ("`estimator'"!="gmm") ereturn scalar F = e(F) * `CorrectDoF' / `WrongDoF'

end


// -------------------------------------------------------------
// Display Regression Table
// -------------------------------------------------------------
 program define Replay, eclass
	syntax , [*]
	Assert e(cmd)=="reghdfe"
	local subcmd = e(subcmd)
	Assert "`subcmd'"!="" , msg("e(subcmd) is empty")

	* Add pretty names for AvgE variables
	tempname b
	matrix `b' = e(b)
	local backup_colnames : colnames `b'
	matrix colnames `b' = `e(prettynames)'
	local savefirst = e(savefirst)
	local suboptions = e(suboptions)


	local diopts = "`e(diopts)'"
	if ("`options'"!="") { // Override
		_get_diopts diopts /* options */, `options'
	}

	if ("`subcmd'"=="ivregress") {
		* Don't want to display anova table or footnote
		_coef_table_header
		_coef_table, `diopts' bmatrix(`b') vmatrix(e(V)) // plus 
	}
	else if ("`subcmd'"=="ivreg2") {
		* Backup before showing both first and second stage
		tempname hold
		
		if ("`e(ivreg2_firsteqs)'"!="") {
			estimates store `hold'

			local i 0
			foreach fs_eqn in `e(ivreg2_firsteqs)' {
				local instrument  : word `++i' of `e(instd)'
				di as input _n "{title:First stage for `instrument'}"
				estimates replay `fs_eqn' , nohead `diopts'
				if (!`savefirst') estimates drop `fs_eqn'
			}

			ereturn clear
			qui estimates restore `hold'
			di as input _n "{title:Second stage}"
		}

		estimates store `hold'
		ereturn repost b=`b', rename
		ereturn local cmd = "`subcmd'"
		`subcmd' , `diopts'
		ereturn clear // Need this because -estimates restore- behaves oddly
		qui estimates restore `hold'
		assert e(cmd)=="reghdfe"
		estimates drop `hold'


		*ereturn local cmd = "reghdfe"
		*matrix `b' = e(b)
		*matrix colnames `b' = `backup_colnames'
		*ereturn repost b=`b', rename
	}
	else {

		* Regress-specific code, because it doesn't play nice with ereturn
		sreturn clear 

		if "`e(prefix)'" != "" {
			_prefix_display, `diopts'
			exit
		}
		_coef_table_header
		di
		local plus = cond(e(model)=="ols" & e(vce)=="unadjusted", "plus", "")
		_coef_table, `plus' `diopts' bmatrix(`b') vmatrix(e(V))
	}

	reghdfe_footnote
	* Revert AvgE else -predict- and other commands will choke


end


// -------------------------------------------------------------
// Standard version of the HDFE command
// -------------------------------------------------------------
program define AlternativeCMD

    di in ye _n "[Running standard alternative to the HDFE command]"
    if !e(alternative_ok) di as error "Note: the command will be approximate only (e.g. wrong SDs, ommited AvgE)"
    di as text "[CMD] " as result "`e(alternative_cmd)'"
    local rhs `e(indepvars)' `e(endogvars)' `e(AvgE_Ws)'
    `e(alternative_cmd)'
    
    if ("`rhs'"!="") {
        di in ye _n "(Joint test on non-FE regressors)"
        testparm `rhs'
    }

end


// -------------------------------------------------------------
// Simple assertions
// -------------------------------------------------------------
program define Assert
    syntax anything(everything equalok) [, MSG(string asis) RC(integer 198)]
    if !(`anything') {
        di as error `msg'
        exit `rc'
    }
end


// -------------------------------------------------------------
// Simple debugging
// -------------------------------------------------------------
program define Debug

	syntax, [MSG(string asis) Level(integer 1) NEWline COLOR(string)] [tic(integer 0) toc(integer 0)]
	
	cap mata: st_local("VERBOSE",strofreal(VERBOSE)) // Ugly hack to avoid using a global
	if ("`VERBOSE'"=="") {
		di as result "Mata scalar -VERBOSE- not found, setting VERBOSE=3"
		local VERBOSE 3
		mata: VERBOSE = `VERBOSE'
	}


	assert "`VERBOSE'"!=""
	assert inrange(`level',0, 4)
	assert (`tic'>0) + (`toc'>0)<=1

	if ("`color'"=="") local color text
	assert inlist("`color'", "text", "res", "result", "error", "input")

	if (`VERBOSE'>=`level') {

		if (`tic'>0) {
			timer clear `tic'
			timer on `tic'
		}
		if (`toc'>0) {
			timer off `toc'
			qui timer list `toc'
			local time = r(t`toc')
			if (`time'<10) local time = string(`time'*1000, "%tcss.ss!s")
			else if (`time'<60) local time = string(`time'*1000, "%tcss!s")
			else if (`time'<3600) local time = string(`time'*1000, "%tc+mm!m! SS!s")
			else if (`time'<24*3600) local time = string(`time'*1000, "%tc+hH!h! mm!m! SS!s")
			timer clear `toc'
			local time `" as result " `time'""'
		}

		if (`"`msg'"'!="") di as `color' `msg'`time'
		if ("`newline'"!="") di
	}
end

