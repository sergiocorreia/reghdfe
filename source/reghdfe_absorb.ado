//------------------------------------------------------------------------------
// REGHDFE_ABSORB: Runs three steps required to demean wrt FEs
//------------------------------------------------------------------------------
/* TYPICAL USAGE - Five steps, with some user work in between (see Estimate.ado)

 (-)		Call -Parse-
			
 (1)	reghdfe_absorb, step(start) absorb(...)  avge(...) clustervar1(...) weight(..) weightvar(..)
			Parse absorb(), create almost-empty Mata objects
			Parse avge() and store results in Mata string vectors
			RETURN r(N_hdfe) r(N_avge) r(keepvars)

		[Until here, no data has been touched]

 (-) 		Preserve data
			Drop unused vars
			Expand factors and time series in all varlists
			Drop unused base vars of the factors vars
			Drop obs with MVs

 (2)	reghdfe_absorb, step(precompute) keepvars(...)  [depvar(...)  excludeself]
			Transform wrt -avge-
			Drop MVs caused by -excludeself-
			Fill mata objects with counts and means, delete unused vars
			RETURN r(clustervar1)

 (-)		Compute statistics such as TSS, used later on
			Save untransformed variables

 (3)	reghdfe_absorb, step(demean) varlist(...)  [maximize_options..]
			Obtain residuals of varlist wrt the FEs
			Note: can be run multiple times

 (-)		Drop IDs of the FEs if needed
			Optain DoF
			Run regressions, residuals, etc.
			Load untransformed variables
			Use predict to get resid+d

		reghdfe_absorb, step(demean) save_fe(1) var(resid_d) -> get FEs
 (4)	reghdfe_absorb, step(save) original_depvar(..)
			Save required FEs with proper name, labels, etc.

 (5)	reghdfe_absorb, step(stop)
			Clean up all Mata objects

 (-)		Restore, merge sample(), report tables
//----------------------------------------------------------------------------*/


*clear mata
include _mata/reghdfe.mata

//------------------------------------------------------------------------------
cap pr drop reghdfe_absorb
program define reghdfe_absorb, rclass
//------------------------------------------------------------------------------
	local version `clip(`c(version)', 11.2, 13.1)' // 11.2 minimum, 13+ preferred
	qui version `version'

* This allows to dump the information for one FE as locals in the caller namespace
	cap syntax, fe2local(integer) [*]
	local rc = _rc
	 if (`rc'==0) {
		local g `fe2local'
		assert inrange(`g',1,100)
		mata: fe2local(`fe2local')
		exit
	}

* Parallel instance
	cap syntax, instance [*]
	local rc = _rc
	 if (`rc'==0) {
		ParallelInstance, `options'
		exit
	}

* Parse
	syntax, STEP(string) [CORES(integer 1)] [*]

* Sanity checks
	local numstep = ("`step'"=="start") + 2*("`step'"=="precompute") + ///
		3*("`step'"=="demean") + 4*("`step'"=="save") + 5*("`step'"=="stop")
	Assert (`numstep'>0), msg("reghdfe_absorb: -`step'- is an invalid step name" _n ///
			"valid steps are: start precompute demean save stop")

	cap mata: st_local("prev_numstep", strofreal(prev_numstep))
	if (_rc) local prev_numstep 0

	Assert (`numstep'==`prev_numstep'+1) | (`numstep'==5) | ///
		(`numstep'==3 & `prev_numstep'==3) ///
		, msg("reghdfe_absorb: expected step `=`prev_numstep'+1' instead of step 	`numstep'")
	mata: prev_numstep = `numstep'
	if (`numstep'<5) Debug, msg(_n as text "{title:Running -reghdfe_absorb- step `numstep'/5 (`step')}") level(3)

* Call subroutine and return results
	if (`numstep'==1) {
		Initialize, `options'
		return local keepvars "`r(keepvars)'"
		return scalar N_hdfe = r(N_hdfe)
		return scalar N_avge = r(N_avge)
	}
	else if (`numstep'==2) {
		Prepare, `options'
		local N = r(N_clustervars)
		forval i = 1/`N' {
			assert "`r(clustervar`i')'"!=""
			return local clustervar`i' "`r(clustervar`i')'"
		}
	}
	else if (`numstep'==3 & `cores'==1) {
		Annihilate, `options'
	}
	else if (`numstep'==3 & `cores'>1) {
		AnnihilateParallel, numcores(`cores') `options'
	}
	else if (`numstep'==4) {
		Save, `options'
		return local keepvars "`r(keepvars)'"
	}
	else if (`numstep'==5) {
		Stop, `options'
	}
	else {
		error 198
	}
end

//------------------------------------------------------------------------------
cap pr drop Initialize
program define Initialize, rclass
//------------------------------------------------------------------------------
syntax, Absorb(string) [AVGE(string)] [CLUSTERVARS(string)] [OVER(varname numeric)] [WEIGHT(string) WEIGHTVAR(varname numeric)]
	Assert !regexm("`absorb'","[*?-]"), ///
		msg("error: please avoid pattern matching in -absorb-")

	if ("`over'"!="") Assert "`avge'"=="", msg("-avge- needs to be empty if -over- is used")

	Assert inlist("`weight'", "", "fweight", "aweight", "pweight")

**** ABSORB PART ****

* First pass to get the true number of FEs
	local i 0
	Debug, level(3) msg(_n "Fixed effects:")
	foreach var of local absorb {
		ParseOneAbsvar, absvar(`var')
		local i = `i' + cond(r(is_bivariate), 2, 1)
		* Output: r(target) cvars ivars is_interaction is_cont_interaction is_bivariate
		Assert `i'>1 | "`r(cvars)'"=="" | `r(is_bivariate)', ///
			msg("error parsing absorb : first absvar cannot be continuous interaction" ///
			_n "solution: i) reorder absvars, ii) replace # with ##, iii) add a constant as first absvar (as a workaround)")

		if ("`over'"!="") {
			local ivars r(ivars)
			local dupe : list ivars & over
			Assert ("`dupe'"==""), msg("-over- cannot be part of any absvar")
		}
	}

	if ("`over'"!="") {
		local ++i // We'll add -over- as the first FE
		local pre_absorb `absorb'
		local absorb `over' `absorb'
	}

* Create vector of structures with the FEs
	Assert inrange(`i',1,100), msg("error: too many absorbed variables (do not include the dummies, just the variables)")
	Debug, msg(`"(`i' absorbed fixed `=plural(`i',"effect")': "' as result "`absorb'" as text ")")
	mata: weightexp = ""
	mata: weightvar = ""
	if ("`weightvar'"!="") {
		Debug, msg(`"(`weight': "' as result "`weightvar'" as text ")")
		mata: weightexp = "[`weight'=`weightvar']"
		mata: weightvar = "`weightvar'"
		**qui cou if `fweight'<=0 | `fweight'>=. | (`fweight'!=int(`fweight'))
		** Move this somewhere else.. else it will fail needlesly if some excluded obs. have missing weights
		**Assert (`r(N)'==0), msg("fweight -`fweight'- can only have strictly positive integers (no zero, negative, MVs, or reals)!")
	}
	mata: G = `i'
	mata: initialize()

* Second pass to save the values
	local i 0
	foreach var of local absorb {
		qui ParseOneAbsvar, absvar(`over_prefix'`var')
		local keepvars `keepvars' `r(ivars)' `r(cvars)'
		local varlabel = "i." + subinstr("`r(ivars)'", " ", "#i.", .)
		if (`r(is_cont_interaction)' & !`r(is_bivariate)') local varlabel "`varlabel'#c.`r(cvars)'"
		
		local args `" "`r(target)'", "`r(ivars)'", "`r(cvars)'", `r(is_interaction)', `r(is_cont_interaction)', `r(is_bivariate)', "`weightvar'" "'
		mata: add_fe(`++i', "`varlabel'", `args', 0)
		if (`r(is_bivariate)') {
			local varlabel "`varlabel'#c.`r(cvars)'"
			mata: add_fe(`++i', "`varlabel'", `args', 1)
		}

		if ("`over'"!="") local over_prefix "i.`over'#" // Not for the first one
	}
	local N_hdfe = `i'

	if ("`over'"!="") Debug, msg("absvars expanded due to over: `pre_absorb' -> `absorb'")

**** AVGE PART ****

* First pass to get the true number of FEs
local N_avge = 0
if ("`avge'"!="") {
	local i 0
	foreach var of local avge {
		Debug, level(3) msg(_n "AvgE effects:")
		ParseOneAbsvar, absvar(`var')
		local ++i
		* Output: r(target) cvars ivars is_interaction is_bivariate
		Assert ("`r(cvars)'"=="" & `r(is_bivariate)'==0), ///
			msg("error parsing avge : continuous interactions not allowed")
	}

* Create vectors
	Assert inrange(`i',1,100), msg("error: too many avge variables (do not include the dummies, just the variables)")
	Debug, msg(`"(`i' avge `=plural(`i',"effect")': "' as result "`avge'" as text ")")
}

* Always save this to avoid not-found errors
	mata: avge_ivars = J(1, `i', "")
	mata: avge_target = J(1, `i', "")
	mata: avge_varlabel = J(1, `i', "")

* Second pass to save the values
if ("`avge'"!="") {
	local i 0
	foreach var of local avge {
		qui ParseOneAbsvar, absvar(`var')
		local ++i
		local varlabel = "i." + subinstr("`r(ivars)'", " ", "#i.", .)
		mata: avge_ivars[`i'] = "`r(ivars)'"
		mata: avge_target[`i'] = "`r(target)'"
		mata: avge_varlabel[`i'] = "`varlabel'"
		local keepvars `keepvars' `r(ivars)'
	}
	local N_avge = `i'
}
	mata: avge_num = `N_avge'

*** CLUSTER PART ****
* EG: If clustervar1=foreign, absorb=foreign, then clustervar1 -> __FE1__
	local N_clustervars : word count `clustervars'
	mata: N_clustervars = `N_clustervars'
	forv i=1/`N_clustervars' {
		local clustervar`i' : word `i' of `clustervars'
		mata: ivars_clustervar`i' = ""
		if ("`clustervar`i''"!="") {
			Debug, level(3) msg(_n "Cluster by:")
			ParseOneAbsvar, absvar(`clustervar`i'')
			Assert "`r(cvars)'"=="", msg("clustervar cannot contain continuous interactions")
			local ivars_clustervar`i' "`r(ivars)'"
			local keepvars `keepvars' `r(ivars)'
			mata: ivars_clustervar`i' = "`ivars_clustervar`i''"
		}
	}
	
**** Returns ****
	Debug, level(3) newline
	local keepvars : list uniq keepvars
	return local keepvars `keepvars'
	return scalar N_hdfe = `N_hdfe'
	return scalar N_avge = `N_avge'
end


//------------------------------------------------------------------------------
cap pr drop ParseOneAbsvar
program define ParseOneAbsvar, rclass
//------------------------------------------------------------------------------
syntax, ABSVAR(string)

	Assert !strpos("`absvar'","###"), msg("error parsing <`absvar'> : ### is invalid")
	Assert regexm("`absvar'", "^[a-zA-Z0-9_=.#]+$"), msg("error parsing <`absvar'> : illegal characters ")
	Assert !regexm("`absvar'", "##([^c]|(c[^.]))"), msg("error parsing <`absvar'> : expected c. after ##")
	local original_absvar `absvar'

* Split at equal sign
	local equalsign = strpos("`absvar'","=")
	local target = substr("`absvar'",1,`equalsign'-1)
	local absvar = substr("`absvar'",`equalsign'+1, .)
	if ("`target'"!="") conf new var `target'

	local is_interaction = strpos("`absvar'", "#")>0
	local is_bivariate = strpos("`absvar'", "##")>0

* Split interactions
	mata: st_local("vars", subinstr("`absvar'", "#", " ") )
	foreach var of local vars {

		local dot = strpos("`var'", ".")
		local root = substr("`var'", `dot'+1, .)
		unab root : `root' , max(1)
		conf numeric var `root'
		
		local prefix = substr("`var'", 1, `dot'-1)
		local prefix = lower( cond("`prefix'"=="", "i", "`prefix'") ) // -i.- is default prefix

		Assert inlist("`prefix'", "i", "c") , msg("error parsing <`absvar'><`var'> : only i. and c. allowed, not `prefix'.")
		Assert !strpos("`root'", ".") , msg("error parsing <`absvar'><`var'> : no time series operators allowed")
		
		if ("`prefix'"=="i") {
			local ivars `ivars' `root'
		}
		else {
			Assert "`cvars'"=="", msg("error: can't have more than one continuous variable in the interaction")
			local cvars `cvars' `root'
		}
	}
	local tab  "        "
	Debug, level(3) msg(as text "    Parsing " as result "`original_absvar'")
	Debug, level(3) msg(as text "`tab'ivars = " as result "`ivars'")
	if ("`cvars'"!="") Debug, level(3) msg(as text "`tab'cvars = " as result "`cvars'")
	if ("`target'"!="") Debug, level(3) msg(as text "`tab'target = " as result "`target'")
	Debug, level(3) msg(as text "`tab'is_interaction = " as result "`is_interaction'")
	Debug, level(3) msg(as text "`tab'is_bivariate = " as result "`is_bivariate'")
	// Debug, level(3) newline

	return scalar is_interaction = `is_interaction'
	return scalar is_cont_interaction = `is_interaction' & ("`cvars'"!="")
	return scalar is_bivariate = `is_bivariate'
	if ("`target'"!="") return local target "`target'"
	if ("`cvars'"!="") return local cvars "`cvars'"
	return local ivars "`ivars'"
end


//------------------------------------------------------------------------------
cap pr drop Prepare
program define Prepare, rclass
//------------------------------------------------------------------------------
syntax, KEEPvars(varlist) [DEPVAR(varname numeric) EXCLUDESELF]

**** AVGE PART ****
mata: st_local("N_avge", strofreal(avge_num))
if (`N_avge'>0) {
	forv g=1/`N_avge' {
		Assert ("`depvar'"!=""), msg("reghdfe_absorb: depvar() required")
		mata: st_local("ivars", avge_ivars[`g'])
		mata: st_local("varlabel", avge_varlabel[`g'])
		mata: st_local("target", avge_target[`g'])
		local W __W`g'__

		local note = cond("`excludeself'"=="",""," (excluding obs. at hand)")
		local original_depvar = cond(substr("`depvar'",1,2)=="__", "`: var label `depvar''", "`depvar'")
		Debug, level(2) msg(" - computing AvgE(`original_depvar') wrt (`varlabel')`note'")

		* Syntax: by ... : AverageOthers varname , Generate(name) EXCLUDESELF
		qui AverageOthers `depvar', by(`ivars') gen(`W') `excludeself'
		char `W'[target] `target'
	}

	* Marked obs should have been previously deleted
	tempvar touse
	mark `touse'
	markout `touse' __W*__
	qui keep if `touse'
	drop `touse'
	local keepvars `keepvars' __W*__
}

	Assert c(N)>0, rc(2000) msg("Empty sample, check for missing values or an always-false if statement")
	Assert c(N)>1, rc(2001)

**** ABSORB PART ****

* Create compact IDs and store them into the __FE1__ (and so on)
* Also structures the IDs in the 1..K sense, which is needed by mata:prepare
* (b/c it uses the largest value of the ID as the number of categories)
	mata: st_local("G", strofreal(G))
	mata: st_local("N_clustervars", strofreal(N_clustervars))
	forval i = 1/`N_clustervars' {
		mata: st_local("ivars_clustervar`i'", ivars_clustervar`i')
		local clustervar`i' `ivars_clustervar`i'' // Default if cluster not a FE
	}

	* Get list of cvars to avoid this bug:
	* i.t i.zipcode##c.t -> -t- is missing
	forv g=1/`G' {
		reghdfe_absorb, fe2local(`g') // Dump FEs[g] data into locals
		local all_cvars `all_cvars' `cvars'
	}

	
	forv g=1/`G' {
		reghdfe_absorb, fe2local(`g') // Dump FEs[g] data into locals
		if (`is_mock') continue

		Debug, level(2) msg(" - creating compact IDs for categories of `varlabel' -> " as result "__FE`g'__")

		local ivar_is_cvar : list ivars & all_cvars
		local ivar_is_cvar = ("`ivar_is_cvar'"!="")

		if (`is_interaction' | `ivar_is_cvar') {
			GenerateID `ivars',  gen(`varname')
		}
		else {
			* Can't rename it now b/c it may be used in other absvars
			GenerateID `ivars' , replace
			local rename`g' rename `ivars' `varname'
		}

		local is_cluster 0
		forval i = 1/`N_clustervars' {
			if ("`ivars_clustervar`i''"!="") {
				local is_cluster : list ivars_clustervar`i' === ivars
				*di in red "`ivars_clustervar`i'' === `ivars'"
				if (`is_cluster') {
				 	Debug, level(3) msg(" - clustervar`i': " as result "`ivars_clustervar`i''" as text " -> " as result "`varname'")
				 	local clustervar`i' `varname'
				}			
			}
		}
		local all_cvars `all_cvars' `cvars'
	}

* Create clustervar if needed
	forval i = 1/`N_clustervars' {
		if (!`is_cluster' & `: word count `ivars_clustervar`i'''>1) {
			*local clustervar1 = subinstr("`ivars_clustervar1'", " ", "_", .)
			*local clustervar1 : permname `clustervar1'
			local clustervar`i' __clustervar`i'__
			GenerateID `ivars_clustervar`i'',  gen(`clustervar`i'')
			// Optional: Use GenerateID even with one ivar; will likely save space
		}
	}

* Rename all absvars into __FE1__ notation
	forv g=1/`G' {
		reghdfe_absorb, fe2local(`g')
		if (`is_mock') continue

		`rename`g''
		qui su __FE`g'__, mean
		local K`g' = r(max)
		Assert `K`g''>0
		local name : char __FE`g'__[name]
		local summarize_fe `"`summarize_fe' as text " `name'=" as result "`K`g''" "'
	}

	forval i = 1/`N_clustervars' {
		assert "`clustervar`i''"!=""
		return local clustervar`i' "`clustervar`i''"
		local clustervars `clustervars' `clustervar`i''
	}
	return scalar N_clustervars = `N_clustervars'

	keep __FE*__ `all_cvars' `clustervars' `keepvars' `weightvar'
	Debug, level(1) msg("(number of categories by fixed effect:" `summarize_fe' as text ")") newline

* Fill in auxiliary Mata structures
	Debug, level(2) tic(20)
	mata: prepare()
	Debug, level(2) toc(20) msg("mata:prepare took")
end


//------------------------------------------------------------------------------
cap pr drop program Annihilate
program Annihilate
//------------------------------------------------------------------------------
syntax , VARlist(varlist numeric) ///
	[TOLerance(real 1e-7) MAXITerations(integer 1000) ACCELerate(integer 1) /// See reghdfe.Parse
	CHECK(integer 0) SAVE_fe(integer 0) /// Runs regr of FEs
	NUM_fe(integer -1)] /// Regress only against the first Nth FEs (used in nested Fstats)
	[bad_loop_threshold(integer 1) stuck_threshold(real 5e-3) pause_length(integer 20) ///
	accel_freq(integer 3) accel_start(integer 6)] /// Advanced options

	assert inrange(`tolerance', 1e-20, 1) // However beyond 1e-16 we reach the limits of -double-
	assert inrange(`maxiterations',1,.)
	assert inlist(`accelerate',0,1)
	assert inlist(`check',0,1)
	assert inlist(`save_fe',0,1)
	assert inrange(`num_fe',1,100) | `num_fe'==-1 // -1 ==> Use all FEs

	assert `bad_loop_threshold'>0
	assert `stuck_threshold'>0 & `stuck_threshold'<=1
	assert `pause_length'>=0
	assert `accel_freq'>=0
	assert `accel_start'>0

	* We need to recast everything to -double- (-float- is not good enough)
	Debug, level(2) msg("(recasting variables as -double-)")
	recast double `varlist'

	* We can't save the FEs if there is more than one variable
	cap unab _ : `varlist', max(1)
	Assert (_rc==0 | `save_fe'==0) , rc(`=_rc') ///
		msg("reghdfe_absorb: cannot save FEs of more than one variable at a time")

	tempvar resid
	local save = `save_fe' | `check' // check=1 implies save_fe=1
	local base_args `" "`resid'", `tolerance', `maxiterations', `save', `accelerate', `num_fe'  "'
	local adv_args `" `bad_loop_threshold', `stuck_threshold', `pause_length', `accel_freq', `accel_start' "'
	local args `" `base_args' , `adv_args' "'
	* di in red `"<`args'>"'

	Debug, level(2) tic(30)
	mata: st_local("weightexp", weightexp)
	
	foreach var of varlist `varlist' {
		cap drop __Z*__
		Assert !missing(`var'), msg("reghdfe_absorb: `var' has missing values and cannot be transformed")
		
		* Syntax: MAKE_RESIDUAL(var, newvar, tol, maxiter | , save=0 , accel=1, first_n=`num_fe')
		* Note: summarize doesn't allow pweight ( see http://www.stata.com/support/faqs/statistics/weights-and-summary-statistics/ )
		* Since we only want to compute means, replace with [aw]
		local tmpweightexp = subinstr("`weightexp'", "[pweight=", "[aweight=", 1)
		qui su `var' `tmpweightexp', mean
		char define `var'[mean] `r(mean)'
		mata: make_residual("`var'", `args')
		assert !missing(`resid')

		* Check that coefs are approximately 1
		if (`check') {
			unab _ : __Z*__, min(1)
			local backup = ("`e(cmd)'"!="")
			if (`backup') {
				tempname backup_results
				est store `backup_results', nocopy // nocopy needed to avoid having e(_estimates_name)
			}
			qui _regress `var' __Z*__
			local label : var label `var'
			if ("`label'"=="") local label `var'
			di as text "FE coefficients for `label':{col 36}" _continue
			foreach z of varlist __Z*__ {
				assert !missing(`z')
				di as text " `=string(_b[`z'], "%9.7f")'"  _continue
			}
			di
			
			if (`backup') qui est restore `backup_results'
			if (!`save_fe') cap drop __Z*__
		}

		* If the tol() is not high enough (e.g. 1e-14), we may fail to detect variables collinear with the absorbed categories
		* Again, we can't use pweight with summarize, but in this case it's just for debugging purposes so use [aw]
		qui su `resid' `tmpweightexp'
		local prettyvar `var'
		if (substr("`var'", 1, 2)=="__") local prettyvar : var label `var'
		if inrange(r(sd), 1e-20 , epsfloat()) di in ye "(warning: variable `prettyvar' is probably collinear, maybe try a tighter tolerance)"

		qui replace `var' = `resid' // This way I keep labels and so on
		drop `resid'
		Assert !missing(`var'), msg("REGHDFE.Annihilate: `var' has missing values after transformation")
	}
	Debug, level(2) toc(30) msg("(timer for calls to mata:make_residual)")
end


//------------------------------------------------------------------------------
cap pr drop Save
program Save, rclass
//------------------------------------------------------------------------------
// Run this after -Annihilate .. , save_fe(1)-
// For each FE, if it has a -target-, add label, chars, and demean or divide
syntax , original_depvar(string)

	mata: st_local("G", strofreal(G))
	mata: st_local("weightexp", weightexp)
	forv g=1/`G' {

		// ivars cvars target varname varlabel is_interaction is_cont_interaction is_bivariate is_mock levels
		reghdfe_absorb, fe2local(`g')
		if ("`target'"=="") continue

		* Rename, add label and chars
		rename __Z`g'__ `target'
		local label `varlabel'
		la var `target' "Fixed effects of `label' on `original_depvar'"
		char `target'[label] `label'
		char `target'[levels] `levels'

		* Substract mean, or divide by cvar (fixing division by zero errors)
		if ("`cvars'"!="" & !(`is_bivariate' & !`is_mock')) {
			char `target'[cvars] `cvars'
			qui replace `target' = cond(abs(`cvars')<epsfloat(), 0,  `target'/`cvars')
			// BUGBUG BUGBUG float(`target'/`cvars')) -> this makes them have the same FE but loses precision!
		}
		else {
			qui su `target' `weightexp', mean
			// qui replace `target' = `target' - r(mean)
			// BUGBUG BUGBUG
		}

		local keepvars `keepvars' `target'
	}

	cap drop __Z*__
	return local keepvars " `keepvars'" // the space prevents MVs
end


//------------------------------------------------------------------------------
cap pr drop Stop
program Stop
//------------------------------------------------------------------------------
	cap mata: mata drop prev_numstep // Created at step 1
	cap mata: mata drop VERBOSE // Created before step 1
	cap mata: mata drop G // Num of absorbed FEs
	cap mata: mata drop FEs // Main Mata structure
	cap mata: mata drop betas // Temporary matrices used to store bi/multivariate regr coefs
	cap mata: mata drop varlist_cache // Hash table with the names of the precomputed residuals
	cap mata: mata drop avge_* // Drop AvgE structures
	cap mata: mata drop weightexp weightvar
	
	local N_clustervars 0
	cap mata: st_local("N_clustervars", strofreal(N_clustervars))
	forval i = 1/`N_clustervars' {
		cap mata: mata drop ivars_clustervar`i'
	}
	cap mata: mata drop N_clustervars


	if ("${reghdfe_pwd}"!="") {
		qui cd "${reghdfe_pwd}"
		global reghdfe_pwd
	}

	* PARALLEL SPECIFIC CLEANUP
	cap mata: st_local("path", parallel_path)
	if ("`path'"!="") {
		mata: st_local("cores", strofreal(parallel_cores))
		assert "`cores'"!=""
		local path "`path'"
		cap erase `"`path'hdfe_mata.mo"'
		forv core=1/`cores' {
			cap erase `"`path'`core'_done.txt"'
			cap erase `"`path'`core'_ok.txt"'
			cap erase `"`path'`core'_error.txt"'
			cap erase `"`path'`core'_output.dta"'
			cap erase `"`path'`core'_log.log"'
		}
		cap rmdir `"`path'"'
		cap mata: mata drop parallel_cores
		cap mata: mata drop parallel_dta
		cap mata: mata drop parallel_vars
		cap mata: mata drop parallel_opt
		cap mata: mata drop parallel_path
	}
end

//------------------------------------------------------------------------------
cap pr drop ParallelInstance
program ParallelInstance
//------------------------------------------------------------------------------
	syntax, core(integer) code(string asis)
	set more off
	assert inrange(`core',1,32)
	local path "`c(tmpdir)'reghdfe_`code'`c(dirsep)'"
	cd "`path'"
	set processors 1

	file open fh using "`core'_started.txt" , write text all
	file close _all

	cap noi {
		set linesize 120
		log using `core'_log.log, text

		mata: mata matuse "hdfe_mata.mo"
		mata: st_local("cores",strofreal(parallel_cores))
		assert `core' <= `cores'
		mata: st_local("usedta",parallel_dta)
		mata: st_local("vars",parallel_vars[`core'])
		mata: st_local("weightvar",weightvar)
		mata: st_local("opt",parallel_opt)
		Debug, msg(" - This is core `core'/`cores'")
		sleep 100
	
		local outfn "`core'_output.dta"
		conf new file "`outfn'"

		use `vars' `weightvar' using "`usedta'"
		de, full
		reghdfe_absorb, step(demean) varlist(`vars') `opt'
		keep `vars'
		save `"`outfn'"'
		log close _all
	}

	local rc = _rc
	sleep 100

	if `rc'>0 {
		di in red "ERROR: `rc'"
		file open fh using "`core'_error.txt" , write text all
		file close _all
	}
	else {
		file open fh using "`core'_ok.txt" , write text all
		file close _all
	}

	file open fh using "`core'_done.txt" , write text all
	file close _all
	exit, STATA
end

//------------------------------------------------------------------------------
cap pr drop AnnihilateParallel
program AnnihilateParallel
//------------------------------------------------------------------------------
// Notes:
// First cluster is taking by this stata instance, to save HDD/memory/merge time
// Also, this cluster should have more obs than the other ones so we let it have
// the default number of processes
// (the other start with 1 proc allowed, which should be fine)
// Thus it will usually finish faster, to start waiting for the 2nd fastest  to merge

syntax, VARlist(varlist numeric) FILEname(string) UID(varname numeric) Numcores(integer) [*]

	local varlist : list uniq varlist
	local K : list sizeof varlist
	local numcores = min(`numcores',`K')
	local size = c(N) * c(width) / 2^30
	local wait = int(100 + 1000 * `size') // each gb wait 1 sec

	* Deal each variable like cards in Poker
	local core 1
	foreach var of local varlist {
		local varlist`core' `varlist`core'' `var'
		local ++core
		if (`core'>`numcores') local core 1
	}

	* Folder name.. need some entropy.. use varlist + time
	mata: st_local("hash", strofreal(hash1("`varlist'"), "%20.0f"))
	local seed = real(subinstr(c(current_time),":","",.)) + `hash'
	local seed = mod(`seed',2^30) // Needs to be < 2^31-1
	set seed `seed'
	local code = string(int( uniform() * 1e6 ), "%08.0f")

	* Prepare
	local path "`c(tmpdir)'reghdfe_`code'`c(dirsep)'"
	Debug, level(1) msg(" - tempdir will be " as input "`path'")
	mata: parallel_cores = `numcores'
	mata: parallel_dta = `"`filename'"'
	mata: parallel_vars = J(`numcores',1,"")
	mata: parallel_vars = J(`numcores',1,"")
	mata: parallel_opt = `"`options'"'
	mata: parallel_path = `"`path'"'
	forv i=1/`numcores' {
		mata: parallel_vars[`i'] = "`varlist`i''"
	}

	local dropvarlist : list varlist - varlist1
	drop `dropvarlist' // basically, keeps UID and clustervar
	mata: st_global("reghdfe_pwd",pwd())
	mkdir "`path'"
	qui cd "`path'"

	local objects VERBOSE G FEs betas prev_numstep parallel_* weightexp weightvar
	qui mata: mata matsave "`path'hdfe_mata.mo" `objects' , replace

	* Call -parallel-
	Debug, level(1) msg(" - running parallel instances")
	qui mata: parallel_setstatadir("")
	local binary `"$PLL_DIR"'
	global PLL_DIR

	cap mata: st_local("VERBOSE",strofreal(VERBOSE))
	if (`VERBOSE'==0) local qui qui
	`qui' di as text _n 44 * "_" + "/ PARALLEL \" + 44 * "_"

	* Create instances
	forv i=2/`numcores' {
		local cmd `"winexec `binary' /q  reghdfe_absorb, instance core(`i') code(`code') "'
		Debug, level(1) msg(" - Executing " in ye `"`cmd' "')
		`cmd'
		Debug, level(1) msg(" - Sleeping `wait'ms")
		if (`i'!=`numcores') sleep `wait'
	}
	reghdfe_absorb, step(demean) varlist(`varlist1') `options' // core=1

	* Wait until all instances have started
	local timeout 20
	local elapsed 0
	forv i=2/`numcores' {
		local ok 0
		while !`ok' {
			sleep 100
			local fn "`path'`i'_started.txt"
			cap conf file "`fn'"
			local rc = _rc
			if (`rc'==0) {
				local ok 1
				Debug, level(1) msg(" - process `i' started")
				erase "`fn'"
			}
			else {
				local elapsed = `elapsed' + 0.1
				Assert `elapsed'<`timeout', msg("Failed to start subprocess `i'")
			}
		}
		local cores `cores' `i' // Will contain remaining cores
	}

	* Wait for termination and merge
	while ("`cores'"!="") {
		foreach core of local cores {
			local donefn "`path'`core'_done.txt"
			local okfn "`path'`core'_ok.txt"
			local errorfn "`path'`core'_error.txt"
			local dtafn "`path'`core'_output.dta"
			local logfile "`path'`core'_log.log"


			cap conf file "`donefn'"
			local rc = _rc

			if (`rc'==0) {
				Debug, level(1) msg(" - process `core' finished")
				erase "`donefn'"

				cap conf file "`okfn'"
				if (`=_rc'>0) {
					type "`logfile'"
					//di as error "<`dtafn'> not found"
					Assert 0, msg("Call to subprocess `core' failed, see logfile")
				}

				erase "`okfn'"
				Debug, level(1) msg(" - Subprocess `core' done")
				local cores : list cores - core
				mata: st_local("VERBOSE",strofreal(VERBOSE))
				
				if (`VERBOSE'>=3) {
					type "`logfile'"
				}
				erase "`logfile'"

				* Merge file
				Debug, level(1) msg(" - Merging dta #`core'")
				merge 1:1 _n using "`dtafn'", nogen nolabel nonotes noreport sorted assert(match)
				erase "`dtafn'"
			}
			else {
				sleep 500 // increase this
			}
		}
	}

	* Cleanup
	qui cd "${reghdfe_pwd}"
	erase "`path'hdfe_mata.mo"
	cap rmdir `"`path'"'
	`qui' di as text 44 * "_" + "\ PARALLEL /" + 44 * "_"

end

//------------------------------------------------------------------------------

include "_reghdfe_absorb/GenerateID.ado"
include "_reghdfe_absorb/AverageOthers.ado"
include "_common/Assert.ado"
include "_common/Debug.ado"
