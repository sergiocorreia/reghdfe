*! version 6.12.3 08aug2023
program define reghdfe
	
	* Intercept storing alphas
	cap syntax, store_alphas
	if (!c(rc)) {
		Store_Alphas
		exit
	}

	* Super secret option
	cap syntax, shrug
	if (!c(rc)) {
		di as text _n `"    {browse "https://www.theawl.com/2014/05/the-life-and-times-of-%C2%AF_%E3%83%84_%C2%AF/":¯\_(ツ)_/¯}"'
		exit
	}

	* Intercept calls to parallel workers
	* NOTE: EXPERIMENTAL FEATURE
	cap syntax, worker [*]
	if (!c(rc)) {
		ParallelWorker, `options'
		exit
	}

	* Intercept previous versions of reghdfe
	cap syntax anything(everything) [fw aw pw/], [*] VERSION(integer) [noWARN]
	if !c(rc) {

		_assert inlist(`version', 3, 5)
		* reghdfe 3: 3.2.9 (21feb2016)
		* reghdfe 5: 5.9.0 (03jun2020)

		if ("`warn'" != "nowarn") di as error "(running historical version of reghdfe: `version')"
		if ("`weight'"!="") local weightexp [`weight'=`exp']

		if (`version' == 3) {
			reghdfe3 `anything' `weightexp', `options'
		}
		else {
			reghdfe5 `anything' `weightexp', `options'
		}
		exit
	}

	* Intercept replay of results
	* replay() is True if the first non-blank character of `0' is a comma
	if replay() {
		Replay `0'
		exit
	}

	* Main use case
	loc keep_mata 0
	Cleanup 0 `keep_mata'
	qui which ftools // ms_get_version
	ms_get_version ftools, min_version("2.49.1")
	cap noi Estimate `0' // Estimate will change the value of `keep_mata' if it receives noregress
	Cleanup `c(rc)' `keep_mata'
end


program Cleanup
	* Clean up globals, variables, mata objects, etc. Raise error if required
	args rc keep_mata

	* Cleanup after parallel
	loc cleanup_folder = !`keep_mata' & ("$LAST_PARALLEL_DIR"!="")
	if (`cleanup_folder') cap mata: unlink_folder(HDFE.parallel_dir, 0) // (folder, verbose)
	global LAST_PARALLEL_DIR
	global pids
	
	* Main cleanup
	if (!`keep_mata') cap mata: mata drop HDFE
	cap mata: mata drop hdfe_*
	cap drop __temp_reghdfe_resid__

	if (`rc') exit `rc'
end


program Replay, rclass
	syntax [, * noHEADer noTABLE noFOOTnote]
	if (`"`e(cmd)'"' != "reghdfe") error 301
	_get_diopts options, `options'

	if ("`header'" == "") reghdfe_header // ReghdfeHeader // _coef_table_header // reghdfe_header
	if ("`header'" == "" & "`table'" == "") di ""
	if ("`table'" == "") _coef_table, `options'
	return add // adds r(level), r(table), etc. to ereturn (before the footnote deletes them)
	if ("`footnote'" == "") ReghdfeFootnote
end


program Estimate, eclass

// --------------------------------------------------------------------------
// Parse syntax
// --------------------------------------------------------------------------
	* We try to keep consistency with "help maximize" and "help estimation options"

	#delimit ;
	syntax varlist(fv ts numeric) [if] [in] [fw aw pw/] [ ,
		
		/* Standard FEs */
		Absorb(string)
	
		/* Team FEs */
		Group_id(varname numeric)
		Individual_id(varname numeric)
		AGgregation(string)
		
		/* Model */
		VCE(string) CLuster(string)
		RESiduals(name) RESiduals2 	/* if no name, residuals saved as _reghdfe_resid */

		/* Degrees-of-freedom Adjustments */
		DOFadjustments(string)
		GROUPVar(name) 				/* var with the first connected group between FEs */

		/* Optimization */
		TEChnique(string) 			/* lsmr lsqr map gt */
		TOLerance(real 1e-8)
		ITERATE(real 16000)			/* */
		
		/* Optimization: MAP */
		TRAnsform(string)
		ACCELeration(string)

		/* Optimization: LSMR and LSQR */
		PREConditioner(string)
		
		/* Optimization: GT */
		PRUNE /* (currently disabled) */

		/* Memory usage */
		NOSAMPle					/* do not save e(sample) */
		COMPACT						/* use as little memory as possible but is slower */
		POOLsize(integer 10) 		/* Process variables in batches of # ; 0 turns it off */

		/* Multiprocessing (parallel computing) */
		PARallel(string asis)

		/* Extra display options (inspired on regress options) */
		noHEader noTABle noFOOTnote

		/* Debugging */
		Verbose(integer 0) noWARN
		TIMEit
		
		/* Undocumented */
		KEEPSINgletons

		/* Programming options (for other packages, etc) */
		noPARTIALout /* stop before partialling out; will keep the mata object */
		varlist_is_touse /* is used by ivreghdfe together with nopartialout; UNDOCUMENTED */
		noREGress /* will stop before regressing; will keep the mata HDFE object, will allow extra options that will be kept in HDFE.extra_options */
		KEEPMATA /* keep the mata HDFE object */
		FASTREGress /* Use a much faster but less numerically accurate regression method */

		/* Kept for backward-compatibility with reghdfe v5.9 */
		noCONstant /* report constant; enabled by default as otherwise -margins- fails */
		noAbsorb2 /* doesn't do anything; suffix it with "2" so it doesn't override absorb() */

		/* Options parsed later:
			a) display options: _get_diopts
			b) advanced options that correspond to properties of the FixedEffects Mata class
			EG: miniter(#) abort(0|1) ; min_ok(#) accel_start(#) accel_freq(#)
		*/
		] [*]
		;
	#delimit cr


// --------------------------------------------------------------------------
// Validate options
// --------------------------------------------------------------------------

	loc timeit = ("`timeit'"!="")
	if (`timeit') timer on 20 // total time

	if (`verbose' >= 2) di _n `"{txt}{bf:[CMD]} {inp}reghdfe `0'"'

	cap drop __hdfe* // previously saved alphas (fixed effect coefs)

	if (`verbose' > 0) di as text "{title:Parsing and validating options:}" _n
	
	* Store display options in `diopts', keep the rest in `options'
	_get_diopts diopts options, `options'

	* Convert some options to boolean
	loc drop_singletons = ("`keepsingletons'" == "")
	loc compact = ("`compact'" != "")
	loc has_standard_fe = (`"`absorb'"' != "")
	loc report_constant = "`constant'" != "noconstant"
	loc has_teams = (`"`group_id'"' != "")
	loc has_individual_fe = (`"`individual_id'"' != "")

	* Convert programmer options to boolean
	loc stop_before_partial_out = ("`partialout'" == "nopartialout")
	loc stop_before_regression = ("`regress'" == "noregress")
	loc fast_regression = ("`fastregress'" == "fastregress")

	if (`has_individual_fe') _assert `has_teams', msg("cannot set the individual() identifiers without the group() identifiers") rc(198)

	* Default values
	if ("`technique'" == "") loc technique = cond("`individual_id'"=="", "map", "lsmr") // note: map doesn't yet support indiv FEs
	if ("`transform'" == "") loc transform "symmetric_kaczmarz"
	if ("`acceleration'" == "") loc acceleration "conjugate_gradient"
	if ("`preconditioner'" == "") loc preconditioner "block_diagonal"
	if (`poolsize' == 0) loc poolsize = .

	* Warn if keeping singletons
	if (`verbose'>-1 & "`keepsingletons'"!="" & "`warn'" != "nowarn") {
		loc url "http://scorreia.com/reghdfe/nested_within_cluster.pdf"
		loc msg "WARNING: Singleton observations not dropped; statistical significance is biased"
		di as error `"`msg' {browse "`url'":(link)}"'
	}
	
	* Allow cluster(vars) as a shortcut for vce(cluster vars)
	if ("`cluster'"!="") {
		_assert ("`vce'"==""), msg("only one of cluster() and vce() can be specified") rc(198)
		loc vce cluster `cluster'
	}

	* How are individual FEs added up
	if ("`aggregation'" == "") loc aggregation mean
	if ("`aggregation'" == "average" | "`aggregation'" == "avg") loc aggregation mean
	_assert inlist("`aggregation'", "mean", "sum")
	loc function_individual "`aggregation'" // this is the name used in HDFE()


// --------------------------------------------------------------------------
// Parse all options except absorb()
// --------------------------------------------------------------------------

	* Split varlist into <depvar> and <indepvars>
	if (`verbose' > 0) di as text "# Parsing varlist: {res}`varlist'" _c
	ms_parse_varlist `varlist'
	if (`verbose' > 0) return list
	loc depvar `r(depvar)'
	loc indepvars `r(indepvars)'
	loc fe_format "`r(fe_format)'"
	loc basevars `r(basevars)'

	* Parse weights (weight type saved in `weight'; weight var in `exp')
	if ("`weight'"!="") unab exp : `exp', min(1) max(1) // simple weights only, no expressions

	* Parse VCE
	if (`verbose' > 0) di as text _n "# Parsing vce({res}`vce'{txt})" _c
	ms_parse_vce, vce(`vce') weighttype(`weight')
	if (`verbose' > 0) sreturn list
	loc vcetype `s(vcetype)'
	loc clustervars `s(clustervars)' // e.g. "exporter#importer"
	loc base_clustervars `s(base_clustervars)' // e.g. "exporter importer"
	loc num_clusters = `s(num_clusters)'
	confirm /*numeric*/ variable `base_clustervars', exact // usually clustervar can be strings, but it's way easier to disallow it

	if (`stop_before_partial_out' & "`varlist_is_touse'" != "") {
		* We will just create and fill the HDFE object
		loc touse `varlist'
		loc varlist
		markout `touse' `base_clustervars' `group_id' `individual_id'
	}
	else {
		* Set sample based on if + in + weight + group + individual + cluster variables + regression varlist (i.e. excluding absvars)
		loc varlist `depvar' `indepvars' `base_clustervars' `group_id' `individual_id'
		marksample touse, strok
		la var `touse' "[touse]"
	}
	
	if (`stop_before_partial_out') loc varlist // clear out the varlist as it is only used to set -touse-

	* Algorithm
	loc valid_techniques map cg lsmr lsqr
	_assert (`: list technique in valid_techniques'), msg("invalid technique: `technique'")

	* Transforms: allow abbreviations (cim --> cimmino)
	loc transform = lower("`transform'")
	loc valid_transforms cimmino kaczmarz symmetric_kaczmarz rand_kaczmarz
	foreach x of local valid_transforms {
		if (strpos("`x'", "`transform'")==1) loc transform `x'
	}
	_assert (`: list transform in valid_transforms'), msg("invalid transform: `transform'")

	* Accelerations
	loc acceleration = lower("`acceleration'")
	if ("`acceleration'"=="cg") loc acceleration conjugate_gradient
	if ("`acceleration'"=="sd") loc acceleration steepest_descent
	if ("`acceleration'"=="off") loc acceleration none
	loc valid_accelerations conjugate_gradient steepest_descent aitken none
	foreach x of local valid_accelerations {
		if (strpos("`x'", "`acceleration'")==1) loc acceleration `x'
	}
	_assert (`: list acceleration in valid_accelerations'), msg("invalid acceleration: `acceleration'")

	* Preconditioners
	loc preconditioner = lower("`preconditioner'")
	if ("`preconditioner'"=="off") loc preconditioner none
	loc valid_preconditioners none diagonal block_diagonal
	foreach x of local valid_preconditioners {
		if (strpos("`x'", "`preconditioner'")==1) loc preconditioner `x'
	}
	_assert (`: list preconditioner in valid_preconditioners'), msg("invalid preconditioner: `preconditioner'")


	 * Parse DoF Adjustments
	if (`verbose' > 0) di as text _n `"# Parsing dof({res}`dofadjustments'{txt})"' _c
	ParseDOF, `dofadjustments' // s(dofadjustments)
	loc dofadjustments `s(dofadjustments)'
	if (`verbose' > 0) sreturn list

	* Residuals
	 opts_exclusive "`residuals' `residuals2'" residuals
	if ("`residuals2'" != "") {
		cap drop _reghdfe_resid // destructive!
		loc residuals _reghdfe_resid
	}
	else if ("`residuals'"!="") {
		conf new var `residuals'
	}

	* Parallel
	if (`"`parallel'"' != "") {
		if (`verbose' > 0) di as text _n `"# Parsing parallel options: {inp}`parallel'"' _c
		ParseParallel `parallel'
		if (`verbose' > 0) sreturn list
		loc parallel_maxproc `s(parallel_maxproc)' // If this is zero, no parallel code will be run
		loc parallel_dir `"`s(parallel_dir)'"'
		loc parallel_force `s(parallel_force)'
		loc parallel_opts `"`s(parallel_opts)'"' // options passed to -parallel_map-
	}
	else {
		loc parallel_maxproc 0
		loc parallel_force 0
	}


	* 1) Assert that all variables (depvar, regressors, weights, etc.) have the same value with a given group
	* 2) Assert that the same individual does not appear twice within a group
	* 3) Restrict 'touse' to only include one observation per value of group_id
	* 4) Generate a 'touse_individual' to be used with the individual FE variables
	* We do these together to avoid sorting multiple times
	if (`has_teams') {
		tempvar indiv_tousevar
		ValidateGroups `basevars' `base_clustervars' `exp', group_id(`group_id') touse(`touse') indivtouse(`indiv_tousevar') individual(`individual_id')
		_assert ("`weight_type'"=="fweight") + ("`indiv_tousevar'" != "") < 2, msg("fweights are incompatible with individual ids as there cannot be two observations for a given group-individual touple")
	}

	* Parse absorb(), construct Mata object "fixed_effects", update -touse- accordingly
	* Syntax: absorb(..., NOIisly SAVEfe|Generate Replace)
	mata: HDFE = FixedEffects()
	
	if (`verbose' > 0) di as text _n `"# Passing main options to Mata"' _n

	* Pass string options to HDFE object
	loc absvars `"`absorb'"' // trick
	loc tousevar `"`touse'"' // trick
	loc weight_type `"`weight'"' // trick
	loc weight_var `"`exp'"' // trick
	loc optim_options absvars tousevar weight_type weight_var technique transform acceleration preconditioner parallel_dir parallel_opts
	if (`has_teams') loc optim_options `optim_options' group_id individual_id indiv_tousevar function_individual
	foreach opt of local optim_options {
		if (`verbose' > 0) di as text `"    - HDFE.`opt' = {res}`"``opt''"' "'
		mata: HDFE.`opt' = `"``opt''"'
	}

	* Pass numeric options to HDFE object
	loc maxiter = `iterate' // prefer to call it maxiter within Mata
	loc optim_options drop_singletons tolerance maxiter compact poolsize verbose parallel_maxproc parallel_force timeit // prune 
	foreach opt of local optim_options {
		if (`verbose' > 0) di as text `"    - HDFE.`opt' = {res}``opt''"'
		mata: HDFE.`opt' = ``opt''
	}

	if (`verbose' > 0) di as text _n `"# Parsing absorb({res}`absorb'{txt}) and initializing FixedEffects() object"'
	if (`timeit') timer on 21
	mata: HDFE.init()
	if (`timeit') timer off 21

	* Pass undocumented numeric options to HDFE object
	mata: add_undocumented_options("HDFE", `"`options'"', `verbose')

	* Preserve memory
	if (`compact') {
		loc panelvar "`_dta[_TSpanel]'"
		loc timevar "`_dta[_TStvar]'"

		cap conf var `panelvar', exact
		if (c(rc)) loc panelvar
		mata: HDFE.panelvar = "`panelvar'"

		cap conf var `timevar', exact
		if (c(rc)) loc timevar
		mata: HDFE.timevar = "`timevar'"

		if (`verbose' > 0) di as text "## Preserving dataset"
		preserve
		novarabbrev keep `basevars' `base_clustervars' `panelvar' `timevar' `touse' // `exp'
	}

	* Fill out VCE information
	mata: HDFE.vcetype = "`vcetype'"
	mata: HDFE.num_clusters = `num_clusters' // used when computing degrees-of-freedom
	mata: HDFE.clustervars = tokens("`clustervars'")
	mata: HDFE.base_clustervars = tokens("`base_clustervars'")
	
	* Compute degrees-of-freedom (either exact or conservative lower bound)
	if (`timeit') timer on 22
	mata: estimate_dof(HDFE, tokens("`dofadjustments'"), "`groupvar'")
	if (`timeit') timer off 22

	if (`stop_before_partial_out') {
		if (`verbose' > 0) di as text "{title:Stopping reghdfe without partialling out}" _n
		c_local keep_mata 1
		exit
	}

	if (`verbose' > 0) di as text "{title:Working on varlist: partialling out and regression}" _n

	* Expand varlists
	if (`verbose' > 0) di as text "# Parsing and expanding indepvars: {res}`indepvars'" _c
	if (`timeit') timer on 23
	ms_expand_varlist `indepvars' if `touse'
	if (`timeit') timer off 23
	if (`verbose' > 0) return list
	loc indepvars			"`r(varlist)'"
	loc fullindepvars		"`r(fullvarlist)'"
	loc fullindepvars_bn	"`r(fullvarlist_bn)'"
	loc not_omitted			"`r(not_omitted)'"

	* Partial out variables
	* Syntax:
	* 	partial_out(varnames | ,
	* 		Save TSS if HDFE.tss is missing? [0] ,
	*		Standardize data? [1],  // Note: standardize=2 will standardize, partial out, and return the data standardized!
	*		First col is depvar? [1]
	* 	)
	if (`timeit') timer on 24
	mata: HDFE.partial_out("`depvar' `indepvars'", 1, 1) // saves results in HDFE.solution
	if (`timeit') timer off 24
	// NOTE: I would prefer to do "solution = HDFE.partial_out()" instead of attaching solution to HDFE,
	// but Mata's compiler is unable to do type inference with "class.method()" and can only do it with "function()"

	if (`parallel_maxproc' > 0) {
		if (`timeit') timer on 27
		ParallelBoss // call the parallel processes (workers), wait for them, and finalize by grouping back the data
		if (`timeit') timer off 27
	}

	mata: HDFE.solution.depvar = "`depvar'"
	mata: HDFE.solution.indepvars = tokens("`indepvars'")
	mata: HDFE.solution.fullindepvars = tokens("`fullindepvars'")
	mata: HDFE.solution.fullindepvars_bn = tokens("`fullindepvars_bn'")
	mata: HDFE.solution.indepvar_status = !strtoreal(tokens("1 `not_omitted'")) // indepvar_status[i]=1 for variables omitted due to being a basevar (hack: the first element is the depvar!)
	mata: HDFE.solution.collinear_tol = min(( 1e-6 , HDFE.tolerance / 10))
	mata: HDFE.solution.check_collinear_with_fe(`verbose') // mark variables collinear with the FEs (set HDFE.solution.indepvar_status[i]=2)
	mata: HDFE.solution.fast_regression = `fast_regression'

	if (`stop_before_regression') {
		if (`verbose' > 0) di as text "{title:Stopping reghdfe without running regression}" _n
		c_local keep_mata 1
		exit
	}

	* Regress
	mata: HDFE.solution.report_constant = HDFE.has_intercept & `report_constant' // must set this BEFORE solve_ols()
	if ("`keepmata'" != "") mata: hdfe_data = HDFE.solution.data // uses memory!
	if ("`keepmata'" != "") mata: hdfe_tss = HDFE.solution.tss
	if (`timeit') timer on 25
	mata: reghdfe_solve_ols(HDFE, HDFE.solution, "vce_small")
	if (`timeit') timer off 25
	mata: HDFE.solution.cmdline = HDFE.solution.cmd + " " + st_local("0")

	* Restore
	if (`compact') {
		if (`verbose' > 0) di as text "## Restoring dataset"
		restore
	}

	if (`verbose' > 0) di as text "{title:Posting results to e() and displaying them}" _n
	
	* Post regression results
	* 1) Expand 'b' and 'V' to add base/omitted vars; expand fullindepvars/indepvars to add _cons
	tempname b V
	mata: HDFE.solution.expand_results("`b'", "`V'", HDFE.verbose)

	* 2) Run "ereturn post"
	loc store_sample = ("`nosample'"=="")
	EreturnPost `touse' `b' `V' `store_sample'

	* 4) Store resids if needed (run this before solution.post())
	* Need to save resids if saving FEs, even if temporarily
	mata: st_local("save_any_fe", strofreal(HDFE.save_any_fe))
	if ("`residuals'" == "" & `save_any_fe') loc residuals "__temp_reghdfe_resid__"
	if ("`residuals'" != "") {
		if (`verbose' > 0) di as text "# Storing residuals in {res}`residuals'{txt}" _n
		mata: HDFE.save_variable("`residuals'", HDFE.solution.resid, "Residuals")
		mata: HDFE.solution.residuals_varname = "`residuals'"
	}

	* 3) Post the remaining elements
	mata: HDFE.solution.post()
	mata: HDFE.post_footnote()

	* 4) Store alphas if needed
	if (`timeit') timer on 28
	reghdfe, store_alphas
	if (`timeit') timer off 28

	* 4) View estimation tables
	Replay, `diopts' `header' `table' `footnote'

	if ("`keepmata'" != "") mata: swap(HDFE.solution.data, hdfe_data)
	if ("`keepmata'" != "") mata: swap(HDFE.solution.tss, hdfe_tss)
	if ("`keepmata'" != "") c_local keep_mata 1

	if (`timeit') {
		timer off 20
		ViewTimer, title("Top-level") percent range(20 29) legend(21 "HDFE" 22 "DoF" 23 "Expand factors/Lags" 24 "Partial out" 25 "Solve OLS" 27 "Parallel Boss" 28 "Store alphas")
		ViewTimer, title("Partial-out") percent range(40 49) legend(41 "Load data" 42 "Standardize/etc" 46 "MAP/LSMR" 47 "Data assign" 49 "Parallel Save")
	}
end


program define ParseDOF, sclass
	syntax, [ALL NONE] [FIRSTpair PAIRwise] [CLusters] [CONTinuous]

	* Select all methods by default
	if ("`none'`firstpair'`pairwise'`clusters'`continuous'" == "") local all "all"

	opts_exclusive "`all' `none' `firstpair' `pairwise'" dofadjustments
	opts_exclusive "`all' `none' `clusters'" dofadjustments
	opts_exclusive "`all' `none' `continuous'" dofadjustments

	local opts `pairwise' `firstpair' `clusters' `continuous'
	if ("`all'" != "") local opts pairwise clusters continuous
	sreturn loc dofadjustments "`opts'"
end


program define ValidateGroups, sortpreserve
	syntax varlist, Group_id(varname numeric) TOUSE(name) [INDIVIDUAL(string) INDIVTOUSE(name)]

	* What if I exclude some obs within a group? then the sort order might not work well...
	* EG: junior board members, small neighboring countries, etc.
	* Solution: don't use "if" in syntax, and explictly use -touse- in "by:"

	sort `touse' `group_id' `individual_id'
	
	* 1) Verify that regression variables are constant within group_id
	foreach var of local varlist {
		loc msg "variable `var' is not constant within `group_id'"
		by `touse' `group_id': _assert2 `var' == `var'[1] if `touse', msg("`msg'") rc(498)
	}

	* 2) Validate that individuals are unique within group
	*if ("`individual_id'" != "") {
	*	loc msg "individual identifier {bf:`individual_id'} is not unique within `group_id'"
	*	by `touse' `group_id' `individual_id': _assert2 _N == 1 if `touse', msg("`msg'") error(498)
	*}

	* 2) Even more strict, validate that group and individuals are unique identifiers of the data
	if ("`individual_id'" != "") {
		loc msg "identifiers for group (`group_id') and individual (`individual_id') do not uniquely identify the observations'" 
		by `touse' `group_id' `individual_id': _assert2 _n == 1 if `touse', msg("`msg'") rc(459)
	}

	* 3) Create touse variable that is 1 only once per group
	rename `touse' `indivtouse'
	by `indivtouse' `group_id': gen byte `touse' = (_n == 1) & (`indivtouse') // tag=0 if touse=0
	la var `touse' "[touse]"
	la var `indivtouse' "[touse_individual]"

end


program EreturnPost, eclass
	ereturn clear
	args touse b V store_sample
	mata: st_local("depvar", HDFE.solution.depvar)
	mata: st_local("indepvars", invtokens(HDFE.solution.fullindepvars))
	if (`store_sample') loc esample "esample(`touse')"
	if ("`indepvars'" != "") {
		matrix colnames `b' = `indepvars'
		matrix colnames `V' = `indepvars'
		matrix rownames `V' = `indepvars'
		_ms_findomitted `b' `V'
		ereturn post `b' `V', `esample' buildfvinfo depname(`depvar') 
	}
	else {
		ereturn post, `esample' buildfvinfo depname(`depvar')
	}
end


program define UpdateTouseWithTag
	assert 0 // NOT USED?
	args touse tag group
	tempvar touse_update
	gegen byte `touse_update' = max(`tag'), by(`group_id')
	*tab `touse', m
	qui replace `touse' = 0 if `touse_update' == 0
	*tab `touse', m
end


program Store_Alphas, eclass
	* Explanation:
	* Recall that residuals of regression with dummies are the same as with partialled-out variables
	* Thus, (omitting hats):
	*		e = y - x β - Dα = ytilde - xtilde β
	*		d := Dα = e - (y - xβ)
	* We can compute "d" using known variables (the RHS of above) and then solve the system again

	mata: st_local("save_any_fe", strofreal(HDFE.save_any_fe))
	assert inlist(`save_any_fe', 0, 1)
	if (`save_any_fe') {
		_assert e(depvar) != "", msg("e(depvar) is empty")
		_assert e(resid) != "", msg("e(resid) is empty")
		// we can't use -confirm var- because it might have TS operators
		fvrevar `e(depvar)', list
		confirm numeric var `e(resid)', exact
		tempvar d
		if (e(rank)) {
			qui _predict double `d' if e(sample), xb
		}
		else if (e(report_constant)) {
			gen double `d' = _b[_cons] if e(sample)
		}
		else {
			gen double `d' = 0 if e(sample)
		}
		qui replace `d' = `e(depvar)' - `d' - `e(resid)' if e(sample)

		mata: HDFE.store_alphas("`d'")
		drop `d'

		// Drop resid if we don't want to save it; and update e(resid)
		cap drop __temp_reghdfe_resid__
		if (!c(rc)) ereturn local resid
	}
end


program define ViewTimer, rclass
	syntax, range(numlist min=2 max=2 integer >0 <=100 ascending) LEGend(string asis) [Title(string)] [percent]
	loc show_percent = ("`percent'" != "")
	
	* Process legend
	forval i = 1/100 {
		loc msg`i' "`i':"
	}
	while (`"`legend'"' != "") {
		gettoken key legend : legend
		gettoken val legend : legend
		loc msg`key' `"`val'"'
	}

	* Fill matrix row names based on legend or default values
	qui timer list
	gettoken start end : range
	loc index 0

	* If we show percents we assume that the first value in the range is NOT missing and has the totals
	if (`show_percent') {
		loc total_time = r(t`start') / r(nt`start')
		loc ++start
	}

	forval i = `start'/`end' {
		loc t = r(t`i')
		if (!mi(`t')) {
			loc ++index
			loc t`index' = r(t`i') // / r(nt`i')
			loc rownames `"`rownames' "`msg`i''" "'
		}
	}
	loc num_rows `index'

	* Create and fill matrix
	tempname timer

	if (`show_percent') {
		loc sum_time 0
		matrix `timer' = J(`num_rows'+1, 2, .)
		forval i = 1/`num_rows' {
			loc sum_time = `sum_time' + `t`i''
			matrix `timer'[`i', 1] = `t`i''
			matrix `timer'[`i', 2] = 100 * `t`i'' / `total_time'
		}
		matrix `timer'[`num_rows'+1, 1] = `total_time' - `sum_time'
		matrix `timer'[`num_rows'+1, 2] = 100 * (`total_time' - `sum_time') / `total_time'
		matrix rownames `timer' = `rownames' "(Remainder)"
		matrix colnames `timer' = "Time" "(% Total)"

		di as text _n `"{bf: Timer results:} `title'"'
		loc spaces = (`: rowsof `timer'' - 1) * "&"
		matlist `timer', noblank cspec(& %25s | %10.2fc | %9.1f &) rspec(||`spaces'|) rowtitle(Step)
	}
	else {
		matrix `timer' = J(`num_rows', 1, .)
		forval i = 1/`num_rows' {
			matrix `timer'[`i', 1] = `t`i''
		}
		matrix rownames `timer' = `rownames'
		matrix colnames `timer' = "Time"

		di as text _n `"{bf: Timer results:} `title'"'
		loc spaces = (`: rowsof `timer'' - 1) * "&"
		matlist `timer', noblank cspec(& %25s | %10.2fc &) rspec(||`spaces'|) rowtitle(Step)
	}

	return matrix timer = `timer'
end

program _assert2, byable(onecall)
* Improved version of _assert:
* - allows if/in
* - fails fast
* - can use by:
* TODO: move to -ftools-

        syntax anything(everything) [ , msg(str) rc(str) ]

        if _by() {
        	capture by `_byvars': assert `anything', fast
        }
        else {
        	capture assert `anything', fast
        }

        local rcc = _rc
        if `rcc' {
                if `"`msg'"' != "" {
                        dis as err `"`msg'"'
                }
                else {
                        dis as err `"assert failed: `anything'"'
                }

                if "`rc'" != "" {
                        exit `rc'
                }
                else {
                        exit `rcc'
                }
        }
end


program define ReghdfeFootnote
	reghdfe_footnote
	*di as text "TODO: add info on absorbed individual dummies"
end


program define ParseParallel, sclass
	
	* We intercept three options from -parallel_map-
	syntax anything(name=maxprocesses id="number of worker processes"), ///
		[MAXprocesses2(integer 0) ID(integer 0) TMP_path(string) FORCE] [*]

	* maxprocesses2 -> DISCARDED
	* tmp_path		-> Where will we store the Mata objects
	* id			-> Subfolder within the tmp path where we actually store them

	_assert inrange(`maxprocesses', 0, 1000)
	if (`maxprocesses' == 1) & ("`force'" == "") {
		di as text "(ignoring parallel(1) as it's slower than parallel(0)"
		loc maxprocesses 0 // Launching one extra instance is slower than using the current one
	}

	if (`id' <= 0) {
		loc seed_time = mod(clock("$S_DATE $S_TIME", "DMYhms")/1000, 24*3600)
		loc seed_data = c(N) * c(k)
		loc seed_rand = runiformint(1, 1e9-1)
		loc seed_rand = runiformint(1, 1e9-1)
		loc id = mod((`seed_time' + `seed_data' + `seed_rand'), 1e9-1)
		assert inrange(`id', 1, 1e9-1)
	}

	if ("`tmp_path'" == "") loc tmp_path = c(tmpdir)
	* Workaround for inputs that don't end with "/"
	loc last_char = substr("`tmp_path'", strlen("`tmp_path'"), 1)
	if (!inlist("`last_char'", "/", "\")) loc tmp_path = "`tmp_path'`c(dirsep)'"
	loc padded_caller_id = string(`id', "%09.0f")
	loc parallel_dir = "`tmp_path'PARALLEL_`padded_caller_id'"

	loc options `"maxproc(`maxprocesses') id(`id') tmp_path("`tmp_path'") `options' `force'"'

	loc force_settings = ("`force'" != "")
	sreturn clear
	sreturn loc parallel_force			`force_settings'
	sreturn loc parallel_maxproc 		`maxprocesses'
	sreturn loc parallel_dir 			`"`parallel_dir'"'
	sreturn loc parallel_opts 			`"`options'"'

	* Run the syntax line of parallel_map just to be sure there are no syntax errors
	* This must be done as the end as it overwrites the 'options' local
	loc 0 `", `options'"'
	syntax, [MAXprocesses(integer 0) COREs_per_process(integer 0) FORCE ID(integer 0) ///
			METHOD(string) STATA_path(string) TMP_path(string) PRograms(string) Verbose]	
end


program define ParallelBoss

	mata: st_local("n", strofreal(HDFE.parallel_numproc))
	mata: st_local("opts", HDFE.parallel_opts)
	mata: st_local("parallel_dir", HDFE.parallel_dir)
	mata: st_local("verbose", strofreal(HDFE.verbose))
	
	if (`verbose' <= 0) loc logtable "nologtable"
	if (`verbose' > 0) loc verbose_string "verbose"
	if (`verbose' > 0) di as text "{title:Starting parallel processes:}" _n
	
	loc cmd `"parallel_map, val(1/`n') `verbose_string' `logtable' `opts': reghdfe, worker parallel_path("`parallel_dir'")"'
	if (`verbose' > 0) di as text `"command: {inp}`cmd'"'
	`cmd'

	mata: parallel_combine(HDFE)
end


program ParallelWorker
	syntax, parallel_path(string)
	_assert "${task_id}" != "", msg("global -task_id- cannot be missing")
	_assert (${task_id} == int(${task_id})) & inrange(${task_id}, 1, 1000), msg("global -task_id- must be an integer between 1 and 100")

	loc hdfe_object "`parallel_path'`c(dirsep)'data0.tmp"
	loc vars_object "`parallel_path'`c(dirsep)'data${task_id}.tmp"
	conf file "`hdfe_object'"
	conf file "`vars_object'"

	di as text "Files to load: `hdfe_object'"
	di as text "Files to load: `vars_object'"
	
	mata: worker_partial_out("`hdfe_object'", "`vars_object'")
	di as text "exiting worker thread"
end


include "reghdfe.mata", adopath

exit
