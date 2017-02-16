program reghdfe_parse
	loc OPT = "HDFE.options"
	loc OUT = "HDFE.output"

* Trim whitespace (caused by "///" line continuations; aesthetic only)
	mata: st_local("0", stritrim(st_local("0")))
	mata: `OUT'.cmdline = "reghdfe " + st_local("0")

* Main syntax
	#d;
	syntax anything(id=varlist equalok) [if] [in] [aw pw fw/] , [

		/* Model */
		Absorb(string) NOAbsorb
		RESiduals(name) RESiduals2 /* use _reghdfe_resid */
		SUmmarize SUmmarize2(string asis) /* simulate implicit options */

		/* Standard Errors */
		VCE(string) CLuster(string)

		/* Diagnostic */
		Verbose(numlist min=1 max=1 >=-1 <=5 integer)
		TIMEit

		/* Optimization (defaults are handled within Mata) */
		TOLerance(real -1)
		MAXITerations(real -1)
		ALGorithm(string) /* map gt lsmr cg */
		TRAnsform(string)
		ACCELeration(string)
		SLOPEmethod(string)
		noPRUNE
		KEEPSINgletons

		/* Degrees-of-freedom Adjustments */
		DOFadjustments(string)
		GROUPVar(name) /* var with the first connected group between FEs */

		/* Undocumented */
		OLD /* use latest v3 */
		NOTES(string) /* NOTES(key=value ...), will be stored on e() */
		] [*] /* capture display options, etc. */
		;
	#d cr

* Unused
		/* Speedup and memory Tricks */
*		SAVEcache
*		USEcache
*		CLEARcache
*		COMPACT /* use as little memory as possible but is slower */
*		NOSAMPle /* do not save e(sample) */

* Convert options to boolean
	if ("`verbose'" == "") loc verbose 0
	mata: HDFE.verbose = `verbose'
	mata: HDFE.timeit = ("`timeit'"!="")
	mata: `OPT'.drop_singletons = ("`keepsingletons'"=="")

	if (`verbose'>-1 & "`keepsingletons'"!="") {
		di as error `"WARNING: Singleton observations not dropped; statistical significance is biased {browse "http://scorreia.com/reghdfe/nested_within_cluster.pdf":(link)}"'
	}

* Sanity checks
	if ("`cluster'"!="") {
		_assert ("`vce'"==""), msg("cannot specify both cluster() and vce()")
		loc vce cluster `cluster'
		loc cluster // clear it to avoid bugs in subsequent lines
	}


* Optimization
	loc maxiterations = int(`maxiterations')
	if (`tolerance' > 0) mata: HDFE.tolerance = `tolerance'
	if (`maxiterations' > 0) mata: HDFE.maxiter = `maxiterations'

	* Transforms: allow abbreviations (cim --> cimmino)
	if ("`transform'" != "") {
		loc transform = lower("`transform'")
		loc valid_transforms cimmino kaczmarz symmetric_kaczmarz rand_kaczmarz
		foreach x of local valid_transforms {
			if (strpos("`x'", "`transform'")==1) loc transform `x'
		}
		_assert (`: list transform in valid_transforms'), msg("invalid transform: `transform'")
		mata: HDFE.transform = "`transform'"
	}

	* Accelerations
	if ("`acceleration'" != "") {
		loc acceleration = lower("`acceleration'")
		if ("`acceleration'"=="cg") loc acceleration conjugate_gradient
		if ("`acceleration'"=="sd") loc acceleration steepest_descent
		if ("`acceleration'"=="off") loc acceleration none
		loc valid_accelerations conjugate_gradient steepest_descent aitken none hybrid
		foreach x of local valid_accelerations {
			if (strpos("`x'", "`acceleration'")==1) loc acceleration `x'
		}
		_assert (`: list acceleration in valid_accelerations'), msg("invalid acceleration: `acceleration'")
		mata: HDFE.acceleration = "`acceleration'"
	}

	* Disable prune of degree-1 edges
	if ("`prune'"=="noprune") mata: HDFE.prune = 0


* Parse Varlist
	ms_fvunab `anything'
	// mata: HDFE.base_varlist = "`s(basevars)'"
	ms_parse_varlist `s(varlist)'
	foreach cat in depvar indepvars endogvars instruments {
		loc `cat' "`s(`cat')'"
		mata: `OPT'.`cat' = `OPT'.original_`cat' = "``cat''"
	}
	loc model = cond("`s(instruments)'" == "", "ols", "iv")
	mata: `OPT'.model = "`model'"
	mata: `OPT'.original_varlist = "`s(varlist)'" // as before but without parens or equal


* Parse Weights
	if ("`weight'"!="") {
		unab exp : `exp', min(1) max(1) // simple weights only
		mata: `OPT'.weight_type = "`weight'"
		mata: `OPT'.weight_var = "`exp'"
		mata: `OPT'.weight_exp = `"[`weight'=`exp']"'
	}


* Parse Absvars (only used for validation, run again from Mata!)
	_assert  ("`absorb'`noabsorb'" != ""), msg("option {bf:absorb()} or {bf:noabsorb} required")
	if ("`noabsorb'" != "") {
		_assert ("`absorb'" == ""), msg("{bf:absorb()} and {bf:noabsorb} are mutually exclusive")
		tempvar c
		gen byte `c' = 1
		loc absorb `c'
	}
	ms_parse_absvars `absorb'
	loc extended_absvars `"`s(extended_absvars)'"'
	mata: `OPT'.absorb = `"`absorb'"'


* Parse summarize
	if ("`summarize'" != "") {
		_assert ("`summarize2'" == ""), msg("summarize() syntax error")
		loc summarize2 mean min max  // default values
	}
	ParseSummarize `summarize2'
	mata: `OPT'.summarize_stats = "`s(stats)'"
	mata: `OPT'.summarize_quietly = `s(quietly)'


* Parse VCE
	ms_parse_vce, vce(`vce') weighttype(`weight_type')
	mata: `OPT'.vcetype = "`s(vcetype)'"
	mata: `OPT'.num_clusters = `s(num_clusters)'
	mata: `OPT'.clustervars = tokens("`s(clustervars)'")
	mata: `OPT'.base_clustervars = tokens("`s(base_clustervars)'")
	mata: `OPT'.vceextra = "`s(vceextra)'"
	loc base_clustervars "`s(base_clustervars)'"


* Parse DoF Adjustments
	if ("`dofadjustments'"=="") local dofadjustments all
	ParseDOF , `dofadjustments'
	if ("`groupvar'"!="") conf new var `groupvar'
	mata: `OPT'.dofadjustments = tokens("`s(dofadjustments)'")
	mata: `OPT'.groupvar = "`s(groupvar)'"


* Parse residuals
	if ("`residuals2'" != "") {
		_assert ("`residuals'" == ""), msg("residuals() syntax error")
		loc residuals _reghdfe_resid
		cap drop `residuals' // destructive!
	}
	else if ("`residuals'"!="") {
		conf new var `residuals'
	}
	mata: `OPT'.residuals = "`residuals'"


* Parse misc options
	mata: `OUT'.notes = `"`notes'"'
	mata: `OPT'.suboptions = `"`suboptions'"'


* Parse Coef Table Options (do this last!)
	_get_diopts diopts options, `options' // store in `diopts', and the rest back to `options'
	_assert (`"`options'"'==""), msg(`"invalid options: `options'"')
	if ("`hascons'"!="") di in ye "(option ignored: `hascons')"
	if ("`tsscons'"!="") di in ye "(option ignored: `tsscons')"
	mata: `OPT'.diopts = `"`diopts'"'


* Inject back locals
	c_local varlist `"`depvar' `indepvars' `endogvars' `instruments' `extended_absvars' `base_clustervars'"'
	c_local if `"`if'"'
	c_local in `"`in'"'
	c_local weight `"`weight'"'
	c_local exp `"`exp'"'
end


program ParseSummarize, sclass
	sreturn clear
	syntax [namelist(name=stats)] , [QUIetly]
	local quietly = ("`quietly'"!="")
	sreturn loc stats "`stats'"
	sreturn loc quietly = `quietly'
end


program ParseDOF, sclass
	sreturn clear
	syntax, [ALL NONE] [FIRSTpair PAIRwise] [CLusters] [CONTinuous]
	local opts `pairwise' `firstpair' `clusters' `continuous'
	local n : word count `opts'
	local first_opt : word 1 of `opt'

	opts_exclusive "`all' `none'" dofadjustments
	opts_exclusive "`pairwise' `firstpair'" dofadjustments
	opts_exclusive "`all' `first_opt'" dofadjustments
	opts_exclusive "`none' `first_opt'" dofadjustments

	if ("`none'" != "") local opts
	if ("`all'" != "") local opts pairwise clusters continuous

	//if (`: list posof "three" in opts') {
	//	cap findfile group3hdfe.ado
	//	_assert !_rc , msg("error: -group3hdfe- not installed, please run {stata ssc install group3hdfe}")
	//}

	sreturn local dofadjustments "`opts'"
end
