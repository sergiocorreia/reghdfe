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
		Absorb(string)
		NOAbsorb
		RESiduals(name) RESiduals2 /* use _reghdfe_resid */
		SUmmarize SUmmarize2(string asis) /* simulate implicit options */
		SUBOPTions(string) /* gets passed to the e.g regress or ivreg2 */

		/* Standard Errors */
		VCE(string)
		CLuster(string) /* undocumented alternative to vce(cluster ...) */

		/* IV/2SLS/GMM */
		ESTimator(string) /* 2SLS GMM2s CUE LIML */
		STAGEs(string) /* iv (always on) first reduced ols acid (and all) */
		FFirst /* save first-stage stats (only with ivreg2) */
		IVsuite(string) /* ivreg2 ivregress */

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

		/* Speedup and memory Tricks */
		SAVEcache
		USEcache
		CLEARcache
		COMPACT /* use as little memory as possible but is slower */
		NOSAMPle /* do not save e(sample) */

		/* Degrees-of-freedom Adjustments */
		DOFadjustments(string)
		GROUPVar(name) /*var with the first connected group between FEs*/

		/* Undocumented */
		OLD /* use latest v3 */
		KEEPSINgletons
		NOTES(string) /* NOTES(key=value ...), will be stored on e() */
		] [*] /* capture display options, etc. */
		;
	#d cr


* Convert options to boolean
	if ("`verbose'" == "") loc verbose 0
	mata: HDFE.verbose = `verbose'
	mata: HDFE.timeit = ("`timeit'"!="")
	mata: `OPT'.drop_singletons = ("`keepsingletons'"=="")
	loc ffirst = ("`ffirst'"!="")
	mata: `OPT'.ffirst = `ffirst'

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
	mata: `OPT'.fe_format = "`s(fe_format)'"
	loc model = cond("`s(instruments)'" == "", "ols", "iv")
	mata: `OPT'.model = "`model'"
	mata: `OPT'.original_varlist = "`s(varlist)'" // as before but without parens or equal
	// TODO:
	// fvexpand `indepvars'
	// local cnames `r(varlist)'


* Parse Estimator (picks the estimation subcommand)
	if ("`model'" == "iv") {
		if ("`ivsuite'"=="") local ivsuite ivreg2 // default
		_assert inlist("`ivsuite'","ivreg2","ivregress") , msg("error: wrong IV routine (`ivsuite'), valid options are -ivreg2- and -ivregress-")
		cap findfile `ivsuite'.ado
		_assert !_rc , msg("error: -`ivsuite'- not installed, please run {stata ssc install `ivsuite'} or change the ivsuite() option")
		local subcmd `ivsuite'

		if ("`estimator'"=="") local estimator 2sls // Set default
		if (substr("`estimator'", 1, 3)=="gmm") local estimator gmm2s
		_assert inlist("`estimator'", "2sls", "gmm2s", "liml", "cue"), msg("reghdfe error: invalid estimator `estimator'")
		if ("`estimator'"=="cue") _assert ("`ivsuite'"=="ivreg2"), msg("reghdfe error: estimator `estimator' only available with the ivreg2 command, not ivregress")
		if ("`estimator'"=="cue") di as error "(WARNING: -cue- estimator is not exact, see help file)"
	}
	else {
		local subcmd regress
		_assert "`estimator'"=="", msg("estimator() requires an instrumental-variable regression")
	}
	mata: `OPT'.estimator = "`estimator'"
	mata: `OPT'.ivsuite = "`ivsuite'"
	mata: `OUT'.subcmd = "`subcmd'"


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


* Parse stages
	ParseStages, stages(`stages') model("`model'")
	mata: `OPT'.stages = "`s(stages)'"
	mata: `OPT'.stages_save = `s(savestages)'
	mata: `OPT'.stages_opt = "`s(stage_suboptions)'"


* Parse VCE
	ms_parse_vce, model(`model') vce(`vce') weighttype(`weight_type') ivsuite(`ivsuite')
	//sreturn list
	mata: `OPT'.vceoption = "`s(vceoption)'"
	mata: `OPT'.vcetype = "`s(vcetype)'"
	mata: `OPT'.vcesuite = "`s(vcesuite)'"
	mata: `OPT'.vceextra = "`s(vceextra)'"
	mata: `OPT'.vce_is_hac = `s(vce_is_hac)'
	mata: `OPT'.num_clusters = `s(num_clusters)'
	mata: `OPT'.clustervars = tokens("`s(clustervars)'")
	mata: `OPT'.base_clustervars = tokens("`s(base_clustervars)'")
	mata: `OPT'.bw = `s(bw)'
	mata: `OPT'.kernel = "`s(kernel)'"
	mata: `OPT'.dkraay = `s(dkraay)'
	mata: `OPT'.twicerobust = `s(twicerobust)'
	mata: `OPT'.kiefer = "`s(kiefer)'"
	loc base_clustervars "`s(base_clustervars)'"


* Sanity checks for -ffirst-
	if (`ffirst') {
		_assert ("`model'" != "ols"), msg("ols does not support {cmd}ffirst")
		_assert ("`ivsuite'" == "ivreg2"), msg("option {bf:ffirst} requires ivreg2")
	}


* DoF Adjustments
	if ("`dofadjustments'"=="") local dofadjustments all
	ParseDOF , `dofadjustments'
	if ("`groupvar'"!="") conf new var `groupvar'
	mata: `OPT'.dofadjustments = "`s(dofadjustments)'"
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
	_assert (`"`options'"'==""), ///
		msg(`"invalid options: `options'"')
	if ("`hascons'"!="") di in ye "(option ignored: `hascons')"
	if ("`tsscons'"!="") di in ye "(option ignored: `tsscons')"
	mata: `OPT'.diopts = `"`diopts'"'


* Mark sample
	c_local varlist `"`depvar' `indepvars' `endogvars' `instruments' `extended_absvars' `base_clustervars'"'
	c_local if `"`if'"'
	c_local in `"`in'"'
	c_local weight `"`weight'"'
	c_local exp `"`exp'"'

	// if (`verbose' > 0) ViewMataOptions
end



cap pr drop ParseSummarize
pr ParseSummarize, sclass
	sreturn clear
	syntax [namelist(name=stats)] , [QUIetly]
	local quietly = ("`quietly'"!="")
	sreturn loc stats "`stats'"
	sreturn loc quietly = `quietly'
end


cap pr drop ParseStages
pr ParseStages, sclass
	syntax, model(string) [stages(string)]
	local 0 `stages'
	syntax [namelist(name=stages)], [noSAVE] [*]
	
	if ("`stages'"=="") local stages none
	if ("`stages'"=="all") local stages iv first ols reduced acid

	if ("`stages'"!="none") {
		_assert ("`model'" == "iv"), ///
			msg("{cmd:stages(`stages')} not allowed with ols")
		local special iv none
		local valid_stages first ols reduced acid
		local stages : list stages - special
		local wrong_stages : list stages - valid_stages
		_assert "`wrong_stages'"=="", msg("Error, invalid stages(): `wrong_stages'")
		* The "iv" stage will be always on for IV-type regressions
		local stages `stages' iv // put it last so it does the restore
	}

	sreturn local stages `stages'
	sreturn local stage_suboptions `options'
	sreturn local savestages = ("`save'"!="nosave" & "`stages'"!="none")
end


cap pr drop ParseDOF
pr ParseDOF, sclass
	sreturn clear
	syntax, [ALL NONE] [PAIRwise TWO THREE] [CLusters] [CONTinuous]
	local opts `pairwise' `two' `three' `clusters' `continuous'
	local n : word count `opts'
	local first_opt : word 1 of `opt'

	opts_exclusive "`all' `none'" dofadjustments
	opts_exclusive "`pairwise' `two' `three'" dofadjustments
	opts_exclusive "`all' `first_opt'" dofadjustments
	opts_exclusive "`none' `first_opt'" dofadjustments

	if ("`none'" != "") local opts
	if ("`all'" != "") local opts pairwise clusters continuous

	if (`: list posof "three" in opts') {
		cap findfile group3hdfe.ado
		Assert !_rc , msg("error: -group3hdfe- not installed, please run {stata ssc install group3hdfe}")
	}

	sreturn local dofadjustments "`opts'"
end
