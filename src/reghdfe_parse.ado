* This program should only be called by fixed_effects()
program reghdfe_parse, sclass

* Parse absorb
	cap drop __hdfe* // destructive!
	ms_parse_absvars `0'
	loc extended_absvars `"`s(extended_absvars)'"'
	mata: st_local("unquoted_absvars", subinstr(st_global("s(absvars)"), `"""', ""))
	loc 0, `s(options)'
	loc G = `s(G)'

* Main syntax
	#d;
	syntax, [

		/* Model */
		RESiduals(name) RESiduals2 /* use _reghdfe_resid */

		/* Optimization (defaults are handled within Mata) */
		TOLerance(real -1)
		MAXITerations(real -1)
		ALGorithm(string) /* map gt lsmr cg */
		TRAnsform(string)
		ACCELeration(string)
		SLOPEmethod(string)
		PRUNE
		PRECONDition /* always compute LSMR preconditioner */

		/* Memory usage (also see -compact- option) */ 
		POOLsize(integer 0) /* Process variables in batches of # ; 0 turns it off */

		/* Degrees-of-freedom Adjustments */
		DOFadjustments(string)
		GROUPVar(name) /* var with the first connected group between FEs */

		CONDition // Report finite condition number; SLOW!
		RRE(varname) // Report relative residual error
		noCONstant // Report constant; enabled by default as otherwise -margins- fails

		/* Duplicated options */
		KEEPSINgletons
		Verbose(numlist min=1 max=1 >=-1 <=5 integer)

		] [*] /* capture display options, etc. */
		;
	#d cr

	if ("`keepsingletons'"!="") sreturn loc drop_singletons = 0
	if ("`verbose'"!="") sreturn loc verbose = `verbose'
	sreturn loc report_constant = "`constant'" != "noconstant"

	sreturn loc options `"`options'"'

	assert "$reghdfe_touse" != ""
	cap conf var $reghdfe_touse
	if (c(rc)) gen byte $reghdfe_touse = 1
	markout $reghdfe_touse `unquoted_absvars', strok

* Optimization
	loc maxiterations = int(`maxiterations')
	if (`tolerance' > 0) sreturn loc tolerance = `tolerance'
	if (`maxiterations' > 0) sreturn loc maxiter = `maxiterations'

	* Transforms: allow abbreviations (cim --> cimmino)
	if ("`transform'" != "") {
		loc transform = lower("`transform'")
		loc valid_transforms cimmino kaczmarz symmetric_kaczmarz rand_kaczmarz
		foreach x of local valid_transforms {
			if (strpos("`x'", "`transform'")==1) loc transform `x'
		}
		_assert (`: list transform in valid_transforms'), msg("invalid transform: `transform'")
		sreturn loc transform "`transform'"
	}

	* Accelerations
	if ("`acceleration'" != "") {
		loc acceleration = lower("`acceleration'")
		if ("`acceleration'"=="cg") loc acceleration conjugate_gradient
		if ("`acceleration'"=="sd") loc acceleration steepest_descent
		if ("`acceleration'"=="off") loc acceleration none
		loc valid_accelerations conjugate_gradient steepest_descent aitken none hybrid lsmr
		foreach x of local valid_accelerations {
			if (strpos("`x'", "`acceleration'")==1) loc acceleration `x'
		}
		_assert (`: list acceleration in valid_accelerations'), msg("invalid acceleration: `acceleration'")
		sreturn loc acceleration "`acceleration'"
	}

	* Disable prune of degree-1 edges
	if ("`prune'" == "prune") sreturn loc prune = 1

* Parse DoF Adjustments
	if ("`dofadjustments'"=="") local dofadjustments all
	loc 0 , `dofadjustments'
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
	if ("`groupvar'"!="") conf new var `groupvar'
	sreturn local dofadjustments "`opts'"
	sreturn loc groupvar "`s(groupvar)'"

* Residuals
	if ("`residuals2'" != "") {
		_assert ("`residuals'" == ""), msg("residuals() syntax error")
		cap drop _reghdfe_resid // destructive!
		sreturn loc residuals _reghdfe_resid
	}
	else if ("`residuals'"!="") {
		conf new var `residuals'
		sreturn loc residuals `residuals'
	}

* Misc
	if ("`condition'"!="") {
		_assert `G'==2, msg("Computing finite condition number requires two FEs")
		sreturn loc finite_condition 1
	}

	sreturn loc compute_rre = ("`rre'" != "")
	if ("`rre'" != "") {
		sreturn loc rre `rre'
	}

	if (`poolsize' < 1) loc poolsize .
	sreturn loc poolsize `poolsize'

	sreturn loc precondition = "`precondition'" != ""
end
