*! version 4.0.5 dev - 16feb2017

program reghdfe, eclass
	* Intercept old+version
	cap syntax, version old
	if !c(rc) {
		reghdfe_old, version
		exit
	}

	* Intercept old
	cap syntax anything(everything) [fw aw pw/], [*] old
	if !c(rc) {
		di as error "(running historical version of reghdfe)"
		if ("`weight'"!="") local weightexp [`weight'=`exp']
		reghdfe_old `anything' `weightexp', `options'
		exit
	}

	* Aux. subcommands
	cap syntax, [*]
	if inlist("`options'", "check", "compile", "reload", "update", "version") {
		if ("`options'"=="compile") loc args force
		if ("`options'"=="check") loc options compile
		if ("`options'"=="update") {
			loc args 1
			loc options reload
		}
		loc subcmd = proper("`options'")
		`subcmd' `args'
	}
	else if replay() {
		Replay `0'
	}
	else {
		Cleanup 0
		Compile // takes 0.01s to run this useful check (ensures .mlib exists)
		cap noi Estimate `0'
		Cleanup `c(rc)'
	}
end

// --------------------------------------------------------------------------

program Compile
	args force
	ftools, check // in case lftools.mlib does not exist or is outdated
	ms_get_version reghdfe // from moresyntax package; save local package_version

	loc list_objects "FixedEffects() fixed_effects() FE_Options() FE_Output() BipartiteGraph()"
	loc list_functions "reghdfe_*() transform_*() accelerate_*() panelmean() panelsolve_*()"
	loc list_misc "weighted_quadcolsum() safe_divide() check_convergence()"
	// TODO: prefix everything with reghdfe_*

	ms_compile_mata, ///
		package(reghdfe) ///
		version(`package_version') ///
		fun("`list_objects' `list_functions' `list_misc'") ///
		verbose ///
		`force'
end

// --------------------------------------------------------------------------

program Reload
	* Internal debugging tool.
	* Updates dependencies and reghdfe from local path or from github
	* Usage:
	* 	reghdfe, update // from c:\git\..
	* 	reghdfe, reload // from github

	args online
	if ("`online'" == "") loc online 0

	di as text _n "{bf:reghdfe: updating required packages}"
	di as text "{hline 64}"

	* -moresyntax- https://github.com/sergiocorreia/moresyntax/
	cap ado uninstall moresyntax
	if (`online') net install moresyntax, from("https://github.com/sergiocorreia/moresyntax/raw/master/src/")
	if (!`online') net install moresyntax, from("c:\git\moresyntax\src")
	di as text "{hline 64}"

	* -ftools- https://github.com/sergiocorreia/ftools/
	cap ado uninstall ftools
	if (`online') net install ftools, from("https://github.com/sergiocorreia/ftools/raw/master/src/")
	if (!`online') net install ftools, from("c:\git\ftools\src")
	di as text "{hline 64}"
	ftools, compile // requires moresyntax
	di as text "{hline 64}"

	* Update -reghdfe-
	di as text _n  _n "{bf:reghdfe: updating self}"
	di as text "{hline 64}"
	qui ado uninstall reghdfe
	if (`online') net install reghdfe, from("https://github.com/sergiocorreia/reghdfe/raw/version-4/src/")
	if (!`online') net install reghdfe, from("c:\git\reghdfe\src")
	qui which reghdfe
	di as text "{hline 64}"
	reghdfe, compile
	di as text "{hline 64}"

	* Cleaning up
	di as text _n "{bf:Note:} You need to run {stata program drop _all} now."
end

// --------------------------------------------------------------------------

program Version
	which reghdfe

	di as text _n "Dependencies installed?"
	local dependencies moresyntax ftools
	// ivreg2 avar tuples group3hdfe boottest
	foreach dependency of local dependencies {
		loc fn `dependency'.ado
		if ("`dependency'"=="moresyntax") loc fn ms_get_version.ado
		cap findfile `fn'
		if (_rc) {
			di as text "{lalign 20:- `dependency'}" as error "not"
		}
		else {
			di as text "{lalign 20:- `dependency'}" as result "yes"
		}
	}
end

// --------------------------------------------------------------------------

program Cleanup
	args rc
	cap mata: mata drop HDFE
	cap mata: mata drop hdfe_*
	cap matrix drop reghdfe_statsmatrix
	if (`rc') exit `rc'
end

// --------------------------------------------------------------------------

program Parse
	loc OPT = "HDFE.options"
	loc OUT = "HDFE.output"

* Trim whitespace (caused by "///" line continuations; aesthetic only)
	mata: st_local("0", stritrim(st_local("0")))

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

		/* Undocumented */
		KEEPSINgletons
		OLD /* use latest v3 */
		NOTES(string) /* NOTES(key=value ...), will be stored on e() */
		
		] [*] /* capture optimization options, display options, etc. */
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
	loc timeit = ("`timeit'"!="")
	loc drop_singletons = ("`keepsingletons'"=="")

	if (`verbose'>-1 & "`keepsingletons'"!="") {
		di as error `"WARNING: Singleton observations not dropped; statistical significance is biased {browse "http://scorreia.com/reghdfe/nested_within_cluster.pdf":(link)}"'
	}


* Sanity checks
	if ("`cluster'"!="") {
		_assert ("`vce'"==""), msg("cannot specify both cluster() and vce()")
		loc vce cluster `cluster'
		loc cluster // clear it to avoid bugs in subsequent lines
	}


* Parse Varlist
	ms_fvunab `anything'
	// mata: HDFE.base_varlist = "`s(basevars)'"
	ms_parse_varlist `s(varlist)'
	foreach cat in depvar indepvars endogvars instruments {
		loc `cat' "`s(`cat')'"
	}
	loc model = cond("`s(instruments)'" == "", "ols", "iv")
	//mata: `OPT'.model = "`model'"
	loc original_varlist = "`s(varlist)'" // no parens or equal // bugbug add to HDFE.options


* Parse Weights
	if ("`weight'"!="") {
		unab exp : `exp', min(1) max(1) // simple weights only
	}


* Parse VCE
	ms_parse_vce, vce(`vce') weighttype(`weight_type')
	loc vcetype = "`s(vcetype)'"
	loc num_clusters = `s(num_clusters)'
	loc clustervars = "`s(clustervars)'"
	loc base_clustervars = "`s(base_clustervars)'"
	loc vceextra = "`s(vceextra)'"


* Select sample
	loc varlist `original_varlist' `base_clustervars'
	tempvar touse
	marksample touse, strok // based on varlist + cluster + if + in + weight


* Parse noabsorb (only used for validation, run again from Mata!)
	_assert  ("`absorb'`noabsorb'" != ""), msg("option {bf:absorb()} or {bf:noabsorb} required")
	if ("`noabsorb'" != "") {
		_assert ("`absorb'" == ""), msg("{bf:absorb()} and {bf:noabsorb} are mutually exclusive")
		tempvar c
		gen byte `c' = 1
		loc absorb `c'
	}


* Construct HDFE object
	// SYNTAX: fixed_effects(absvars | , touse, weighttype, weightvar, dropsing, verbose)
	mata: st_local("comma", strpos(`"`absorb'"', ",") ? "" : ",")
	mata: HDFE = fixed_effects(`"`absorb' `comma' `options'"', "`touse'", "`weight'", "`exp'", `drop_singletons', `verbose')
	mata: `OUT'.cmdline = "reghdfe " + st_local("0")
	loc options `s(options)'

	mata: st_local("N", strofreal(HDFE.N))
	if (`N' == 0) error 2000

* Fill out HDFE object
	mata: `OPT'.clustervars = tokens("`clustervars'")
	...


* Parse summarize
	if ("`summarize'" != "") {
		_assert ("`summarize2'" == ""), msg("summarize() syntax error")
		loc summarize2 mean min max  // default values
	}
	ParseSummarize `summarize2'
	mata: `OPT'.summarize_stats = "`s(stats)'"
	mata: `OPT'.summarize_quietly = `s(quietly)'


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

// --------------------------------------------------------------------------

program ParseSummarize, sclass
	sreturn clear
	syntax [namelist(name=stats)] , [QUIetly]
	local quietly = ("`quietly'"!="")
	sreturn loc stats "`stats'"
	sreturn loc quietly = `quietly'
end

// --------------------------------------------------------------------------

program Estimate, eclass
	ereturn clear

* Parse and fill out HDFE object
	Parse `0'
	loc OPT HDFE.options
	loc OUT HDFE.output


* Compute degrees-of-freedom
	mata: HDFE.estimate_dof()

* Save updated e(sample) (singletons reduce sample);
* required to parse factor variables to partial out
	mata: HDFE.save_touse("`touse'", 1)

* Expand varlists
	foreach cat in varlist depvar indepvars endogvars instruments {
		mata: st_local("vars", `OPT'.original_`cat')
		if ("`vars'" == "") continue
		// HACK: addbn replaces 0.foreign with 0bn.foreign , to prevent st_data() from loading a bunch of zeros
		ms_fvstrip `vars' if `touse', expand dropomit addbn onebyone
		// If we don't use onebyone, then 1.x 2.x ends up as 2.x
		loc vars "`r(varlist)'"
		mata: `OPT'.`cat' = "`vars'"
	}

* Stats
	mata: st_local("stats", `OPT'.summarize_stats)
	if ("`stats'" != "") Stats

* Partial out; and save TSS of depvar
	mata: hdfe_variables = HDFE.partial_out(`OPT'.varlist, 1) // 1=Save TSS of first var if HDFE.output.tss is missing

* Regress
	mata: assert(`OPT'.model=="ols")
	Regress `touse'
	mata: st_local("diopts", `OPT'.diopts)
	Replay, `diopts'


// ~~ Preserve relevant dataset ~~
	// 
	// if (`compact') {
	// 	preserve
	// 	keep `varlist' `absorb'
	// 	loc N = c(N)
	// }
	// else {
	// 	qui cou if `touse'
	// 	loc N = r(N)
	// }
end


// -------------------------------------------------------------
// Run OLS regressions with partialled-out variables
// -------------------------------------------------------------
program Regress, eclass
	args touse

	tempname b V N rank df_r
	// Carefully split y from X, without duplicating the memory used
	mata: hdfe_y = hdfe_variables[., 1]
	mata: hdfe_variables = cols(hdfe_variables)==1 ? J(rows(hdfe_variables), 0, .) : hdfe_variables[., 2..cols(hdfe_variables)]
	mata: reghdfe_post_ols(HDFE, hdfe_y, hdfe_variables, "`b'", "`V'", "`N'", "`rank'", "`df_r'")
	mata: st_local("indepvars", HDFE.options.indepvars)
	mata: hdfe_y = hdfe_variables = .
	
	if ("`indepvars'" != "") {
		matrix colnames `b' = `indepvars'
		matrix colnames `V' = `indepvars'
		matrix rownames `V' = `indepvars'
		ereturn post `b' `V', esample(`touse') buildfvinfo depname(`depvar') 
	}
	else {
		ereturn post, esample(`touse') buildfvinfo depname(`depvar')
	}

	ereturn scalar N       = `N'
	ereturn scalar rank    = `rank'
	ereturn scalar df_r    = `df_r'
	ereturn local  cmd     "reghdfe"
	mata: HDFE.output.post(HDFE.options)

	* Post stats
	cap conf matrix reghdfe_statsmatrix
	if (!c(rc)) {
		ereturn matrix summarize = reghdfe_statsmatrix
		mata: st_local("summarize_quietly", strofreal(HDFE.options.summarize_quietly))
		ereturn scalar summarize_quietly = `summarize_quietly'
	}
end



program Replay
	syntax [, *]
	_get_diopts options, `options'
	reghdfe_header // _coef_table_header
	di ""
	_coef_table, `options' // ereturn display, `options'
	reghdfe_footnote

	* Replay stats
	if (e(summarize_quietly)==0) {
		di as text _n "{sf:Regression Summary Statistics:}" _c
		matlist e(summarize)', border(top bottom) rowtitle(Variable) // twidth(18) 
	}
end



// -----------------------------------------------------------------------------
// Matrix of summary statistics
// -----------------------------------------------------------------------------
program Stats
	* Optional weights
	mata: st_local("weight", sprintf("[%s=%s]", HDFE.weighttype, HDFE.weightvar))
	assert "`weight'" != ""
	if ("`weight'" == "[=]") loc weight
	loc weight : subinstr local weight "[pweight" "[aweight"

	mata: st_local("stats", HDFE.options.summarize_stats)
	mata: st_local("varlist", HDFE.options.varlist)
	mata: st_local("cvars", invtokens(HDFE.cvars))
	loc full_varlist `varlist' `cvars'

	qui tabstat `full_varlist' `weight' , stat(`stats') col(stat) save
	matrix reghdfe_statsmatrix = r(StatTotal)
end
