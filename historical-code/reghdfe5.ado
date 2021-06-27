*! version 5.9.0 03jun2020

program reghdfe5, eclass

	* Aux. subcommands
	cap syntax, [*]

	//if "$fake"!="0" & !inlist("`options'", "check", "compile", "reload", "update", "version", "requirements") {
	//	groupreg `0'
	//	exit
	//}

	if inlist("`options'", "check", "compile", "reload", "update", "version", "requirements", "store_alphas") {
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
		ms_get_version ftools, min_version("2.39.1")
		cap noi Estimate `0'
		Cleanup `c(rc)'
	}
end


program Compile
	_assert 0, msg("reghdfe5 compile disabled (is it needed?)")
	args force
	
	* Check dependencies
	ftools, check // in case lftools.mlib does not exist or is outdated
	ms_get_version ftools, min_version("2.39.1")
	ms_get_version reghdfe // save local package_version
	loc list_objects "FixedEffects() fixed_effects() BipartiteGraph()"
	loc list_functions "reghdfe_*() transform_*() accelerate_*() panelmean() panelsolve_*() lsmr()"
	loc list_misc "weighted_quadcolsum() safe_divide() check_convergence() precompute_inv_xx() _st_data_wrapper()"
	// TODO: prefix everything with reghdfe_*

	ms_compile_mata, ///
		package(reghdfe) ///
		version(`package_version') ///
		fun("`list_objects' `list_functions' `list_misc'") ///
		verbose ///
		`force'
end


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

	* -ftools- https://github.com/sergiocorreia/ftools/
	cap ado uninstall ftools
	if (`online') net install ftools, from("https://github.com/sergiocorreia/ftools/raw/master/src/")
	if (!`online') net install ftools, from("c:\git\ftools\src")
	di as text "{hline 64}"
	ftools, compile
	di as text "{hline 64}"

	* Update -reghdfe-
	di as text _n  _n "{bf:reghdfe: updating self}"
	di as text "{hline 64}"
	qui ado uninstall reghdfe
	if (`online') net install reghdfe, from("https://github.com/sergiocorreia/reghdfe/raw/master/src/")
	if (!`online') net install reghdfe, from("c:\git\reghdfe\src")
	qui which reghdfe
	di as text "{hline 64}"
	reghdfe, compile
	di as text "{hline 64}"

	* Cleaning up
	di as text _n "{bf:Note:} You need to run {stata program drop _all} now."
end


program Version
	which reghdfe
	Requirements
end


program Requirements
	di as text _n "Required packages installed?"
	loc reqs ftools
	// ivreg2 avar tuples group3hdfe
	if (c(stata_version)<13) loc reqs `reqs' boottest

	loc ftools_github "https://github.com/sergiocorreia/ftools/raw/master/src/"

	loc error 0

	foreach req of local reqs {
		loc fn `req'.ado
		cap findfile `fn'
		if (_rc) {
			loc error 1
			di as text "{lalign 20:- `req'}" as error "not" _c
			di as text "    {stata ssc install `req':install from SSC}" _c
			if inlist("`req'", "ftools") {
				loc github ``req'_github'
				di as text `"    {stata `"net install `req', from(`"`github'"')"':install from github}"'
			}
			else {
				di as text // newline
			}
		}
		else {
			di as text "{lalign 20:- `req'}" as text "yes"
		}
	}

	if (`error') exit 601
end


program Store_Alphas, eclass
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


program Cleanup
	args rc
	cap mata: mata drop HDFE
	cap mata: mata drop hdfe_*
	cap drop __temp_reghdfe_resid__
	cap matrix drop reghdfe_statsmatrix
	if (`rc' == 132) {
		di as text "- If you got the {it:parentheses unbalanced} error, note that IV/2SLS was moved to {help ivreghdfe}"
		di as smcl `"- Latest version: {browse "https://github.com/sergiocorreia/ivreghdfe":https://github.com/sergiocorreia/ivreghdfe}"'
		di as smcl `"- SSC version: {stata "net describe ivreghdfe, from(http://fmwww.bc.edu/RePEc/bocode/i)"}"'
		di as smcl `"- Note: the older functionality can still be accessed through the {it:old} option"'
	}
	if (`rc') exit `rc'
end


program Parse
	* Trim whitespace (caused by "///" line continuations; aesthetic only)
	mata: st_local("0", stritrim(st_local("0")))

	* Main syntax
	#d;
	syntax varlist(fv ts numeric) [if] [in] [aw pw fw/] , [

		/* Model */
		Absorb(string) NOAbsorb
		SUmmarize SUmmarize2(string asis) /* simulate implicit options */

		/* Standard Errors */
		VCE(string) CLuster(string)

		/* Diagnostic */
		Verbose(numlist min=1 max=1 >=-1 <=5 integer)
		TIMEit

		/* Speedup and memory Tricks */
		NOSAMPle /* do not save e(sample) */
		COMPACT /* use as little memory as possible but is slower */

		/* Extra display options (based on regress) */
		noHEader noTABle noFOOTnote
		
		/* Undocumented */
		KEEPSINgletons
		OLD /* use latest v3 */
		NOTES(string) /* NOTES(key=value ...), will be stored on e() */

		] [*] /* capture optimization options, display options, etc. */
		;
	#d cr

	* Unused
	* SAVEcache
	* USEcache
	* CLEARcache

	* Convert options to boolean
	if ("`verbose'" == "") loc verbose 0
	loc timeit = ("`timeit'"!="")
	loc drop_singletons = ("`keepsingletons'" == "")
	loc compact = ("`compact'" != "")

	if (`timeit') timer on 29

	* Sanity checks
	if (`verbose'>-1 & "`keepsingletons'"!="") {
		loc url "http://scorreia.com/reghdfe/nested_within_cluster.pdf"
		loc msg "WARNING: Singleton observations not dropped; statistical significance is biased"
		di as error `"`msg' {browse "`url'":(link)}"'
	}
	if ("`cluster'"!="") {
		_assert ("`vce'"==""), msg("cannot specify both cluster() and vce()")
		loc vce cluster `cluster'
		loc cluster // clear it to avoid bugs in subsequent lines
	}

	* Split varlist into <depvar> and <indepvars>
	ms_parse_varlist `varlist'
	if (`verbose' > 0) {
		di as text _n "## Parsing varlist: {res}`varlist'"
		return list
	}
	loc depvar `r(depvar)'
	loc indepvars `r(indepvars)'
	loc fe_format "`r(fe_format)'"
	loc basevars `r(basevars)'

	* Parse Weights
	if ("`weight'"!="") {
		unab exp : `exp', min(1) max(1) // simple weights only
	}

	* Parse VCE
	ms_parse_vce, vce(`vce') weighttype(`weight')
	if (`verbose' > 0) {
		di as text _n "## Parsing vce({res}`vce'{txt})"
		sreturn list
	}
	loc vcetype = "`s(vcetype)'"
	loc num_clusters = `s(num_clusters)'
	loc clustervars = "`s(clustervars)'"
	loc base_clustervars = "`s(base_clustervars)'"
	loc vceextra = "`s(vceextra)'"

	* Select sample (except for absvars)
	loc varlist `depvar' `indepvars' `base_clustervars'
	tempvar touse
	marksample touse, strok // based on varlist + cluster + if + in + weight

	* Parse noabsorb
	_assert  ("`absorb'`noabsorb'" != ""), msg("option {bf:absorb()} or {bf:noabsorb} required")
	if ("`noabsorb'" != "") {
		_assert ("`absorb'" == ""), msg("{bf:absorb()} and {bf:noabsorb} are mutually exclusive")
	}

	if (`timeit') timer off 29

	* Construct HDFE object
	// SYNTAX: fixed_effects(absvars | , touse, wtype, wtvar, dropsing, verbose)
	ms_add_comma, loc(absorb) cmd(`"`absorb'"') opt(`"`options'"')
	if (`timeit') timer on 20
	mata: HDFE = fixed_effects(`"`absorb'"', "`touse'", "`weight'", "`exp'", `drop_singletons', `verbose')
	if (`timeit') timer off 20
	mata: HDFE.cmdline = "reghdfe " + st_local("0")
	loc options `s(options)'

	mata: st_local("N", strofreal(HDFE.N))
	if (`N' == 0) error 2000

	* Fill out HDFE object
	* mata: HDFE.varlist = "`base_varlist'"
	mata: HDFE.depvar = "`depvar'"
	mata: HDFE.indepvars = "`indepvars'"
	mata: HDFE.vcetype = "`vcetype'"
	mata: HDFE.num_clusters = `num_clusters'
	mata: HDFE.clustervars = tokens("`clustervars'")
	mata: HDFE.base_clustervars = tokens("`base_clustervars'")
	mata: HDFE.vceextra = "`vceextra'"

	* Preserve memory
	mata: HDFE.compact = `compact'
	if (`compact') {
		loc panelvar "`_dta[_TSpanel]'"
		loc timevar "`_dta[_TStvar]'"

		cap conf var `panelvar', exact
		if (c(rc)) loc panelvar

		cap conf var `timevar', exact
		if (c(rc)) loc timevar

		mata: HDFE.panelvar = "`panelvar'"
		mata: HDFE.timevar = "`timevar'"
		c_local keepvars `basevars' `base_clustervars' `panelvar' `timevar' // `exp'
	}

	* Parse summarize
	if ("`summarize'" != "") {
		_assert ("`summarize2'" == ""), msg("summarize() syntax error")
		loc summarize2 mean min max  // default values
	}
	ParseSummarize `summarize2'
	mata: HDFE.summarize_stats = "`s(stats)'"
	mata: HDFE.summarize_quietly = `s(quietly)'


	* Parse misc options
	mata: HDFE.notes = `"`notes'"'
	mata: HDFE.store_sample = ("`nosample'"=="")
	mata: HDFE.timeit = `timeit'


	* Parse Coef Table Options (do this last!)
	_get_diopts diopts options, `options' // store in `diopts', and the rest back to `options'
	loc diopts `diopts' `header' `table' `footnote'
	_assert (`"`options'"'==""), msg(`"invalid options: `options'"')
	if ("`hascons'"!="") di in ye "(option ignored: `hascons')"
	if ("`tsscons'"!="") di in ye "(option ignored: `tsscons')"
	mata: HDFE.diopts = `"`diopts'"'
end


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
	mata: st_local("timeit", strofreal(HDFE.timeit))
	mata: st_local("compact", strofreal(HDFE.compact))
	mata: st_local("verbose", strofreal(HDFE.verbose))

	* Compute degrees-of-freedom
	if (`timeit') timer on 21
	mata: HDFE.estimate_dof()
	if (`timeit') timer off 21

	* Save updated e(sample) (singletons reduce sample);
	* required to parse factor variables to partial out
	if (`timeit') timer on 29
	tempvar touse
	mata: HDFE.save_touse("`touse'")
	if (`timeit') timer off 29

	* Expand varlists
	if (`timeit') timer on 22
	mata: st_local("depvar", HDFE.depvar)
	mata: st_local("indepvars", HDFE.indepvars)
	if (`verbose' > 0) di as text _n "## Parsing and expanding indepvars: {res}`indepvars'"
	ms_expand_varlist `indepvars' if `touse'
	if (`verbose' > 0) return list
	mata: HDFE.fullindepvars = "`r(fullvarlist)'"
	mata: HDFE.indepvars = "`r(varlist)'"
	mata: HDFE.not_basevar = strtoreal(tokens("`r(not_omitted)'"))
	mata: HDFE.varlist = "`depvar' `r(varlist)'"
	if (`timeit') timer off 22

	* Stats
	mata: st_local("stats", HDFE.summarize_stats)
	if ("`stats'" != "") Stats `touse'

	* Condition number
	mata: HDFE.estimate_cond()

	* Preserve
	if (`compact') {
		if (`verbose' > 0) di as text "## Preserving dataset"
		preserve
		novarabbrev keep `keepvars'
	}

	* Partial out; save TSS of depvar
	if (`timeit') timer on 23
	// SYNTAX: partial_out(Varlist/Matrix | , Save TSS if HDFE.tss is missing? [0], Standardize data? [1], First col is depvar? [1])
	// Note: standardize=2 will standardize, partial out, and return the data standardized!
	mata: hdfe_variables = HDFE.partial_out(HDFE.varlist, 1, 2, .)
	if (`timeit') timer off 23

	* Regress
	if (`timeit') timer on 24
	tempname b V N rank df_r
	mata: reghdfe_post_ols(HDFE, hdfe_variables, "`b'", "`V'", "`N'", "`rank'", "`df_r'")
	mata: hdfe_variables = .
	* Restore
	if (`compact') {
		if (`verbose' > 0) di as text "## Restoring dataset"
		restore
		mata: st_local("residuals", HDFE.residuals)
		if ("`residuals'" != "") mata: HDFE.save_variable(HDFE.residuals, HDFE.residuals_vector, "Residuals")
	}
	RegressOLS `touse' `b' `V' `N' `rank' `df_r'
	if (`timeit') timer off 24

	* (optional) Store FEs
	if (`timeit') timer on 29
	reghdfe5, store_alphas
	if (`timeit') timer off 29

	* View estimation tables
	mata: st_local("diopts", HDFE.diopts)
	Replay, `diopts'

	if (`timeit') {
		di as text _n "{bf: Timer results:}"
		timer list
		di as text "Legend: 20: Create HDFE object; 21: Estimate DoF; 22: expand varlists; 23: partial out; 24: regress; 29: rest"
		di
	}
end


program RegressOLS, eclass
	args touse b V N rank df_r

	mata: st_local("store_sample", strofreal(HDFE.store_sample))
	if (`store_sample') loc esample "esample(`touse')"
	
	mata: st_local("indepvars", HDFE.fullindepvars)
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

	ereturn scalar N       = `N'
	ereturn scalar rank    = `rank'
	ereturn scalar df_r    = `df_r'
	ereturn local  cmd     "reghdfe"
	mata: HDFE.post()

	* Post stats
	cap conf matrix reghdfe_statsmatrix
	if (!c(rc)) {
		ereturn matrix summarize = reghdfe_statsmatrix
		mata: st_local("summarize_quietly", strofreal(HDFE.summarize_quietly))
		ereturn scalar summarize_quietly = `summarize_quietly'
	}
end


program Replay, rclass
	syntax [, noHEader noTABle noFOOTnote *]

	if `"`e(cmd)'"' != "reghdfe"  {
	        error 301
	}

	_get_diopts options, `options'
	if ("`header'" == "") {
		reghdfe5_header // _coef_table_header
		di ""
	}
	if ("`table'" == "") {
		_coef_table, `options' // ereturn display, `options'
		return add // adds r(level), r(table), etc. to ereturn (before the footnote deletes them)
	}
	if ("`footnote'" == "") {
		reghdfe5_footnote
	}

	* Replay stats
	if (e(summarize_quietly)==0) {
		di as text _n "{sf:Regression Summary Statistics:}" _c
		matlist e(summarize)', border(top bottom) rowtitle(Variable) // twidth(18) 
	}
end


program Stats
	args touse
	* Optional weights
	mata: st_local("weight", sprintf("[%s=%s]", HDFE.weight_type, HDFE.weight_var))
	assert "`weight'" != ""
	if ("`weight'" == "[=]") loc weight
	loc weight : subinstr local weight "[pweight" "[aweight"

	mata: st_local("stats", HDFE.summarize_stats)
	mata: st_local("varlist", HDFE.varlist)
	mata: st_local("cvars", invtokens(HDFE.cvars))
	loc full_varlist `varlist' `cvars'

	* quick workaround b/c -tabstat- does not support factor variables
	fvrevar `full_varlist', list
	loc full_varlist `r(varlist)'

	qui tabstat `full_varlist' if `touse' `weight' , stat(`stats') col(stat) save
	matrix reghdfe_statsmatrix = r(StatTotal)
end

include "reghdfe5.mata", adopath

exit
