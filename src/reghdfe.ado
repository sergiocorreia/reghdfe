*! version 4.0.4 dev - 15feb2017

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
	if inlist("`options'", "check", "compile", "requirements", "setup", "update", "version") {
		if ("`options'"=="compile") loc args force
		if ("`options'"=="check") loc options compile
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

program Requirements
	di as text _n "{bf:[reghdfe]} - updating required packages"
	di as text "{dup 64:=}"
	* -moresyntax- https://github.com/sergiocorreia/moresyntax/
	cap ado uninstall moresyntax
	//net install moresyntax, from(https://github.com/sergiocorreia/moresyntax/raw/master/source/)
	net install moresyntax, from(c:\git\moresyntax\src)
	di as text "{dup 64:-}"

	* -ftools- https://github.com/sergiocorreia/ftools/
	cap ado uninstall ftools
	//net install ftools, from(https://github.com/sergiocorreia/ftools/raw/master/source/)
	net install ftools, from(c:\git\ftools\src)
	di as text "{dup 64:-}"
	ftools, compile // requires moresyntax
	di as text "{dup 64:-}"
end


program Setup
	* Update requirements/dependencies
	Requirements

	* Compile -reghdfe-
	loc source "c:\git\reghdfe\src"
	di as text _n  _n "{bf:[reghdfe]} - compiling Mata code for {res}reghdfe"
	di as text "{dup 64:=}"
	reghdfe, compile
	di as text "{dup 64:-}"
end


program Update
	* Update requirements/dependencies
	Requirements

	* Update -reghdfe-
	loc source "c:\git\reghdfe\src"
	di as text _n  _n "{bf:[reghdfe]} - updating self from {res}`source'"
	di as text "{dup 64:=}"
	qui ado uninstall reghdfe
	net install reghdfe, from("`source'") // TODO: Change this
	qui which reghdfe
	di as text "{dup 64:-}"
	reghdfe, compile
	di as text "{dup 64:-}"
	di as text "{bf:Note:} You need to run {stata program drop _all} now."
end


program Version
	which reghdfe

	di as text _n "Dependencies installed?"
	local dependencies moresyntax ftools ivreg2 avar tuples
	foreach dependency of local dependencies {
		cap findfile `dependency'.ado
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


program Estimate, eclass
	ereturn clear

* Parse and fill out HDFE object
	mata: HDFE = FixedEffects()
	reghdfe_parse `0' // store results on HDFE, HDFE.options, and HDFE.output
	loc OPT HDFE.options
	loc OUT HDFE.output

* Select sample
	marksample touse, strok // based on locals for varlist+cluster+absvars+if+in+weight (injected by -reghdfe_parse-)
	// -strok- because clustervars can be string

* Construct HDFE object
	// SYNTAX: fixed_effects(absvars | , touse, weighttype, weightvar, dropsing, verbose, object)
	mata: HDFE = fixed_effects(`OPT'.absorb, "`touse'", `OPT'.weight_type, `OPT'.weight_var, `OPT'.drop_singletons, HDFE.verbose, HDFE)
	mata: st_local("N", strofreal(HDFE.N))
	if (`N' == 0) error 2000

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
