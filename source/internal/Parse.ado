// -------------------------------------------------------------
// Parsing and basic sanity checks for REGHDFE.ado
// -------------------------------------------------------------

cap pr drop Parse
program define Parse

* Remove extra spacing from cmdline (just for aesthetics)
	mata: st_local("cmdline", stritrim(`"reghdfe `0'"') )

* Parse the broad syntax (also see map_init(), ParseAbsvars.ado, ParseVCE.ado, etc.)
	syntax anything(id="varlist" name=0 equalok) [if] [in] [aw pw fw/] , ///
		/// Model ///
		Absorb(string) [ ///
		RESiduals(name) ///
		SUBOPTions(string) /// Options to be passed to the estimation command (e.g . to regress)
		/// Standard Errors ///
		VCE(string) CLuster(string) /// cluster() is an undocumented alternative to vce(cluster ...)
		/// IV/2SLS/GMM ///
		ESTimator(string) /// 2SLS GMM2s CUE LIML
		STAGEs(string) /// besides iv (always on), first reduced ols acid (and all)
		FFirst /// Save first-stage stats (only with ivreg2)
		IVsuite(string) /// ivreg2 or ivregress
		/// Diagnostic ///
		Verbose(string) ///
		TIMEit ///
		/// Optimization /// Defaults are handled within Mata		
		TOLerance(string) ///
		MAXITerations(string) ///
		POOLsize(string) /// Process variables in batches of #
		ACCELeration(string) ///
		TRAnsform(string) ///
		/// Speedup Tricks ///
		CACHE(string) ///
		FAST ///
		/// Degrees-of-freedom Adjustments ///
		DOFadjustments(string) ///
		GROUPVar(name) /// Variable that will contain the first connected group between FEs
		/// Undocumented ///
		KEEPSINgletons /// (UNDOCUMENTED) Will keep singletons
		NOTES(string) /// NOTES(key=value ...), will be stored on e()
		] [*] // Captures i) display options, ii) SUmmarize|SUmmarize(...)

	local allkeys cmdline if in timeit

* Do this early
	local timeit = "`timeit'"!=""
	local fast = "`fast'"!=""
	local ffirst = "`ffirst'"!=""
	
	if ("`cluster'"!="") {
		Assert ("`vce'"==""), msg("cannot specify both cluster() and vce()")
		local vce cluster `cluster'
		local cluster // Set it to empty to avoid bugs in subsequent lines
	}

* Also early
	ParseCache, cache(`cache') ifin(`if'`in') absorb(`absorb') vce(`vce')
	local keys savecache keepvars usecache
	foreach key of local keys {
		local `key' "`s(`key')'"
	}
	local allkeys `allkeys' `keys'

* Parse varlist: depvar indepvars (endogvars = iv_vars)
	ParseIV `0', estimator(`estimator') ivsuite(`ivsuite')
	local keys subcmd model ivsuite estimator depvar indepvars endogvars instruments fe_format basevars
	foreach key of local keys {
		local `key' "`s(`key')'"
	}
	local allkeys `allkeys' `keys'

* Weights
	if ("`weight'"!="") {
		local weightvar `exp'
		local weighttype `weight'
		local weightexp [`weight'=`weightvar']
		confirm var `weightvar', exact // just allow simple weights

		* Check that weights are correct (e.g. with fweight they need to be integers)
		local num_type = cond("`weight'"=="fweight", "integers", "reals")
		local basenote "{txt}weight {res}`weightvar'{txt} can only contain strictly positive `num_type', but"
		qui cou if `weightvar'<0
		Assert (`r(N)'==0), msg("`basenote' `r(N)' negative values were found!")  rc(402)
		qui cou if `weightvar'==0
		if (`r(N)'>0) di as text "`basenote' `r(N)' zero values were found (will be dropped)"
		qui cou if `weightvar'>=.
		if (`r(N)'>0) di as text "`basenote' `r(N)' missing values were found (will be dropped)"
		if ("`weight'"=="fweight") {
			qui cou if mod(`weightvar',1) & `weightvar'<.
			Assert (`r(N)'==0), msg("`basenote' `r(N)' non-integer values were found!" "{err} Stopping execution") rc(401)
		}
	}
	local allkeys `allkeys' weightvar weighttype weightexp

* Parse Absvars and optimization options
if (!`usecache') {
	ParseAbsvars `absorb' // Stores results in r()
		if (inlist("`verbose'", "4", "5")) return list
		local absorb_keepvars `r(all_ivars)' `r(all_cvars)'
		local N_hdfe `r(G)'
		local has_intercept = `r(has_intercept)'
		assert inlist(`has_intercept', 0, 1)

	mata: HDFE_S = map_init() // Reads results from r()
		local will_save_fe = `r(will_save_fe)' // Returned from map_init()
		local original_absvars "`r(original_absvars)'"
		local extended_absvars "`r(extended_absvars)'"
		local equation_d "`r(equation_d)'"
}
else {
	local will_save_fe 0
	local original_absvars : char _dta[original_absvars]
	local extended_absvars : char _dta[extended_absvars]
	local equation_d
	local N_hdfe : char _dta[N_hdfe]
	local has_intercept : char _dta[has_intercept]
}
	local allkeys `allkeys' absorb_keepvars N_hdfe will_save_fe original_absvars extended_absvars equation_d has_intercept

	* Tell Mata what weightvar we have
	if ("`weightvar'"!="" & !`usecache') mata: map_init_weights(HDFE_S, "`weightvar'", "`weighttype'")

	* Time/panel variables (need to give them to Mata)
	local panelvar `_dta[_TSpanel]'
	local timevar `_dta[_TStvar]'
	if ("`panelvar'"!="") {
		cap conf var `panelvar'
		if (c(rc)==111) local panelvar // if the var doesn't exist, set it empty
	}
	if ("`timevar'"!="") {
		cap conf var `timevar'
		if (c(rc)==111) local timevar // if the var doesn't exist, set it empty
	}

	* Parse optimization options (pass them to map_init_*)
	* String options
	local optlist transform acceleration panelvar timevar
	foreach opt of local optlist {
		if ("``opt''"!="" & !`usecache') mata: map_init_`opt'(HDFE_S, "``opt''")
	}
	local allkeys `allkeys' `optlist'

	* This allows changing the groupvar name with -usecache-
	if ("`groupvar'"!="") mata: map_init_groupvar(HDFE_S, "`groupvar'")

	* Numeric options
	local keepsingletons = ("`keepsingletons'"!="")
	local optlist poolsize verbose tolerance maxiterations keepsingletons timeit
	foreach opt of local optlist {
		if ( "``opt''"!="" & (!`usecache' | "`opt'"=="verbose") ) mata: map_init_`opt'(HDFE_S, ``opt'')
	}
	local allkeys `allkeys' `optlist'

	* Return back default value of -verbose-
	mata: verbose2local(HDFE_S, "verbose")
	local allkeys `allkeys' verbose

* Stages (before vce)
	ParseStages, stages(`stages') model(`model')
	local stages "`s(stages)'"
	local stage_suboptions "`s(stage_suboptions)'"
	local savestages = `s(savestages)'
	local allkeys `allkeys' stages stage_suboptions savestages

* Parse VCE options (after stages)
	local keys vceoption vcetype vcesuite vceextra num_clusters clustervars bw kernel dkraay kiefer twicerobust
	if (!`usecache') {
		mata: st_local("hascomma", strofreal(strpos("`vce'", ","))) // is there a commma already in `vce'?
		local vcetmp `vce'
		if (!`hascomma') local vcetmp `vce' ,
		ParseVCE `vcetmp' weighttype(`weighttype') ivsuite(`ivsuite') model(`model')
		foreach key of local keys {
			local `key' "`s(`key')'"
		}
	}
	else {
		foreach key of local keys {
			local `key' : char _dta[`key']
		}
	}

	local allkeys `allkeys' `keys'

* Parse FFIRST (save first stage statistics)
	local allkeys `allkeys' ffirst
	if (`ffirst') Assert "`model'"!="ols", msg("ols does not support {cmd}ffirst")
	if (`ffirst') Assert "`ivsuite'"=="ivreg2", msg("option {cmd}ffirst{err} requires ivreg2")
	
* Update Mata
	if ("`clustervars'"!="" & !`usecache') mata: map_init_clustervars(HDFE_S, "`clustervars'")
	if ("`vceextra'"!="" & !`usecache') mata: map_init_vce_is_hac(HDFE_S, 1)

* DoF Adjustments
	if ("`dofadjustments'"=="") local dofadjustments all
	ParseDOF , `dofadjustments'
	local dofadjustments "`s(dofadjustments)'"
	* Mobility groups
	if ("`groupvar'"!="") conf new var `groupvar'
	local allkeys `allkeys' dofadjustments groupvar

* Parse residuals
	if ("`residuals'"!="") {
		Assert !`will_save_fe', msg("option residuals() is mutually exclusive with saving fixed effects")
		Assert !`savecache', msg("option residuals() is mutually exclusive with -savecache-")
		conf new var `residuals'
		local allkeys `allkeys' residuals
	}

* Parse summarize option: [summarize | summarize( stats... [,QUIetly])]
	* Note: ParseImplicit deals with "implicit" options and fills their default values
	local default_stats mean min max
	ParseImplicit, opt(SUmmarize) default(`default_stats') input(`options') syntax([namelist(name=stats)] , [QUIetly]) inject(stats quietly)
	local summarize_quietly = ("`quietly'"!="")
	if ("`stats'"=="" & "`quietly'"!="") local stats `default_stats'
	local allkeys `allkeys' stats summarize_quietly

* Parse speedups
	if (`fast' & ("`groupvar'"!="" | `will_save_fe'==1 | "`residuals'"!="")) {
		di as error "(warning: option -fast- disabled; not allowed when saving variables: saving fixed effects, mobility groups, residuals)"
		local fast 0
	}
	local allkeys `allkeys' fast level

* Nested
	local nested = cond("`nested'"!="", 1, 0) // 1=Yes
	if (`nested' & !("`model'"=="ols" & "`vcetype'"=="unadjusted") ) {
		di as error "-nested- not implemented currently"
		Debug, level(0) msg("(option nested ignored, only works with OLS and conventional/unadjusted VCE)") color("error")
	}
	local allkeys `allkeys' nested

* Sanity checks on speedups
* With -savecache-, this adds chars (modifies the dta!) so put it close to the end
	if (`savecache') {
		* Savecache "requires" a previous preserve, so we can directly modify the dataset
		Assert "`endogvars'`instruments'"=="", msg("cache(save) option requires a normal varlist, not an iv varlist")
		char _dta[reghdfe_cache] 1
		local chars absorb N_hdfe has_intercept original_absvars extended_absvars vce vceoption vcetype vcesuite vceextra num_clusters clustervars bw kernel dkraay kiefer twicerobust
		foreach char of local  chars {
			char _dta[`char'] ``char''	
		}
	}

* Parse Coef Table Options (do this last!)
	_get_diopts diopts options, `options' // store in `diopts', and the rest back to `options'
	Assert `"`options'"'=="", msg(`"invalid options: `options'"')
	if ("`hascons'`tsscons'"!="") di in ye "(option `hascons'`tsscons' ignored)"
	local allkeys `allkeys' diopts

* Other keys:
	local allkeys `allkeys' suboptions notes
	// Missing keys: check

* Return values
	Debug, level(3) newline
	Debug, level(3) msg("{title:Parsed options:}")
	foreach key of local allkeys {
		if (`"``key''"'!="") Debug, level(3) msg("  `key' = " as result `"``key''"')
		c_local `key' `"``key''"' // Inject values into caller (reghdfe.ado)
	}

end
