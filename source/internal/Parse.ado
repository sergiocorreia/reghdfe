// -------------------------------------------------------------
// Parsing and basic sanity checks for REGHDFE.ado
// -------------------------------------------------------------

cap pr drop Parse
program define Parse

* Remove extra spacing from cmdline (just for aesthetics)
	mata: st_local("cmdline", stritrim(`"reghdfe `0'"') )

* Parse the broad syntax (also see map_init(), ParseAbsvars.ado, ParseVCE.ado, etc.)
	syntax anything(id="varlist" name=0 equalok) [if] [in] [fw aw pw/] , ///
	/// Main Options ///
		Absorb(string) ///
		[ ///
		VCE(string) ///
		Verbose(string) ///
	/// Seldom Used ///
		DOFadjustments(string) ///
		GROUPVAR(name) /// Variable that will contain the first connected group between FEs
	/// Optimization /// Defaults are handled within Mata		
		FAST /// Fast precludes i) saving FE, ii) running predict, iii) saving groupvar
		GROUPsize(string) /// Process variables in batches of #
		TRAnsform(string) ///
		ACCELeration(string) ///
		TOLerance(string) ///
		MAXITerations(string) ///
		KEEPSINgletons /// (UNDOCUMENTED) Will keep singletons
		CHECK /// TODO: Implement
		TIMEit ///
	/// Regression ///
		ESTimator(string) /// GMM2s CUE LIML
		IVsuite(string) ///
		SAVEFIRST ///
		FIRST ///
		SHOWRAW ///
		VCEUNADJUSTED /// (UNDOCUMENTED) Option when running gmm2s with ivregress; will match results of ivreg2
		SMALL Hascons TSSCONS /// ignored options
		SUBOPTions(string) /// Options to be passed to the estimation command (e.g . to regress)
	/// Multiple regressions in one go ///
		STAGEs(string) ///
		SAVEcache ///
		KEEPvars(varlist) ///
		USEcache ///
		BY(varname numeric) /// Requires savecache or usecache
		LEVEL(string) /// level of by (should be an integer actually), requires usecache
		NESTED /// TODO: Implement
	/// Miscellanea ///
		NOTES(string) /// NOTES(key=value ..)
		RESiduals(name) ///
		] [*] // For display options ; and SUmmarize(stats)

	local allkeys cmdline if in timeit

* Need to do this early
	local timeit = "`timeit'"!=""
	local fast = "`fast'"!=""
	local savecache = "`savecache'"!=""
	local usecache = "`usecache'"!=""
	if (!`usecache') {
		local is_cache : char _dta[reghdfe_cache]
		Assert "`is_cache'"!="1", msg("reghdfe error: data transformed with -savecache- requires -usecache-")
	}

* Sanity checks on usecache
* (this was at the end, but if there was no prev savecache, HDFE_S didn't exist and those lines were never reached)
	if (`usecache') {
		local is_cache : char _dta[reghdfe_cache]
		local by_cache : char _dta[by]
		local cache_obs : char _dta[cache_obs]
		local cache_absorb : char _dta[absorb]
		if ("`by'"!="") Assert "`level'"!="", msg("a previous -savecache by()- requires -usecache by() level()-")
		Assert "`by'"=="`by_cache'", msg("by() needs to be the same as in savecache")
		Assert "`is_cache'"=="1" , msg("usecache requires a previous savecache operation")
		Assert `cache_obs'==`c(N)', msg("dataset cannot change after savecache")
		Assert "`cache_absorb'"=="`absorb'", msg("cache dataset has different absorb()")
		Assert "`if'`in'"=="", msg("cannot use if/in with usecache; data has already been transformed")
	}

* Parse varlist: depvar indepvars (endogvars = iv_vars)
	ParseIV `0', estimator(`estimator') ivsuite(`ivsuite') `savefirst' `first' `showraw' `vceunadjusted' `small'
	local keys subcmd model ivsuite estimator depvar indepvars endogvars instruments fe_format ///
		savefirst first showraw vceunadjusted basevars
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
		local basenote "weight {res}`weightvar'{txt} can only contain strictly positive `num_type', but"
		qui cou if `weightvar'<0
		Assert (`r(N)'==0), msg("`basenote' `r(N)' negative values were found!")  rc(402)
		qui cou if `weightvar'==0
		if (`r(N)'>0) di as text "`basenote' `r(N)' zero values were found (will be dropped)"
		qui cou if `weightvar'>=.
		if (`r(N)'>0) di as text "`basenote' `r(N)' missing values were found (will be dropped)"
		if ("`weight'"=="fweight") {
			qui cou if mod(`weightvar',1) & `weightvar'<.
			Assert (`r(N)'==0), msg("`basenote' `r(N)' non-integer values were found!") rc(401)
		}
	}
	local allkeys `allkeys' weightvar weighttype weightexp

* Parse Absvars and optimization options
if (!`usecache') {
	if ("`by'"!="") {
		mata: st_local("has_comma", strofreal(strpos("`absorb'", ",")>0) )
		local bycomma = cond(`has_comma', "by(`by')", ", by(`by')")
	}
	ParseAbsvars `absorb' `bycomma' // Stores results in r()
		local absorb_keepvars `r(all_ivars)' `r(all_cvars)'
		local N_hdfe `r(G)'

	mata: HDFE_S = map_init("`by'") // Reads results from r()
		local will_save_fe = `r(will_save_fe)' // Returned from map_init()
		local original_absvars = "`r(original_absvars)'"
		local extended_absvars = "`r(extended_absvars)'"
		local equation_d = "`r(equation_d)'"
}
else {
	local will_save_fe 0
	local original_absvars : char _dta[original_absvars]
	local extended_absvars : char _dta[extended_absvars]
	local equation_d
	local N_hdfe : char _dta[N_hdfe]
}
	local allkeys `allkeys' absorb_keepvars N_hdfe will_save_fe original_absvars extended_absvars equation_d

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
	local optlist groupsize verbose tolerance maxiterations keepsingletons timeit
	foreach opt of local optlist {
		if ( "``opt''"!="" & (!`usecache' | "`opt'"=="verbose") ) mata: map_init_`opt'(HDFE_S, ``opt'')
	}
	local allkeys `allkeys' `optlist'

	* Return back default value of -verbose-
	mata: verbose2local(HDFE_S, "verbose")
	local allkeys `allkeys' verbose

* Stages (before vce)
	assert "`model'"!="" // just to be sure this goes after `model' is set
	if ("`stages'"=="all") local stages iv ols first acid reduced
	local iv_stage iv
	local stages : list stages - iv_stage
	local valid_stages ols first acid reduced
	local wrong_stages : list stages - valid_stages
	Assert "`wrong_stages'"=="", msg("Error, invalid stages(): `wrong_stages'")
	if ("`stages'"!="") {
		Assert "`model'"=="iv", msg("Error, stages() only valid with an IV regression")
		local stages `stages' `iv_stage' // Put -iv- *last* (so it does the -restore-; note that we don't need it first to trim MVs b/c that's done earlier)
	}
	else {
		local stages none // So we can loop over stages
	}
	local allkeys `allkeys' stages

* Parse VCE options (after stages)
	mata: st_local("hascomma", strofreal(strpos("`vce'", ","))) // is there a commma already in `vce'?
	local keys vceoption vcetype vcesuite vceextra num_clusters clustervars bw kernel dkraay kiefer twicerobust
	if (!`usecache') {
		local vcetmp `vce'
		if (!`hascomma') local vcetmp `vce' ,
		ParseVCE `vcetmp' weighttype(`weighttype') stages(`stages') ivsuite(`ivsuite') model(`model')
		foreach key of local keys {
			local `key' "`s(`key')'"
		}
	}
	else {
		local cache_vce : char _dta[vce]
		assert_msg "`cache_vce'"=="`vce'", msg("vce() must be the same in savecache and usecache (because of the cluster variables)")
		foreach key of local keys {
			local `key' : char _dta[`key']
		}
	}

	local allkeys `allkeys' `keys'
	
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
 	if ("`by'"!="") {
		unab by : `by', max(1)
	}
	if ("`keepvars'"!="" & !`savecache') di as error "(warning: keepvars() has no effect without savecache)"
	local allkeys `allkeys' fast savecache keepvars usecache by level

* Nested
	local nested = cond("`nested'"!="", 1, 0) // 1=Yes
	if (`nested' & !("`model'"=="ols" & "`vcetype'"=="unadjusted") ) {
		di as error "-nested- not implemented currently"
		Debug, level(0) msg("(option nested ignored, only works with OLS and conventional/unadjusted VCE)") color("error")
	}
	local allkeys `allkeys' nested

* Sanity checks on speedups
* With -savecache-, this adds chars (modifies the dta!) so put it close to the end
	Assert `usecache' + `savecache' < 2, msg("savecache and usecache are mutually exclusive")
	if ("`by'`level'"!="") di as error "(warning: by() and level() are currently incomplete)"
	if ("`by'"!="") Assert `usecache' + `savecache' == 1 , msg("by() requires savecache or usecache")
	if ("`level'"!="") Assert `usecache'==1 & "`by'"!="", msg("level() requires by() and usecache")
	if (`savecache') {
		* Savecache "requires" a previous preserve, so we can directly modify the dataset
		Assert "`endogvars'`instruments'"=="", msg("savecache option requires a normal varlist, not an iv varlist")
		char _dta[reghdfe_cache] 1
		local chars absorb N_hdfe original_absvars extended_absvars by vce vceoption vcetype vcesuite vceextra num_clusters clustervars bw kernel dkraay kiefer twicerobust
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

