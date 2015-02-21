// -------------------------------------------------------------
// Parsing and basic sanity checks for REGHDFE.ado
// -------------------------------------------------------------
// depvar: dependent variable
// indepvars: included exogenous regressors
// endogvars: included endogenous regressors
// instruments: excluded exogenous regressors

cap pr drop Parse
program define Parse

* Remove extra spacing from cmdline (just for aesthetics, run before syntax)
	cap syntax anything(name=indepvars) [if] [in] [fweight aweight pweight/] , SAVEcache(string) [*]
	local savingcache = (`=_rc'==0)

if (`savingcache') {

	* Disable these options
	local fast
	local nested

	syntax anything(name=indepvars) [if] [in] [fweight aweight pweight/] , ///
		Absorb(string) SAVEcache(string) ///
		[Verbose(integer 0) CHECK TOLerance(real 1e-7) MAXITerations(integer 1000) noACCELerate ///
		bad_loop_threshold(integer 1) stuck_threshold(real 5e-3) pause_length(integer 20) ///
		accel_freq(integer 3) accel_start(integer 6) /// Advanced optimization options
		CORES(integer 1) OVER(varname numeric)]

	cap conf file "`savecache'.dta"
	if (`=_rc'!=0) {
		cap conf new file "`savecache'.dta"
		Assert (`=_rc'==0), msg("reghdfe will not be able to save `savecache'.dta")
	}

}
else {
	mata: st_local("cmdline", stritrim(`"reghdfe `0'"') )
	ereturn clear // Clear previous results and drops e(sample)
	syntax anything(id="varlist" name=0 equalok) [if] [in] ///
		[fweight aweight pweight/] , ///
		Absorb(string) ///
		[GROUP(name) VCE(string) Verbose(integer 0) CHECK NESTED FAST] ///
		[avge(string) EXCLUDESELF] ///
		[TOLerance(real 1e-7) MAXITerations(integer 1000) noACCELerate] /// See reghdfe_absorb.Annihilate
		[noTRACK] /// Not used here but in -Track-
		[IVsuite(string) ESTimator(string) SAVEFIRST FIRST SHOWRAW dofminus(string)] ///
		[dofmethod(string)] ///
		[SMALL Hascons TSSCONS] /// ignored options
		[gmm2s liml kiefer cue] ///
		[SUBOPTions(string)] /// Options to be passed to the estimation command (e.g . to regress)
		[bad_loop_threshold(integer 1) stuck_threshold(real 5e-3) pause_length(integer 20) ///
		accel_freq(integer 3) accel_start(integer 6)] /// Advanced optimization options
		[CORES(integer 1)] [USEcache(string)] [OVER(varname numeric)] ///
		[noCONstant] /// Disable adding back the intercept (mandatory with -ivreg2-)
		[*] // For display options
}

* Weight
* We'll have -weight- (fweight|aweight|pweight), -weightvar-, -exp-, and -weightexp-
	if ("`weight'"!="") {
		local weightvar `exp'
		conf var `weightvar' // just allow simple weights
		local weightexp [`weight'=`weightvar']
		local backupweight `weight'
	}

* Cache options
	if ("`usecache'"!="") {
		conf file "`usecache'.dta"
		conf var __uid__
		Assert ("`avge'"==""), msg("option -avge- not allowed with -usecache-")
		Assert ("`avge'"==""), msg("option -nested- not allowed with -usecache-")
	}

* Save locals that will be overwritten by later calls to -syntax-
	local ifopt `if'
	local inopt `in'

* Coef Table Options
if (!`savingcache') {
	_get_diopts diopts options, `options'
	Assert `"`options'"'=="", msg(`"invalid options: `options'"')
	if ("`hascons'`tsscons'"!="") di in ye "(option `hascons'`tsscons' ignored)"
}

* Over
	if ("`over'"!="") {
		unab over : `over', max(1)
		Assert ("`usecache'"!="" | "`savecache'"!=""), msg("-over- needs to be used together with either -usecache- or -savecache-")
	}

* Verbose
	assert inlist(`verbose', 0, 1, 2, 3, 4) // 3 and 4 are developer options
	mata: VERBOSE = `verbose' // Ugly hack to avoid using a -global-

* Show raw output of called subcommand (e.g. ivreg2)
	local showraw = ("`showraw'"!="")

* Model settings
if (!`savingcache') {

	// local model = cond(strpos(`"`0'"', " ("), "iv", "ols") // Fails with long strs in stata 12<
	local model ols
	foreach _ of local 0 {
		if (substr(`"`_'"', 1, 1)=="(") {
			local model iv
			continue, break
		}
	}
	

	if ("`model'"=="iv") {

		// get part before parentheses
		local wrongparens 1
		while (`wrongparens') {
			gettoken tmp 0 : 0 ,p("(")
			local left `left'`tmp'
			* Avoid matching the parens of e.g. L(-1/2) and L.(var1 var2)
			* Using Mata to avoid regexm() and trim() space limitations
			mata: st_local("tmp1", subinstr("`0'", " ", "") ) // wrong parens if ( and then a number
			mata: st_local("tmp2", substr(strtrim("`left'"), -1) ) // wrong parens if dot
			local wrongparens = regexm("`tmp1'", "^\([0-9-]") | ("`tmp2'"==".")
			if (`wrongparens') {
				gettoken tmp 0 : 0 ,p(")")
				local left `left'`tmp'
			}
		}

		// get part in parentheses
		gettoken right 0 : 0 ,bind match(parens)
		Assert trim(`"`0'"')=="" , msg("error: remaining argument: `0'")

		// now parse part in parentheses
		gettoken endogvars instruments : right ,p("=")
		gettoken equalsign instruments : instruments ,p("=")

		Assert "`endogvars'"!="", msg("iv: endogvars required")
		local 0 `endogvars'
		syntax varlist(fv ts numeric)

		Assert "`instruments'"!="", msg("iv: instruments required")
		local 0 `instruments'
		syntax varlist(fv ts numeric)
		
		local 0 `left' // So OLS part can handle it
		Assert "`endogvars'`instruments'"!=""
		
		if ("`ivsuite'"=="") local ivsuite ivreg2
		Assert inlist("`ivsuite'","ivreg2","ivregress") , msg("error: wrong IV routine (`ivsuite'), valid options are -ivreg2- and -ivregress-")
		cap findfile `ivsuite'.ado
		Assert !_rc , msg("error: -`ivsuite'- not installed, please run ssc install `ivsuite' or change the option -ivsuite-")
	}

* OLS varlist
	syntax varlist(fv ts numeric)
	gettoken depvar indepvars : 0
	_fv_check_depvar `depvar'

* Extract format of depvar so we can format FEs like this
	fvrevar `depvar', list
	local fe_format : format `r(varlist)' // The format of the FEs and AvgEs that will be saved

* Variables shouldn't be repeated
* This is not perfect (e.g. doesn't deal with "x1-x10") but still helpful
	local allvars `depvar' `indepvars' `endogvars' `instruments'
	local dupvars : list dups allvars
	Assert "`dupvars'"=="", msg("error: there are repeated variables: <`dupvars'>")

	Debug, msg(_n " {title:REGHDFE} Verbose level = `verbose'")
	*Debug, msg("{hline 64}")

* Add back constants (place this *after* we define `model')
	local addconstant = ("`constant'"!="noconstant") & !("`model'"=="iv" & "`ivsuite'"=="ivreg2")

* VCE
	gettoken vcetype vcerest : vce, parse(" ,")
	if ("`vcetype'"=="") {
		local vcetype unadjusted
		if ("`ivsuite'"=="ivregress" & "`estimator'"=="gmm") local vcetype robust // Default for ivregress gmm
		if ("`backupweight'"=="pweight") local vcetype robust // pweight requires robust
	}

	* Abbreviations:
	if (substr("`vcetype'",1,2)=="un") local vcetype unadjusted
	if (substr("`vcetype'",1,1)=="r") local vcetype robust
	if (substr("`vcetype'",1,2)=="cl") local vcetype cluster
	if (substr("`vcetype'",1,4)=="boot") local vcetype bootstrap
	if (substr("`vcetype'",1,4)=="jack") local vcetype jackknife
	if ("`vcetype'"=="conventional") local vcetype unadjusted // Conventional is the name given in e.g. xtreg
	assert inlist("`vcetype'", "unadjusted", "robust", "cluster", "bootstrap", "jackknife")

	* Cluster vars
	local num_clusters 0
	if ("`vcetype'"=="cluster") {
		local num_clusters = `: word count `vcerest''
		assert inlist(`num_clusters', 1, 2)
		gettoken clustervar1 clustervar2 : vcerest
		cap unab clustervar1 : `clustervar1'
		cap unab clustervar2 : `clustervar2'
		local vcerest
	}

	//local vceoption = trim("`vcetype' `clustervar1' `clustervar2' `vcerest'")
	local temp_clustervar1 = cond("`clustervar1'"=="","","<CLUSTERVAR1>")
	local vceoption = trim("`vcetype' `temp_clustervar1'`vcerest'")
	local vceoption vce(`vceoption')
	
	if ("`model'"=="ols") {
		local ok = inlist("`vcetype'", "unadjusted", "robust", "cluster") & (`num_clusters'<=1)
	}
	else if ("`ivsuite'"=="ivregress") {
		local ok = inlist("`vcetype'", "unadjusted", "robust", "cluster", "bootstrap", "jackknife", "hac") & (`num_clusters'<=1)
	}
	else if ("`ivsuite'"=="ivreg2") {
		local ok = inlist("`vcetype'", "unadjusted", "robust", "cluster", "bw") & (`num_clusters'<=2)
	}
	Assert `ok', msg("invalid vce type or number of clusters")
}

* Optimization
	if (`maxiterations'==0) local maxiterations 1e7
	Assert (`maxiterations'>0)
	local accelerate = cond("`accelerate'"!="", 0, 1) // 1=Yes
	local check = cond("`check'"!="", 1, 0) // 1=Yes
	local fast = cond("`fast'"!="", 1, 0) // 1=Yes
	local tolerance = strofreal(`tolerance', "%9.1e") // Purely esthetic
	Assert `cores'<=32 & `cores'>0 , msg("At most 32 cores supported")
	if (`cores'>1) {
		cap findfile parallel.ado
		Assert !_rc , msg("error: -parallel- not installed, please run {cmd:ssc install parallel}")
	}
	local opt_list tolerance maxiterations check accelerate ///
		bad_loop_threshold stuck_threshold pause_length accel_freq accel_start
	foreach opt of local opt_list {
		if ("``opt''"!="") local maximize_options `maximize_options' `opt'(``opt'')
	}

* IV options
if (!`savingcache') {
	if ("`model'"=="iv" & "`ivsuite'"=="ivregress") {
		if ("`estimator'"=="") local estimator 2sls
		Assert inlist("`estimator'","2sls","gmm"), msg("liml estimator not allowed")
		if ("`estimator'"!="2sls") di as error "WARNING! GMM not fully implemented, robust option gives wrong output"
	}
	else if ("`model'"=="iv" & "`ivsuite'"=="ivreg2") {
		if ("`estimator'"=="") local estimator 2sls
	Assert inlist("`estimator'","2sls","gmm2s") , msg("Estimator is `estimator', instead of empty (2sls/ols) or gmm2s")
	* di in ye "WARNING: IV estimates not fully tested; methods like GMM2S will not work correctly"

		Assert inlist("`dofminus'","","large","small"), msg("option -dofminus- is either -small- or -large-")
		local dofminus = cond("`dofminus'"=="large", "dofminus", "sdofminus") // Default uses sdofminus
	}
	else {
		local estimator
	}
	if ("`small'"!="") di in ye "(note: reghdfe will always use the option -small-, no need to specify it)"

	Assert ("`gmm2s'`liml'`kiefer'`cue'"==""), msg("Please use estimator() option instead of gmm2s/liml/etc")
	
	if ("`model'"=="iv") {
		local savefirst = ("`savefirst'"!="")
		local first = ("`first'"!="")
		if (`savefirst') Assert `first', msg("Option -savefirst- requires -first-")
	}

* DoF Estimation
* For more than 2 FEs, no exact soln
* Can estimate M3 and so on using
* i) bounds using FE1 and FE2, ii) assuming M3=1 (fast), iii) bootstrap
	if ("`dofmethod'"=="") local dofmethod bounds
	assert inlist("`dofmethod'", "naive", "simple", "bounds", "bootstrap")

* Mobility groups
	if ("`group'"!="") conf new var `group'
}

* tsset variables, if any
	cap conf var `_dta[_TStvar]'
	if (!_rc) local timevar `_dta[_TStvar]'
	cap conf var `_dta[_TSpanel]'
	if (!_rc) local panelvar `_dta[_TSpanel]'

* Varnames underlying tsvars and fvvars (e.g. i.foo L(1/3).bar -> foo bar)
	foreach vars in depvar indepvars endogvars instruments {
		if ("``vars''"!="") {
			fvrevar ``vars'' , list
			local basevars `basevars' `r(varlist)'
		}
	}

if (!`savingcache') {
* Nested
	local nested = cond("`nested'"!="", 1, 0) // 1=Yes
	if (`nested' & !("`model'"=="ols" & "`vcetype'"=="unadjusted") ) {
		Debug, level(0) msg("(option nested ignored, only works with OLS and conventional/unadjusted VCE)") color("error")
	}

* How can we do the same regression from a standard stata command?
* (useful for benchmarking and testing correctness of results)
	local subcmd = cond("`model'"=="ols" ,"regress", "`ivsuite'")

* _fv_check_depvar overwrites the local -weight-
	local weight `backupweight'
	Assert inlist( ("`weight'"!="") + ("`weightvar'"!="") + ("`weightexp'"!="") , 0 , 3 ) , msg("not all 3 weight locals are set")

* Return values
	local names cmdline diopts model ///
		ivsuite estimator showraw dofminus ///
		depvar indepvars endogvars instruments savefirst first ///
		vceoption vcetype num_clusters clustervar1 clustervar2 vcerest ///
		if in group dofmethod check fast nested fe_format ///
		tolerance maxiterations accelerate maximize_options ///
		subcmd suboptions ///
		absorb avge excludeself ///
		timevar panelvar basevars ///
		addconstant ///
		weight weightvar exp weightexp /// type of weight (fw,aw,pw), weight var., and full expr. ([fw=n])
		cores savingcache usecache over
}

if (`savingcache') {
	local names maximize_options cores if in timevar panelvar indepvars basevars ///
		absorb savecache savingcache fast nested check over ///
		weight weightvar exp weightexp /// type of weight (fw,aw), weight var., and full expr. ([fw=n])
		tolerance maxiterations // Here just used for -verbose- and cache handshake purposes
}

	local if `ifopt'
	local in `inopt'

	Debug, level(3) newline
	Debug, level(3) msg("Parsed options:")
	foreach name of local names {
		if (`"``name''"'!="") Debug, level(3) msg("  `name' = " as result `"``name''"')
		c_local `name' `"``name''"' // Inject values into caller (reghdfe.ado)
	}
	// Debug, level(3) newline
end
