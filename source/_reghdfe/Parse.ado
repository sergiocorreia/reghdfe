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
		[VCE(string)] ///
		[DOFadjustments(string) GROUP(name)] ///
		[avge(string) EXCLUDESELF] ///
		[Verbose(integer 0) CHECK NESTED FAST] ///
		[TOLerance(real 1e-7) MAXITerations(integer 1000) noACCELerate] /// See reghdfe_absorb.Annihilate
		[noTRACK] /// Not used here but in -Track-
		[IVsuite(string) SAVEFIRST FIRST SHOWRAW] /// ESTimator(string)
		[SMALL Hascons TSSCONS] /// ignored options
		[gmm2s liml kiefer cue] ///
		[SUBOPTions(string)] /// Options to be passed to the estimation command (e.g . to regress)
		[bad_loop_threshold(integer 1) stuck_threshold(real 5e-3) pause_length(integer 20) accel_freq(integer 3) accel_start(integer 6)] /// Advanced optimization options
		[CORES(integer 1)] [USEcache(string)] [OVER(varname numeric)] ///
		[POSTestimation(string) NOTES(string)] /// (Quipu) postestimation([SUmmarize QUIetly]) NOTES(key=value ..)
		[STAGEs(string)] ///
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

* tsset variables, if any
	cap conf var `_dta[_TStvar]'
	if (!_rc) local timevar `_dta[_TStvar]'
	cap conf var `_dta[_TSpanel]'
	if (!_rc) local panelvar `_dta[_TSpanel]'

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
	

	* For this, _iv_parse would have been useful, although I don't want to do factor expansions when parsing
	if ("`model'"=="iv") {
		* get part before parentheses
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

		* get part in parentheses
		gettoken right 0 : 0 ,bind match(parens)
		Assert trim(`"`0'"')=="" , msg("error: remaining argument: `0'")

		* now parse part in parentheses
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
		Assert !_rc , msg("error: -`ivsuite'- not installed, please run {stata ssc install `ivsuite'} or change the option -ivsuite-")
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

* Stages
	assert "`model'"!="" // just to be sure this goes after `model' is set
	local iv_stage iv
	local stages : list stages - iv_stage
	local valid_stages ols first acid reduced
	local wrong_stages : list stages - valid_stages
	Assert "`wrong_stages'"=="", msg("Error, invalid stages(): `wrong_stages'")
	if ("`stages'"!="") {
		Assert "`model'"=="iv", msg("Error, stages() only valid with an IV regression")
		local stages `stages' `iv_stage' // Put -iv- *last* (so it does the -restore-; note that we don't need it first to trim MVs b/c that's done earlier)
		Assert "`avge'"=="", msg("Error, avge not allowed with stages()")
	}
	else {
		local stages none // So we can loop over stages
	}

* Add back constants (place this *after* we define `model')
	local addconstant = ("`constant'"!="noconstant") & !("`model'"=="iv" & "`ivsuite'"=="ivreg2") // also see below

* Parse VCE options:
	
	* Note: bw=1 *usually* means just do HC instead of HAC
	* BUGBUG: It is not correct to ignore the case with "bw(1) kernel(Truncated)"
	* but it's too messy to add -if-s everywhere just for this rare case (see also Mark Schaffer's email)

	local 0 `vce'
	syntax [anything(id="VCE type")] , [bw(integer 1)] [KERnel(string)] [dkraay(integer 1)] [kiefer] [suite(string)]
	if ("`anything'"=="") local anything unadjusted
	Assert `bw'>0, msg("VCE bandwidth must be a positive integer")
	gettoken vcetype clustervars : anything
	* Expand variable abbreviations; but this adds unwanted i. prefixes
	if ("`clustervars'"!="") {
		fvunab clustervars : `clustervars'
		local clustervars : subinstr local clustervars "i." "", all
	}

	* vcetype abbreviations:
	if (substr("`vcetype'",1,3)=="ols") local vcetype unadjusted
	if (substr("`vcetype'",1,2)=="un") local vcetype unadjusted
	if (substr("`vcetype'",1,1)=="r") local vcetype robust
	if (substr("`vcetype'",1,2)=="cl") local vcetype cluster
	if ("`vcetype'"=="conventional") local vcetype unadjusted // Conventional is the name given in e.g. xtreg
	Assert strpos("`vcetype'",",")==0, msg("Unexpected contents of VCE: <`vcetype'> has a comma")

	* Sanity checks on vcetype
	if ("`vcetype'"=="" & "`backupweight'"=="pweight") local vcetype robust
	Assert !("`vcetype'"=="unadjusted" & "`backupweight'"=="pweight"), msg("pweights do not work with unadjusted errors, use a different vce()")
	if ("`vcetype'"=="") local vcetype unadjusted
	Assert inlist("`vcetype'", "unadjusted", "robust", "cluster"), msg("VCE type not supported: `vcetype'")

	* Cluster vars
	local num_clusters : word count `clustervars'
	Assert inlist( (`num_clusters'>0) + ("`vcetype'"=="cluster") , 0 , 2), msg("Can't specify cluster without clustervars and viceversa") // XOR

	* VCE Suite
	local vcesuite `suite'
	if ("`vcesuite'"=="") local vcesuite default
	if ("`vcesuite'"=="default") {
		if (`bw'>1 | `dkraay'>1 | "`kiefer'"!="" | "`kernel'"!="") {
			local vcesuite avar
		}
		else if (`num_clusters'>1) {
			local vcesuite mwc
		}
	}

	Assert inlist("`vcesuite'", "default", "mwc", "avar"), msg("Wrong vce suite: `vcesuite'")
	if (inlist("`vcesuite'", "avar", "mwc")) local addconstant 0 // The constant messes up the VCV

	if ("`vcesuite'"=="mwc") {
		cap findfile tuples.ado
		Assert !_rc , msg("error: -tuples- not installed, please run {stata ssc install tuples} to estimate multi-way clusters.")
	}
	
	if ("`vcesuite'"=="avar") {
		cap findfile `vcesuite'.ado
		Assert !_rc , msg("error: -`vcesuite'- not installed, please run {stata ssc install `vcesuite'} or change the option -vcesuite-")
	}

	* Some combinations are not coded
	Assert !("`ivsuite'"=="ivregress" & (`num_clusters'>1 | `bw'>1 | `dkraay'>1 | "`kiefer'"!="" | "`kernel'"!="") ), msg("option vce(`vce') incompatible with ivregress")
	Assert !("`ivsuite'"=="ivreg2" & (`num_clusters'>2) ), msg("ivreg2 doesn't allow more than two cluster variables")
	Assert !("`model'"=="ols" & "`vcesuite'"=="avar" & (`num_clusters'>2) ), msg("avar doesn't allow more than two cluster variables")
	Assert !("`model'"=="ols" & "`vcesuite'"=="default" & (`bw'>1 | `dkraay'>1 | "`kiefer'"!="" | "`kernel'"!="") ), msg("to use those vce options you need to use -avar- as the vce suite")
	if (`num_clusters'>0) local temp_clustervars " <CLUSTERVARS>"
	if (`bw'==1 & `dkraay'==1 & "`kernel'"!="") local kernel // No point in setting kernel here 
	if (`bw'>1 | "`kernel'"!="") local vceextra `vceextra' bw(`bw') 
	if (`dkraay'>1) local vceextra `vceextra' dkraay(`dkraay') 
	if ("`kiefer'"!="") local vceextra `vceextra' kiefer 
	if ("`kernel'"!="") local vceextra `vceextra' kernel(`kernel')
	if ("`vceextra'"!="") local vceextra , `vceextra'
	local vceoption "`vcetype'`temp_clustervars'`vceextra'" // this excludes "vce(", only has the contents

* DoF Adjustments
	if ("`dofadjustments'"=="") local dofadjustments all
	local 0 , `dofadjustments'
	syntax, [ALL NONE] [PAIRwise FIRSTpair] [CLusters] [CONTinuous]
	opts_exclusive "`all' `none'" dofadjustments
	opts_exclusive "`pairwise' `firstpair'" dofadjustments
	if ("`none'"!="") {
		Assert "`pairwise'`firstpair'`clusters'`continuous'"=="", msg("option {bf:dofadjustments()} invalid; {bf:none} not allowed with other alternatives")
		local dofadjustments
	}
	if ("`all'"!="") {
		Assert "`pairwise'`firstpair'`clusters'`continuous'"=="", msg("option {bf:dofadjustments()} invalid; {bf:all} not allowed with other alternatives")
		local dofadjustments pairwise clusters continuous
	}
	else {
		local dofadjustments `pairwise' `firstpair' `clusters' `continuous'
	}

* Mobility groups
	if ("`group'"!="") conf new var `group'

* IV options
	if ("`small'"!="") di in ye "(note: reghdfe will always use the option -small-, no need to specify it)"

	Assert ("`gmm2s'`liml'`cue'"==""), msg("options gmm2s/liml/cue not allowed")
	
	if ("`model'"=="iv") {
		local savefirst = ("`savefirst'"!="")
		local first = ("`first'"!="")
		if (`savefirst') Assert `first', msg("Option -savefirst- requires -first-")
	}

} // End of !`savingcache'

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
		Assert !_rc , msg("error: -parallel- not installed, please run {stata ssc install parallel}")
	}
	local opt_list tolerance maxiterations check accelerate ///
		bad_loop_threshold stuck_threshold pause_length accel_freq accel_start
	foreach opt of local opt_list {
		if ("``opt''"!="") local maximize_options `maximize_options' `opt'(``opt'')
	}

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
		ivsuite showraw ///
		depvar indepvars endogvars instruments savefirst first ///
		vceoption vcetype vcesuite vceextra num_clusters clustervars bw kernel dkraay /// vceextra
		dofadjustments ///
		if in group check fast nested fe_format ///
		tolerance maxiterations accelerate maximize_options ///
		subcmd suboptions ///
		absorb avge excludeself ///
		timevar panelvar basevars ///
		addconstant ///
		weight weightvar exp weightexp /// type of weight (fw,aw,pw), weight var., and full expr. ([fw=n])
		cores savingcache usecache over ///
		postestimation notes stages
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
