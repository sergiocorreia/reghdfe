capture program drop ParseIV
program define ParseIV, sclass
	syntax anything(id="varlist" name=0 equalok), [ ///
		estimator(string) ivsuite(string) ///
		savefirst first showraw vceunadjusted small]

	* Parses varlist: depvar indepvars [(endogvars = instruments)]
		* depvar: dependent variable
		* indepvars: included exogenous regressors
		* endogvars: included endogenous regressors
		* instruments: excluded exogenous regressors

	* Model: OLS or IV-type?
	local model ols
	foreach _ of local 0 {
		if (substr(`"`_'"', 1, 1)=="(") {
			local model iv
			continue, break
		}
	}

	* IV Suite
	if ("`model'"=="iv") {
		if ("`ivsuite'"=="") local ivsuite ivreg2 // Set default
		Assert inlist("`ivsuite'","ivreg2","ivregress") , msg("error: wrong IV routine (`ivsuite'), valid options are -	ivreg2- and -ivregress-")
		cap findfile `ivsuite'.ado
		Assert !_rc , msg("error: -`ivsuite'- not installed, please run {stata ssc install `ivsuite'} or change the option 	-ivsuite-")

		local savefirst = ("`savefirst'"!="")
		local first = ("`first'"!="")
		if (`savefirst') Assert `first', msg("Option -savefirst- requires -first-")
		local subcmd `ivsuite'
	}
	else {
		local savefirst
		local first
		local subcmd regress
	}

	* Estimator
	if ("`estimator'"=="" & "`model'"=="iv") local estimator 2sls // Set default
	if ("`estimator'"!="") {
		Assert "`model'"=="iv", ///
			msg("reghdfe error: estimator() requires an instrumental-variable regression")
		if (substr("`estimator'", 1, 3)=="gmm") local estimator gmm2s
		Assert inlist("`estimator'", "2sls", "gmm2s", "liml", "cue"), ///
			msg("reghdfe error: invalid estimator `estimator'")
		if ("`estimator'"=="cue") Assert "`ivsuite'"=="ivreg2", ///
			msg("reghdfe error: estimator `estimator' only available with the ivreg2 command, not ivregress")
	}

	* For this, _iv_parse would have been useful, but I don't want to do factor expansions when parsing
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
		local endogvars `varlist'

		Assert "`instruments'"!="", msg("iv: instruments required")
		local 0 `instruments'
		syntax varlist(fv ts numeric)
		local instruments `varlist'
		
		local 0 `left' // So OLS part can handle it
	}

* OLS varlist
	syntax varlist(fv ts numeric)
	gettoken depvar indepvars : varlist
	_fv_check_depvar `depvar'

* Extract format of depvar so we can format FEs like this
	fvrevar `depvar', list
	local fe_format : format `r(varlist)' // The format of the FEs that will be saved

* Variables shouldn't be repeated
* This is not perfect (e.g. doesn't deal with "x1-x10") but still helpful
	local allvars `depvar' `indepvars' `endogvars' `instruments'
	local dupvars : list dups allvars
	Assert "`dupvars'"=="", msg("error: there are repeated variables: <`dupvars'>")

* More IV options
	if ("`small'"!="") di in ye "(note: reghdfe will always use the option -small-, no need to specify it)"

* Parse -showraw- : shows raw output of called subcommand (e.g. ivreg2)
	local showraw = ("`showraw'"!="")

* Parse -unadjusted-
	* If true, will use wmatrix(...) vce(unadjusted) instead of the default of setting vce contents equal to wmatrix
	* This basically undoes the extra adjustment that ivregress does, so it's comparable with ivreg2
	*
	* Note: Cannot match exactly the -ivregress- results without vceunadjusted (see test-gmm.do)
	* Thus, I will set this to true ALWAYS
	local vceunadjusted = 1 // ("`vceunadjusted'"!="")

* Get base variables of time and factor variables (e.g. i.foo L(1/3).bar -> foo bar)
	foreach vars in depvar indepvars endogvars instruments {
		if ("``vars''"!="") {
			fvrevar ``vars'' , list
			local basevars `basevars' `r(varlist)'
		}
	}

	local keys subcmd model ivsuite estimator depvar indepvars endogvars instruments fe_format ///
		savefirst first showraw vceunadjusted basevars
	foreach key of local keys {
		sreturn local `key' ``key''
	}
end 
