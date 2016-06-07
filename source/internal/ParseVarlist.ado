cap pr drop ParseVarlist
pr ParseVarlist, sclass
	sreturn clear
	syntax anything(id="varlist" name=0 equalok), ///
		[estimator(string) ivsuite(string)]

	/* 	Varlist syntax: depvar indepvars [(endogvars = instruments)]
		depvar		: 	dependent variable
		indepvars	: 	included exogenous regressors
		endogvars	: 	included endogenous regressors
		instruments	: 	excluded exogenous regressors

	Note: run this AFTER -_fvunab- (which looks for string vars, etc.) */

	ParseDepvar `0' // s(depvar) s(fe_format) ; s(rest)
	ParseIndepvars `s(rest)' // s(indepvars) ; s(parens) s(rest)
	_assert "`s(rest)'" == "", ///
		msg("couldn't parse the end of the varlist: <`s(rest)'>")
	ParseEndogAndInstruments `s(parens)' // s(endogvars) s(instruments)
	sreturn local parens // clear" it"

	* Sanity checks
	if ("`s(instruments)'" != "") {
		if ("`ivsuite'"=="") local ivsuite ivreg2 // Set default
		_assert inlist("`ivsuite'","ivreg2","ivregress") , ///
			msg("error: wrong IV routine (`ivsuite'), valid options are -ivreg2- and -ivregress-")
		cap findfile `ivsuite'.ado
		_assert !_rc , msg("error: -`ivsuite'- not installed, please run {stata ssc install `ivsuite'} or change the option 	-ivsuite-")
		local subcmd `ivsuite'

		if ("`estimator'"=="") local estimator 2sls // Set default
		if (substr("`estimator'", 1, 3)=="gmm") local estimator gmm2s
		_assert inlist("`estimator'", "2sls", "gmm2s", "liml", "cue"), ///
			msg("reghdfe error: invalid estimator `estimator'")
		if ("`estimator'"=="cue") Assert "`ivsuite'"=="ivreg2", ///
			msg("reghdfe error: estimator `estimator' only available with the ivreg2 command, not ivregress")
		if ("`estimator'"=="cue") di as text "(WARNING: -cue- estimator is not exact, see help file)"
	}
	else {
		local subcmd regress
		_assert "`estimator'"=="", msg("estimator() requires an instrumental-variable regression")
	}

	sreturn local ivsuite `ivsuite'
	sreturn local subcmd `subcmd'
	sreturn local estimator `estimator'
end

cap pr drop ParseDepvar
pr ParseDepvar, sclass
	gettoken depvar 0 : 0, bind
	fvexpand `depvar'
	local depvar `r(varlist)'
	local n : word count `depvar'
	_assert (`n'==1), msg("more than one depvar specified: `depvar'")
	sreturn local depvar `depvar'
	sreturn local rest `0'

* Extract format of depvar so we can format FEs like this
	fvrevar `depvar', list
	local fe_format : format `r(varlist)' // The format of the FEs that will be saved
	sreturn local fe_format `fe_format'
end

cap pr drop ParseIndepvars
pr ParseIndepvars, sclass
	while ("`0'" != "") {
		gettoken _ 0 : 0, bind match(parens)
		if ("`parens'" == "") {
			local indepvars `indepvars' `_'
		}
		else {
			continue, break
		}
	}
	sreturn local indepvars `indepvars'
	if ("`parens'" != "") sreturn local parens "`_'"
	sreturn local rest `0'
end

cap pr drop ParseEndogAndInstruments
pr ParseEndogAndInstruments, sclass
	if ("`0'" == "") exit
	gettoken _ 0 : 0, bind parse("=")
	if ("`_'" != "=") {
		sreturn local endogvars `_'
		gettoken equalsign 0 : 0, bind parse("=")
	}
	sreturn local instruments `0'
end
