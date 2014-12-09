cap pr drop Wrapper_ivreg2
program define Wrapper_ivreg2, eclass
	syntax , depvar(varname) endogvars(varlist) instruments(varlist) ///
		[indepvars(varlist) avgevars(varlist)] ///
		original_depvar(string) original_endogvars(string) original_instruments(string) ///
		[original_indepvars(string) avge_targets(string)] ///
		[original_absvars(string) avge_targets] ///
		estimator(string) vceoption(string asis) KK(integer) ///
		[SHOWRAW(integer 0) dofminus(string)] first(integer) [weightexp(string)] ///
		[SUBOPTions(string)] [*] // [*] are ignored!
	if ("`options'"!="") Debug, level(3) msg("(ignored options: `options')")
	if (`c(version)'>=12) local hidden hidden
	
	* Disable some options
	local 0 , `suboptions'
	syntax , [SAVEFPrefix(name)] [*] // Will ignore SAVEFPREFIX
	local suboptions `options'

	if ("`estimator'"=="2sls") local estimator

	if strpos("`vceoption'","unadj") {
		local vceoption
	}
	else if strpos("`vceoption'","robust")>0 {
		local vceoption "robust" 
	}
	else {
		local vceoption = substr(trim("`vceoption'"), 5, strlen("`vceoption'")-5)
		gettoken vcefirst vceoption : vceoption
		local vceoption = trim("`vceoption'")
		local vceoption `vcefirst'(`vceoption')
	}
	
	mata: st_local("vars", strtrim(stritrim( "`depvar' `indepvars' `avgevars' (`endogvars'=`instruments')" )) )
	
	if (`first') {
		local firstoption "first savefirst"
	}
	Assert inlist("`dofminus'","dofminus","sdofminus")

	local subcmd ivreg2 `vars' `weightexp', `estimator' `vceoption' `firstoption' small `dofminus'(`kk') `suboptions'
	Debug, level(3) msg(_n "call to subcommand: " _n as result "`subcmd'")
	local noise = cond(`showraw', "noi", "qui")
	`noise' `subcmd'
	if ("`noise'"=="noi") di in red "{hline 64}" _n "{hline 64}"

	if !missing(e(ecollin)) {
		di as error "endogenous covariate <`e(ecollin)'> was perfectly predicted by the instruments!"
		error 2000
	}

	if (`first') {
		ereturn `hidden' local first_prefix = "_ivreg2_"
		ereturn `hidden' local ivreg2_firsteqs = e(firsteqs)
		ereturn local firsteqs
	}

	foreach cat in exexog insts instd {
		FixVarnames `e(`cat')'
		ereturn local `cat' = "`r(newnames)'"
	}

	if (`first') {
		* May be a problem if we ran out of space for storing estimates
		local ivreg2_firsteqs "`e(ivreg2_firsteqs)'"
		tempname hold
		estimates store `hold' , nocopy
		foreach fs_eqn in `ivreg2_firsteqs' {
			qui estimates restore `fs_eqn'
			FixVarnames `e(depvar)'
			ereturn local depvar = r(prettynames)
			FixVarnames `e(inexog)'
			ereturn local inexog = r(prettynames)

			tempname b
			matrix `b' = e(b)
			local backup_colnames : colnames `b'
			FixVarnames `backup_colnames'
			matrix colnames `b' = `r(prettynames)' // newnames? prettynames?
			ereturn repost b=`b', rename

			estimates store `fs_eqn', nocopy
		}
		qui estimates restore `hold'
		estimates drop `hold'
	}

	* ereturns specific to this command
	mata: st_local("original_vars", strtrim(stritrim( "`original_depvar' `original_indepvars' `avge_targets' `original_absvars' (`original_endogvars'=`original_instruments')" )) )
	ereturn local alternative_cmd ivreg2 `original_vars', small `vceoption' `options' `estimator'
	***if ("`estimator'"!="gmm") ereturn scalar F = e(F) * `CorrectDoF' / `WrongDoF'

end
