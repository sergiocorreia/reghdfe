cap pr drop Wrapper_ivregress
program define Wrapper_ivregress, eclass
	syntax , depvar(varname) endogvars(varlist) instruments(varlist) ///
		[indepvars(varlist) avgevars(varlist)] ///
		original_depvar(string) original_endogvars(string) original_instruments(string) ///
		[original_indepvars(string) avge_targets(string)] ///
		estimator(string) vceoption(string asis) KK(integer) [weightexp(string)] ///
		[SUBOPTions(string)] [*] // [*] are ignored!

	mata: st_local("vars", strtrim(stritrim( "`depvar' `indepvars' `avgevars' (`endogvars'=`instruments')" )) )
	if ("`estimator'"=="gmm") local vceoption = "`vceoption' " + subinstr("`vceoption'", "vce(", "wmatrix(", .)
	
	* Note: the call to -ivregress- could be optimized.
	* EG: -ivregress- calls ereturn post .. ESAMPLE(..) but we overwrite the esample and its SLOW
	* But it's a 1700 line program so let's not worry about it
	*profiler on

	local subcmd ivregress `estimator' `vars' `weightexp', `vceoption' small `suboptions'
	Debug, level(3) msg("Subcommand: " in ye "`subcmd'")
	local noise qui
	if (strpos("`options'", "first")>0) local noise noi
	`noise' `subcmd'
	
	*profiler off
	*profiler report
	
	* Fix DoF if needed
	local N = e(N)
	local K = e(df_m)
	local WrongDoF = `N' - 1 - `K'
	local CorrectDoF = `WrongDoF' - `kk'
	Assert !missing(`CorrectDoF')
	if ("`estimator'"!="gmm" | 1) {
		tempname V
		matrix `V' = e(V) * (`WrongDoF' / `CorrectDoF')
		ereturn repost V=`V'
	}
	ereturn scalar df_r = `CorrectDoF'

	* ereturns specific to this command
	mata: st_local("original_vars", strtrim(stritrim( "`original_depvar' `original_indepvars' `avge_targets' `original_absvars' (`original_endogvars'=`original_instruments')" )) )
	ereturn local alternative_cmd ivregress `estimator' `original_vars', small `vceoption' `options'
	if ("`estimator'"!="gmm") ereturn scalar F = e(F) * `CorrectDoF' / `WrongDoF'

end
