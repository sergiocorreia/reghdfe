cap pr drop Wrapper_ivregress
program define Wrapper_ivregress, eclass
	syntax , depvar(varname) endogvars(varlist) instruments(varlist) ///
		[indepvars(varlist)] ///
		vceoption(string asis) ///
		KK(integer) ///
		[weightexp(string)] ///
		SHOWRAW(integer) first(integer) vceunadjusted(integer) ///
		[ESTimator(string) TWICErobust(string)] ///
		[SUBOPTions(string)] [*] // [*] are ignored!

	if ("`options'"!="") Debug, level(3) msg("(ignored options: `options')")
	mata: st_local("vars", strtrim(stritrim( "`depvar' `indepvars' (`endogvars'=`instruments')" )) )
	if (`c(version)'>=12) local hidden hidden

	local opt_estimator = cond("`estimator'"=="gmm2s", "gmm", "`estimator'")

	* Convert -vceoption- to what -ivreg2- expects
	local 0 `vceoption'
	syntax namelist(max=2)
	gettoken vceoption clustervars : namelist
	local clustervars `clustervars' // Trim
	Assert inlist("`vceoption'", "unadjusted", "robust", "cluster")
	if ("`clustervars'"!="") local vceoption `vceoption' `clustervars'
	local vceoption "vce(`vceoption')"

	if ("`estimator'"=="gmm2s") {
		local wmatrix : subinstr local vceoption "vce(" "wmatrix("
		if ("`twicerobust'"=="") {
			local vceoption = cond(`vceunadjusted', "vce(unadjusted)", "")			
		}
	}
	
	* Note: the call to -ivregress- could be optimized.
	* EG: -ivregress- calls ereturn post .. ESAMPLE(..) but we overwrite the esample and its SLOW
	* But it's a 1700 line program so let's not worry about it

* Show first stage
	if (`first') {
		local firstoption "first"
	}

* Subcmd
	local subcmd ivregress `opt_estimator' `vars' `weightexp', `wmatrix' `vceoption' small noconstant `firstoption' `suboptions'
	Debug, level(3) msg("Subcommand: " in ye "`subcmd'")
	local noise = cond(`showraw', "noi", "qui")
	`noise' `subcmd'
	qui test `indepvars' `endogvars' // Wald test
	ereturn scalar F = r(F)

	
	* Fix DoF if needed
	local N = e(N)
	local K = e(df_m)
	local WrongDoF = `N' - `K'
	local CorrectDoF = `WrongDoF' - `kk'
	Assert !missing(`CorrectDoF')

	* We should have used M/M-1 instead of N/N-1, but we are making ivregress to do the wrong thing by using vce(unadjusted) (which makes it fit with ivreg2)
	local q 1
	if ("`estimator'"=="gmm2s" & "`clustervars'"!="") {
		local N = e(N)
		tempvar group
		GenerateID `clustervars', gen(`group')
		su `group', mean
		drop `group'
		local M = r(max) // N_clust
		local q = ( `M' / (`M' - 1)) / ( `N' / (`N' - 1) ) // multiply correct, divide prev wrong one
		ereturn scalar df_r = `M' - 1
	}

	tempname V
	matrix `V' = e(V) * (`WrongDoF' / `CorrectDoF') * `q'
	ereturn repost V=`V'
	
	if ("`clustervars'"=="") ereturn scalar df_r = `CorrectDoF'

	* ereturns specific to this command
	ereturn scalar F = e(F) * `CorrectDoF' / `WrongDoF'

	ereturn scalar tss = e(mss) + e(rss) // ivreg2 doesn't report e(tss)
	ereturn `hidden' scalar unclustered_df_r = `CorrectDoF' // Used later in R2 adj
end
