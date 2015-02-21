/* Notes:
- For -cluster- I need to run two regressions as in xtreg; the first one is to get the df_m

- El FTEST es distinto entre areg/reg y xtreg porque xtreg hace un ajuste extra
para pasar de areg a xtreg, multiplicar el F por Q^2
Donde Q =  (e(N) - e(rank)) / (e(N) - e(rank) - e(df_a))
Es decir, en vez de dividir entre N-K-KK, me basta con dividir entre N-K
Asi que me bastaria usar -test- despues de correr la regresion y deberia salir igual que el FTEST ajustado del areg!!!
(tambien igual al del xreg pero eso es mas limitante , aunq igual probar para 1 HDFE creando t=_n a nivel del ID1)
*/
*/

cap pr drop Wrapper_regress
program define Wrapper_regress, eclass
	syntax , depvar(varname) [indepvars(varlist) avgevars(varlist)] ///
		original_absvars(string) original_depvar(string) [original_indepvars(string) avge_targets(string)] ///
		vceoption(string asis) kk(integer) vcetype(string) [weightexp(string)] ///
		addconstant(integer) ///
		[SUBOPTions(string)] [*] // [*] are ignored!

	if ("`options'"!="") Debug, level(3) msg("(ignored options: `options')")

	local vceoption = regexr("`vceoption'", "vce\( *unadjusted *\)", "vce(ols)")
	mata: st_local("vars", strtrim(stritrim( "`depvar' `indepvars' `avgevars'" )) ) // Just for esthetic purposes

* Hide constant
	if (!`addconstant') {
		local nocons noconstant
		local kk = `kk' + 1
	}

* Run regression just to compute true DoF
	local subcmd _regress `vars' `weightexp', noheader notable `suboptions'
	Debug, level(3) msg("Subcommand: " in ye "`subcmd'")
	qui `subcmd'
	local N = e(N)
	local K = e(df_m) // Should also be equal to e(rank)+1
	*** scalar `sse' = e(rss)
	local WrongDoF = `N' - `addconstant' - `K'
	local CorrectDoF = `WrongDoF' - `kk' // kk = Absorbed DoF
	Assert !missing(`CorrectDoF')
	
* Now run intended regression and fix VCV
	qui regress `vars' `weightexp', `vceoption' noheader notable `suboptions' `nocons'
	* Fix DoF
	tempname V
	cap matrix `V' = e(V) * (`WrongDoF' / `CorrectDoF')

	* Avoid corner case error when all the RHS vars are collinear with the FEs
	if (`K'>0) {
		cap ereturn repost V=`V' // Else the fix would create MVs and we can't post
		Assert inlist(_rc,0,504), msg("error `=_rc' when adjusting the VCV")
	}
	else {
		ereturn scalar rank = 1 // Set e(rank)==1 when e(df_m)=0 , due to the constant
		* (will not be completely correct if model is already demeaned?)
	}
	
	*** if ("`vcetype'"!="cluster") ereturn scalar rank = e(rank) + `kk'

* ereturns specific to this command
	ereturn scalar df_r = max(`CorrectDoF', 0)
	mata: st_local("original_vars", strtrim(stritrim( "`original_depvar' `original_indepvars' `avge_targets' `original_absvars'" )) )
	ereturn local alternative_cmd regress `original_vars', `vceoption' `options'

	if ("`vcetype'"!="cluster") { // ("`vcetype'"=="unadjusted")
		ereturn scalar F = e(F) * `CorrectDoF' / `WrongDoF'
		if missing(e(F)) di as error "WARNING! Missing FStat"
	}


	local run_test = ("`vcetype'"=="cluster") | ( e(df_m)+1!=e(rank) )
	if (`run_test') {
		if ("`vcetype'"!="cluster") {
			Debug, level(0) msg("Note: equality df_m+1==rank failed (is there a collinear variable in the RHS?), running -test- to get correct values")
		}
		return clear
		if (`K'>0) {
			qui test `indepvars' `avge' // Wald test
			ereturn scalar F = r(F)
			ereturn scalar df_m = r(df)
			ereturn scalar rank = r(df)+1 // Add constant
		}
		else {
			ereturn scalar F = 0
			ereturn scalar df_m = 0
			ereturn scalar rank = 1
		}
	}

* Fstat
	* _U: Unrestricted, _R: Restricted
	* FStat = (RSS_R - RSS_U) / RSS * (N-K) / q
	*       = (R2_U - R2_R) / (1 - R2_U) * DoF_U / (DoF_R - DoF_U)
	Assert e(df_m)+1==e(rank) , rc(0) /// rc(322)
		msg("Error: expected e(df_m)+1==e(rank), got (`=`e(df_m)'+1'!=`e(rank)')")
end

* Cluster notes (see stata PDFs):
* We don't really have "N" indep observations but "M" (aka `N_clust') superobservations,
* and we are replacing (N-K) DoF with (M-1) (used when computing the T and F tests)
		
* For the VCV matrix, the multiplier (small sample adjustement) is q := (N-1)/(N-K) * M / (M-1)
* Notice that if every obs is its own cluster, M=N and q = N/(N-K) (the usual multiplier for -ols- and -robust-)
		
* Also, if one of the absorbed FEs is nested within the cluster variable, then we don't need to include that variable in K
* (this is the adjustment that xtreg makes that areg doesn't)
