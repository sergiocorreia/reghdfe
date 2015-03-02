cap pr drop Wrapper_regress
program define Wrapper_regress, eclass
	syntax , depvar(varname) [indepvars(varlist) avgevars(varlist)] ///
		original_absvars(string) original_depvar(string) [original_indepvars(string) avge_targets(string)] ///
		vceoption(string asis) ///
		kk(integer) ///
		[weightexp(string)] ///
		addconstant(integer) ///
		[SUBOPTions(string)] [*] // [*] are ignored!

	if ("`options'"!="") Debug, level(3) msg("(ignored options: `options')")
	mata: st_local("vars", strtrim(stritrim( "`depvar' `indepvars' `avgevars'" )) ) // Just for esthetic purposes
	if (`c(version)'>=12) local hidden hidden

* Convert unadjusted into ols so -regress- can handle it
	local vceoption : subinstr local vceoption "unadjusted" "ols"
	local vceoption "vce(`vceoption')"

* Hide constant
	if (!`addconstant') {
		local nocons noconstant
		local kk = `kk' + 1
	}

* Obtain e(df_m)
	_rmcoll `vars' `weightexp', forcedrop
	local varlist = r(varlist)
	local df_m : list sizeof varlist
	
	

	local subcmd regress `vars' `weightexp', noheader notable `suboptions'
	Debug, level(3) msg("Subcommand: " in ye "`subcmd'")
	qui `subcmd'

	local K = e(df_m)
	local WrongDoF = `N' - `addconstant' - `K'
	assert `WrongDoF'==e(df_r) // DoF should equal N-K-1
	
	local CorrectDoF = `WrongDoF' - `kk' // kk = Absorbed DoF
	Assert !missing(`CorrectDoF')
	
* Now run intended regression and fix VCV
	local new_df_r = max(`CorrectDoF', 0)
	local subcmd regress `vars' `weightexp', `vceoption' `suboptions' `nocons' dof(`new_df_r') // noheader notable
	Debug, level(3) msg("Subcommand: " in ye "`subcmd'")
	qui `subcmd'

	* When vcetype -unadjusted-, setting the dof() option is enough
	* But with -robust-, 
	
	* Fix DoF
	tempname V

	cap matrix `V' = e(V) * (`WrongDoF' / `CorrectDoF')
	ereturn scalar F = e(F) / (`WrongDoF' / `CorrectDoF')
	ereturn scalar df_r = `new_df_r'
	if ("`vcetype'"=="cluster")  ereturn scalar df_r = e(N_clust) - 1
	* unadj: dof(new_df_r) ; no ajustar V luego
	* robust: si anhado el ajuste a V() se arregla V ; pero ademas tengo que arreglar el df_r y el F porque EL PUTO REGRESS AJUSTA MAL EL DOF!!!!!!!!!!
	di as error e(df_r)
	di as error `new_df_r'

	* Avoid corner case error when all the RHS vars are collinear with the FEs
	if (`K'>0) {
		cap ereturn repost V=`V' // Else the fix would create MVs and we can't post
		Assert inlist(_rc,0,504), msg("error `=_rc' when adjusting the VCV")
	}
	else {
		ereturn scalar rank = 1 // Set e(rank)==1 when e(df_m)=0 , due to the constant
		* (will not be completely correct if model is already demeaned?)
	}

* ereturns specific to this command
	mata: st_local("original_vars", strtrim(stritrim( "`original_depvar' `original_indepvars' `avge_targets' `original_absvars'" )) )

	if ("`vcetype'"!="cluster") { // ("`vcetype'"=="unadjusted")
		ereturn scalar F = e(F) // * `CorrectDoF' / `WrongDoF'
		if missing(e(F)) di as error "WARNING! Missing FStat"
	}
	else {
		ereturn `hidden' scalar unclustered_df_r = `CorrectDoF'
		assert e(N_clust)<.
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

	ereturn scalar tss = e(mss) + e(rss) // Regress doesn't report e(tss)

* Fstat
	* _U: Unrestricted, _R: Restricted
	* FStat = (RSS_R - RSS_U) / RSS * (N-K) / q
	*       = (R2_U - R2_R) / (1 - R2_U) * DoF_U / (DoF_R - DoF_U)
	Assert e(df_m)+1==e(rank) , rc(0) msg("Error: expected e(df_m)+1==e(rank), got (`=`e(df_m)'+1'!=`e(rank)')")
end

* Cluster notes (see stata PDFs):
* We don't really have "N" indep observations but "M" (aka `N_clust') superobservations,
* and we are replacing (N-K) DoF with (M-1) (used when computing the T and F tests)
		
* For the VCV matrix, the multiplier (small sample adjustement) is q := (N-1)/(N-K) * M / (M-1)
* Notice that if every obs is its own cluster, M=N and q = N/(N-K) (the usual multiplier for -ols- and -robust-)
		
* Also, if one of the absorbed FEs is nested within the cluster variable, then we don't need to include that variable in K
* (this is the adjustment that xtreg makes that areg doesn't)
