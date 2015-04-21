cap pr drop Wrapper_regress
program define Wrapper_regress, eclass
	syntax , depvar(varname) [indepvars(varlist) avgevars(varlist)] ///
		vceoption(string asis)  ///
		kk(integer) ///
		[weightexp(string)] ///
		[SUBOPTions(string)] [*] // [*] are ignored!
	
	if ("`options'"!="") Debug, level(3) msg("(ignored options: `options')")
	mata: st_local("vars", strtrim(stritrim( "`depvar' `indepvars' `avgevars'" )) ) // Just for aesthetic purposes
	if (`c(version)'>=12) local hidden hidden

* Convert -vceoption- to what -regress- expects
	gettoken vcetype clustervars : vceoption
	local clustervars `clustervars' // Trim
	local vceoption : subinstr local vceoption "unadjusted" "ols"
	local vceoption "vce(`vceoption')"

* Note: the dof() option of regress is *useless* with robust errors,
* and overriding e(df_r) is also useless because -test- ignores it,
* so we have to go all the way and do a -post- from scratch

* Obtain K so we can obtain DoF = N - K - kk
* This is already done by regress EXCEPT when clustering
* (but we still need the unclustered version for r2_a, etc.)
	_rmcoll `indepvars' `avgevars' `weightexp', forcedrop
	local varlist = r(varlist)
	if ("`varlist'"==".") local varlist
	local K : list sizeof varlist

* Run -regress-
	local subcmd regress `vars' `weightexp', `vceoption' `suboptions' noconstant noheader notable
	Debug, level(3) msg("Subcommand: " in ye "`subcmd'")
	qui `subcmd'
	
	local N = e(N) // We couldn't just use c(N) due to possible frequency weights
	local WrongDoF = `N' - `K'
	if ("`vcetype'"!="cluster" & e(df_r)!=`WrongDoF') {
		local difference = `WrongDoF' - e(df_r)
		local NewDFM = e(df_m) - `difference'	
		di as result "(warning: regress returned e(df_r)==`e(df_r)', but we expected it to be `WrongDoF')"
		Assert e(df_m)>=0, msg("try removing collinear regressors or setting a higher tol()")
		di as result "(workaround: we will set e(df_m)=`NewDFM' instead of `e(df_m)')"
	}
	local CorrectDoF = `WrongDoF' - `kk' // kk = Absorbed DoF

* Store results for the -ereturn post-
	tempname b V
	matrix `b' = e(b)
	matrix `V' = e(V)
	local N = e(N)
	local marginsok = e(marginsok)
	local rmse = e(rmse)
	local rss = e(rss)
	local tss = e(mss) + e(rss) // Regress doesn't report e(tss)
	local N_clust = e(N_clust)

	local predict = e(predict)
	local cmd = e(cmd)
	local cmdline = e(cmdline)
	local title = e(title)

	* Fix V
	if (`K'>0) matrix `V' = `V' * (`WrongDoF' / `CorrectDoF')

	* DoF
	if ("`vcetype'"=="cluster") {
		Assert e(df_r) == e(N_clust) - 1
		Assert e(N_clust) > `K', msg("insufficient observations (N_clust=`e(N_clust)', K=`K')") rc(2001)
	}
	local df_r = cond( "`vcetype'"=="cluster" , e(df_r) , max( `CorrectDoF' , 0 ) )

	capture ereturn post `b' `V' `weightexp', dep(`depvar') obs(`N') dof(`df_r') properties(b V)
	local rc = _rc
	Assert inlist(_rc,0,504), msg("error `=_rc' when adjusting the VCV") // 504 = Matrix has MVs
	Assert `rc'==0, msg("Error: estimated variance-covariance matrix has missing values")
	ereturn local marginsok = "`marginsok'"
	ereturn local predict = "`predict'"
	ereturn local cmd = "`cmd'"
	ereturn local cmdline = "`cmdline'"
	ereturn local title = "`title'"
	ereturn local clustvar = "`clustervars'"
	ereturn scalar rmse = `rmse'
	ereturn scalar rss = `rss'
	ereturn scalar tss = `tss'
	if (`N_clust'<.) ereturn scalar N_clust = `N_clust'
	if (`N_clust'<.) ereturn scalar N_clust1 = `N_clust'
	ereturn `hidden' scalar unclustered_df_r = `CorrectDoF' // Used later in R2 adj

* Compute model F-test
	if (`K'>0) {
		qui test `indepvars' `avgevars' // Wald test
		ereturn scalar F = r(F)
		ereturn scalar df_m = r(df)
		ereturn scalar rank = r(df) // Not adding constant anymore
		if missing(e(F)) di as error "WARNING! Missing FStat"
	}
	else {
		ereturn scalar F = 0
		ereturn scalar df_m = 0
		ereturn scalar rank = 0 // Not adding constant anymore
	}
end
