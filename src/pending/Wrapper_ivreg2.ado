cap pr drop Wrapper_ivreg2
pr Wrapper_ivreg2, eclass
	syntax , depvar(varname) endogvars(varlist) instruments(varlist) ///
		[indepvars(varlist)] ///
		vceoption(string asis) ///
		KK(integer) ///
		ffirst(integer) ///
		[weightexp(string)] ///
		[ESTimator(string)] ///
		[num_clusters(string) clustervars(string)] ///
		[SUBOPTions(string)] [*] // [*] are ignored!
	if ("`options'"!="") Debug, level(3) msg("(ignored options: `options')")
	if (c(version)>=12) local hidden hidden
	
	* Disable some options
	local 0 , `suboptions'
	syntax , [SAVEFPrefix(name)] [*] // Will ignore SAVEFPREFIX
	local suboptions `options'

	* Convert -vceoption- to what -ivreg2- expects
	local 0 `vceoption'
	syntax namelist(max=3) , [bw(string) dkraay(string) kernel(string) kiefer]
	gettoken vcetype transformed_clustervars : namelist
	local transformed_clustervars `transformed_clustervars' // Trim
	Assert inlist("`vcetype'", "unadjusted", "robust", "cluster")
	local vceoption = cond("`vcetype'"=="unadjusted", "", "`vcetype'")
	if ("`transformed_clustervars'"!="") local vceoption `vceoption'(`transformed_clustervars')
	if ("`bw'"!="") local vceoption `vceoption' bw(`bw')
	if ("`dkraay'"!="") local vceoption `vceoption' dkraay(`dkraay')
	if ("`kernel'"!="") local vceoption `vceoption' kernel(`kernel')
	if ("`kiefer'"!="") local vceoption `vceoption' kiefer
	
	mata: st_local("vars", strtrim(stritrim( "`depvar' `indepvars' (`endogvars'=`instruments')" )) )

	local opt small nocons sdofminus(`kk') `vceoption'  `suboptions'
	if (`ffirst') local opt `opt' ffirst
	if ("`estimator'"!="2sls") local opt `opt' `estimator'
	
	local subcmd ivreg2 `vars' `weightexp', `opt'
	Debug, level(3) msg(_n "call to subcommand: " _n as result "`subcmd'")
	qui `subcmd'
	ereturn scalar tss = e(mss) + e(rss) // ivreg2 doesn't report e(tss)
	ereturn scalar unclustered_df_r = e(N) - e(df_m)

	if ("`e(vce)'"=="robust cluster") ereturn local vce = "cluster"

	if !missing(e(ecollin)) {
		di as error "endogenous covariate <`e(ecollin)'> was perfectly predicted by the instruments!"
		error 2000
	}

	local cats depvar instd insts inexog exexog collin dups ecollin clist redlist ///
		exexog1 inexog1 instd1 
	foreach cat in `cats' {
		FixVarnames `e(`cat')'
		ereturn local `cat' = "`r(newnames)'"
	}
end
