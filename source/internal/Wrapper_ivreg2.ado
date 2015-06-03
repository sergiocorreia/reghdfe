cap pr drop Wrapper_ivreg2
program define Wrapper_ivreg2, eclass
	syntax , depvar(varname) endogvars(varlist) instruments(varlist) ///
		[indepvars(varlist)] ///
		vceoption(string asis) ///
		KK(integer) ///
		[SHOWRAW(integer 0)] first(integer) [weightexp(string)] ///
		[ESTimator(string)] ///
		[num_clusters(string) clustervars(string)] ///
		[SUBOPTions(string)] [*] // [*] are ignored!
	if ("`options'"!="") Debug, level(3) msg("(ignored options: `options')")
	if (`c(version)'>=12) local hidden hidden
	
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
	
	if (`first') {
		local firstoption "first savefirst"
	}

	if ("`estimator'"!="2sls") local opt_estimator `estimator'
	
	* Variables have already been demeaned, so we need to add -nocons- or the matrix of orthog conditions will be singular
	if ("`cue'"=="") {
		local nocons nocons // Exception to get the same results as ivreg2, partial
	}
	else {
		local nocons nocons // partial(cons)
	}

	local subcmd ivreg2 `vars' `weightexp', `vceoption' `firstoption' small sdofminus(`kk') `nocons' `opt_estimator' `suboptions'
	Debug, level(3) msg(_n "call to subcommand: " _n as result "`subcmd'")
	local noise = cond(`showraw', "noi", "qui")
	if ("`noise'"=="noi") di as input _n "{title:Raw Results (without fixing tempnames):}"
	`noise' `subcmd'
	if ("`noise'"=="noi") di in red "{hline 64}" _n "{hline 64}"
	ereturn scalar tss = e(mss) + e(rss) // ivreg2 doesn't report e(tss)
	ereturn scalar unclustered_df_r = e(N) - e(df_m)

	if ("`e(vce)'"=="robust cluster") ereturn local vce = "cluster"

	if !missing(e(ecollin)) {
		di as error "endogenous covariate <`e(ecollin)'> was perfectly predicted by the instruments!"
		error 2000
	}

	if (`first') {
		ereturn `hidden' local first_prefix = "_ivreg2_"
	}

	local cats depvar instd insts inexog exexog collin dups ecollin clist redlist ///
		exexog1 inexog1 instd1 
	foreach cat in `cats' {
		FixVarnames `e(`cat')'
		ereturn local `cat' = "`r(newnames)'"
	}

	if (`first') {
		local firsteqs "`e(firsteqs)'"
		tempname hold
		estimates store `hold' , nocopy
		foreach fs_eqn in `firsteqs' {
			qui estimates restore `fs_eqn'

			foreach cat in `cats' {
				FixVarnames `e(`cat')'
				ereturn local `cat' = "`r(newnames)'"
			}

			* Fix e(clustvar) and e(clustvar#); modified from Post.ado
			if ("`e(clustvar)'"!="") {
				local subtitle = "`e(hacsubtitleV)'"
				if (`num_clusters'>1) {
					local rest `clustervars'
					forval i = 1/`num_clusters' {
						gettoken token rest : rest
						if (strpos("`e(clustvar`i')'", "__")==1) {
							local subtitle = subinstr("`subtitle'", "`e(clustvar`i')'", "`token'", 1)
						}
						ereturn local clustvar`i' `token'
					}
				}
				else {
					local subtitle = subinstr("`subtitle'", "`e(clustvar)'", "`clustervars'", 1)
				}
				ereturn scalar N_clustervars = `num_clusters'
				ereturn local clustvar `clustervars'
				ereturn local hacsubtitleV = "`subtitle'"
			}

			tempname b
			matrix `b' = e(b)
			local backup_colnames : colnames `b'
			FixVarnames `backup_colnames'
			matrix colnames `b' = `r(newnames)'
			ereturn repost b=`b', rename

			estimates store `fs_eqn', nocopy
		}
		qui estimates restore `hold'
		estimates drop `hold'

	}
end
