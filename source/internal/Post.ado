capture program drop Post
program define Post, eclass
	syntax, coefnames(string) ///
		model(string) stage(string) stages(string) subcmd(string) cmdline(string) vceoption(string) original_absvars(string) extended_absvars(string) vcetype(string) vcesuite(string) tss(string) num_clusters(string) ///
		[dofadjustments(string) clustervars(string) timevar(string) r2c(string) equation_d(string) subpredict(string) savefirst(string) diopts(string) weightvar(string) gmm2s(string) cue(string) dkraay(string) liml(string) by(string) level(string)] ///
		[backup_original_depvar(string) original_indepvars(string) original_endogvars(string) original_instruments(string)]

	if (`c(version)'>=12) local hidden hidden // ereturn hidden requires v12+

	Assert e(tss)<., msg("within tss is missing")
	Assert `tss'<., msg("overall tss is missing")
	Assert e(N)<., msg("# obs. missing in e()")

	* Why is this here and not right after FixVarnames?
	* Because of some Stata black magic, if I repost *before* the restore this will not work
	ereturn repost b=`coefnames', rename

	if ("`weightvar'"!="") {
		qui su `weightvar', mean
		ereturn scalar sumweights = r(sum)
	}

* Absorbed-specific returns
	* e(N_hdfe) e(N_hdfe_extended) e(mobility)==M e(df_a)==K-M
	* e(M#) e(K#) e(M#_exact) e(M#_nested) -> for #=1/e(N_hdfe_extended)
	mata: map_ereturn_dof(HDFE_S)
	local N_hdfe = e(N_hdfe)
	Assert e(df_r)<. , msg("e(df_r) is missing")

* MAIN LOCALS
	ereturn local cmd = "reghdfe"
	ereturn local cmdline `"`cmdline'"'
	ereturn local subcmd = cond(inlist("`stage'", "none", "iv"), "`subcmd'", "regress")
	ereturn local model = cond("`gmm2s'"=="", "`model'", "gmm2s")
	ereturn local model = cond("`cue'"=="", "`model'", "cue")
	ereturn local model = cond("`liml'"=="", "`model'", "liml")
	ereturn local dofadjustments = "`dofadjustments'"
	ereturn local title = "HDFE " + e(title)
	ereturn local title2 =  "Absorbing `N_hdfe' HDFE " + plural(`N_hdfe', "group")
	ereturn local predict = "reghdfe_p"
	ereturn local estat_cmd = "reghdfe_estat"
	ereturn local footnote = "reghdfe_footnote"
	ereturn `hidden' local equation_d = "`equation_d'" // The equation used to construct -d- (used to predict)
	ereturn local absvars = "`original_absvars'"
	ereturn `hidden' local extended_absvars = "`extended_absvars'"
	

	ereturn `hidden' local diopts = "`diopts'"
	ereturn `hidden' local subpredict = "`subpredict'"

* CLUSTER AND VCE
	
	ereturn local vcesuite = "`vcesuite'"
	if ("`e(subcmd)'"=="ivreg2") local vcesuite = "avar" // This is what ivreg2 uses
	if ("`e(subcmd)'"=="ivregress") local vcesuite = "default"

	* Replace __CL#__ and __ID#__ from cluster subtitles

	if ("`e(clustvar)'"!="") {
		if ("`e(subcmd)'"=="ivreg2") local subtitle = "`e(hacsubtitleV)'"
		if (`num_clusters'>1) {
			local rest `clustervars'
			forval i = 1/`num_clusters' {
				gettoken token rest : rest
				if ("`e(subcmd)'"=="ivreg2" & strpos("`e(clustvar`i')'", "__")==1) {
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
		if ("`e(subcmd)'"=="ivreg2") ereturn local hacsubtitleV = "`subtitle'"
	}
	if (`dkraay'>1) {
		ereturn local clustvar `timevar'
		ereturn scalar N_clustervars = 1
	}

	
	* Stata uses e(vcetype) for the SE column headers
	* In the default option, leave it empty.
	* In the cluster and robust options, set it as "Robust"
	ereturn local vcetype = proper("`vcetype'") //
	if (e(vcetype)=="Cluster") ereturn local vcetype = "Robust"
	if (e(vcetype)=="Unadjusted") ereturn local vcetype
	if ("`e(vce)'"=="." | "`e(vce)'"=="") ereturn local vce = "`vcetype'" // +-+-
	Assert inlist("`e(vcetype)'", "", "Robust", "Jackknife", "Bootstrap")

* STAGE
	if ("`stage'"!="none") ereturn local iv_depvar = "`backup_original_depvar'"

* VARLISTS
	* Besides each cmd's naming style (e.g. exogr, exexog, etc.) keep one common one
	foreach cat in indepvars endogvars instruments {
		if ("`original_`cat''"=="") continue
		ereturn local `cat' "`original_`cat''"
	}

* MAIN NUMERICS
	ereturn `hidden' scalar tss_within = e(tss)
	ereturn scalar tss = `tss'
	ereturn scalar mss = e(tss) - e(rss)
	ereturn scalar ll   = -0.5 * (e(N)*ln(2*_pi) + e(N)*ln(e(rss)       /e(N)) + e(N))
	ereturn scalar ll_0 = -0.5 * (e(N)*ln(2*_pi) + e(N)*ln(e(tss_within)/e(N)) + e(N))
	ereturn scalar r2 = 1 - e(rss) / e(tss)
	ereturn scalar r2_within = 1 - e(rss) / e(tss_within)

	* ivreg2 uses e(r2c) and e(r2u) for centered/uncetered R2; overwrite first and discard second
	if (e(r2c)!=.) {
		ereturn scalar r2c = e(r2)
		ereturn scalar r2u = .
	}

	* Computing Adj R2 with clustered SEs is tricky because it doesn't use the adjusted inputs:
	* 1) It uses N instead of N_clust
	* 2) For the DoFs, it uses N - Parameters instead of N_clust-1
	* 3) Further, to compute the parameters, it includes those nested within clusters
	
	* Note that this adjustment is NOT PERFECT because we won't compute the mobility groups just for improving the r2a
	* (when a FE is nested within a cluster, we don't need to compute mobilty groups; but to get the same R2a as other estimators we may want to do it)
	* Instead, you can set by hand the dof() argument and remove -cluster- from the list

	if ("`model'"=="ols" & `num_clusters'>0) Assert e(unclustered_df_r)<., msg("wtf-`vcesuite'")
	local used_df_r = cond(e(unclustered_df_r)<., e(unclustered_df_r), e(df_r)) - e(M_due_to_nested)
	ereturn scalar r2_a = 1 - (e(rss)/`used_df_r') / ( e(tss) / (e(N)-1) )
	ereturn scalar rmse = sqrt( e(rss) / `used_df_r' )
	ereturn scalar r2_a_within = 1 - (e(rss)/`used_df_r') / ( e(tss_within) / (`used_df_r'+e(df_m)) )

	if (e(N_clust)<.) Assert e(df_r) == e(N_clust) - 1, msg("Error, `wrapper' should have made sure that N_clust-1==df_r")
	*if (e(N_clust)<.) ereturn scalar df_r = e(N_clust) - 1

	if ("`model'"=="ols" & inlist("`vcetype'", "unadjusted", "ols")) {
		 // -1 b/c we exclude constant for this
		 ereturn scalar F_absorb = (e(r2)-`r2c') / (1-e(r2)) * e(df_r) / (e(df_a)-1)

		//if (`nested') {
		//	local rss`N_hdfe' = e(rss)
		//	local temp_dof = e(N) - e(df_m) // What if there are absorbed collinear with the other RHS vars?
		//	local j 0
		//	ereturn `hidden' scalar rss0 = `rss0'
		//	forv g=1/`N_hdfe' {
		//		local temp_dof = `temp_dof' - e(K`g') + e(M`g')
		//		*di in red "g=`g' RSS=`rss`g'' and was `rss`j''.  dof=`temp_dof'"
		//		ereturn `hidden' scalar rss`g' = `rss`g''
		//		ereturn `hidden' scalar df_a`g' = e(K`g') - e(M`g')
		//		local df_a_g = e(df_a`g') - (`g'==1)
		//		ereturn scalar F_absorb`g' = (`rss`j''-`rss`g'') / `rss`g'' * `temp_dof' / `df_a_g'
		//		ereturn `hidden' scalar df_r`g' = `temp_dof'
		//		local j `g'
		//	}   
		//}
	}

	if ("`savefirst'"!="") ereturn `hidden' scalar savefirst = `savefirst'

	* We have to replace -unadjusted- or else subsequent calls to -suest- will fail
	Subtitle `vceoption' // will set title2, etc. Run after e(bw) and all the others are set!
	if (e(vce)=="unadjusted") ereturn local vce = "ols"

	if ("`stages'"!="none") {
		ereturn local stage = "`stage'"
		ereturn `hidden' local stages = "`stages'"
	}
end
