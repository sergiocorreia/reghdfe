// --------------------------------------------------------------------------
// Regression output
// --------------------------------------------------------------------------
// Important: avoid duplicating information with HDFE and HDFE.options

mata:

class FE_Output
{
	`String'		cmdline
	`String'		subcmd
	`String'		title
	`String'		model
	`Boolean'		converged
	`Integer'       iteration_count // e(ic)
	`Varlist'		absvars
	`Varlist'		extended_absvars
	`Integer'		N_hdfe
	`String'		notes
	`String'        equation_d

	// Fixed effect DoF computations
	`Integer'       G			        // Number of intercepts plus slopes
	`Integer'       G_extended          // Number of intercepts plus slopes
	`Integer'		dof_M				// e(mobility)
	`Integer'		dof_K
	`Integer'		df_a
	`Vector'		doflist_M
	`Vector'		doflist_K
	`Vector'		doflist_M_is_exact
	`Vector'		doflist_M_is_nested
	`Vector'		is_slope
	`Integer'		M_due_to_nested

	`Integer'		df_r
	`Integer'		df_m
	`Integer'		N_clust
	`Integer'		N_clust_list
	`Real'			rss
	`Real'			rmse
	`Real'			F
	`Real'			tss
	`Real'			tss_within
	`Real'			sumweights
	`Real'			r2
	`Real'			r2_within
	`Real'			r2_a
	`Real'			r2_a_within
	`Real'			ll
	`Real'			ll_0


	`Void'			post_footnote()
	`Void'			post()
}



`Void' FE_Output::post_footnote()
{
	`Matrix'				table
	`StringVector'  		rowstripe
	`StringRowVector' 		colstripe
	`String'				text

	text = invtokens(absvars)
	st_global("e(absvars)", text)
	
	text = invtokens(extended_absvars)
	text = subinstr(text, "1.", "")
	st_global("e(extended_absvars)", text)

	st_numscalar("e(df_a)", df_a)

	// Absorbed degrees-of-freedom table
	table = (doflist_K \ doflist_M \ (doflist_K-doflist_M) \ !doflist_M_is_exact \ doflist_M_is_nested)'
	rowstripe = extended_absvars'
	rowstripe = J(rows(table), 1, "") , extended_absvars' // add equation col
	colstripe = "Categories" \ "Redundant" \ "Num. Coefs." \ "Exact?" \ "Nested?"
	colstripe = J(cols(table), 1, "") , colstripe // add equation col
	st_matrix("e(dof_table)", table)
	st_matrixrowstripe("e(dof_table)", rowstripe)
	st_matrixcolstripe("e(dof_table)", colstripe)
}


`Void' FE_Output::post(`Options' options)
{
	`String'		text
	`Integer'		i

	post_footnote()

	// ---- constants -------------------------------------------------------

	st_global("e(predict)", "reghdfe_p")
	st_global("e(estat_cmd)", "reghdfe_estat")
	st_global("e(footnote)", "reghdfe_footnote")
	st_global("e(marginsok)", "")


	// ---- .output properties ----------------------------------------------

	assert(title != "")
	text = sprintf("HDFE %s", title)
	st_global("e(title)", text)
	
	text = sprintf("Absorbing %g HDFE %s", G, plural(G, "group"))
	st_global("e(title2)", text)
	
	st_global("e(model)", model)
	st_global("e(cmdline)", cmdline)

	st_numscalar("e(N_hdfe)", G)
	st_numscalar("e(N_hdfe_extended)", G_extended)
	st_numscalar("e(redundant)", dof_M)
	// st_numscalar("e(df_r)", df_r)
	st_numscalar("e(df_m)", df_m)
	//st_numscalar("e(rank)", df_m)
	st_numscalar("e(ic)", iteration_count)
	st_numscalar("e(rss)", rss)
	st_numscalar("e(rmse)", rmse)
	st_numscalar("e(tss)", tss)
	st_numscalar("e(tss_within)", tss_within)
	st_numscalar("e(mss)", tss - rss)
	st_numscalar("e(F)", F)
	st_numscalar("e(mobility)", dof_M)

	st_numscalar("e(ll)", ll)
	st_numscalar("e(ll_0)", ll_0)

	st_numscalar("e(M_due_to_nested)", M_due_to_nested)

	st_numscalar("e(r2)", r2)
	st_numscalar("e(r2_within)", r2_within)
	st_numscalar("e(r2_a)", r2_a)
	st_numscalar("e(r2_a_within)", r2_a_within)
	
	if (!missing(N_clust)) {
		st_numscalar("e(N_clust)", N_clust)
		for (i=1; i<=options.num_clusters; i++) {
			text = sprintf("e(N_clust%g)", i)
			st_numscalar(text, N_clust_list[i])
		}
		text = "Statistics robust to heteroskedasticity"
		st_global("e(title3)", text)
	}

	if (!missing(sumweights)) st_numscalar("e(sumweights)", sumweights)		


	// ---- .options properties ---------------------------------------------

	st_global("e(depvar)", options.depvar)
	st_global("e(indepvars)", options.indepvars)
	st_global("e(endogvars)", options.endogvars)
	st_global("e(instruments)", options.instruments)

	if (!missing(N_clust)) {
		st_numscalar("e(N_clustervars)", options.num_clusters)
		st_global("e(clustvar)", invtokens(options.clustervars))
		for (i=1; i<=options.num_clusters; i++) {
			text = sprintf("e(clustvar%g)", i)
			st_global(text, options.clustervars[i])
		}
	}

	if (options.residuals != "") {
		st_global("e(resid)", options.residuals)
	}

	// Stata uses e(vcetype) for the SE column headers
	// In the default option, leave it empty.
	// In the cluster and robust options, set it as "Robust"
	text = strproper(options.vcetype)
	if (text=="Cluster") text = "Robust"
	if (text=="Unadjusted") text = ""
	assert(anyof( ("", "Robust", "Jackknife", "Bootstrap") , text))
	if (text!="") st_global("e(vcetype)", text)

	text = options.vcetype
	if (text=="unadjusted") text = "ols"
	st_global("e(vce)", text)

	// Weights
	if (options.weight_type != "") {
		st_global("e(wexp)", "= " + options.weight_var)		
		st_global("e(wtype)", options.weight_type)
	}



}


end
