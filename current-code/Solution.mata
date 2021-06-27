// --------------------------------------------------------------------------
// Regression Solution class
// --------------------------------------------------------------------------

mata:

class Solution
{
	// Created by partial_out()
	`RowVector'				tss
	`RowVector'				means
	`RowVector'				stdevs
	`RowVector'				norm2
	`Boolean'				is_standardized
	`Boolean'				has_intercept

	// Created by solver (LSMR, MAP, etc)
	`Integer'				stop_code // stop codes 1-9 are valid, >10 are errors
	`Matrix'				alphas
	`Matrix'				data

	// Core estimation results (saved by solve_ols)
	`Vector'				b
	`Matrix'				V
	`Integer'				N
	`Integer'				K
	`Integer'				rank // same as .K
	`Integer'				df_r
	`Integer'				df_a
	`Vector'				resid
	`RowVector'				kept
	
	// Estimation results - macros
	`String'				cmd
	`String'				cmdline
	`String'				model
	`String'				title
	
	// Estimation results - others
	`Boolean'				converged
	`Integer'				iteration_count // e(ic)
	`Varlist'				extended_absvars
	`Integer'				df_m
	`Integer'				N_clust
	`Integer'				N_clust_list
	`Real'					rss
	`Real'					rmse
	`Real'					F
	`Real'					tss_within
	`Real'					sumweights
	`Real'					r2
	`Real'					r2_within
	`Real'					r2_a
	`Real'					r2_a_within
	`Real'					ll
	`Real'					ll_0
	`Real'					accuracy

	// Copied from FixedEffects so we can add them with .post()
	// TODO: just save this to HDFE and run this from HDFE.post_footnote()
	`String'				weight_var          // Weighting variable
	`String'				weight_type         // Weight type (pw, fw, etc)
	`String'				vcetype
	`Integer'				num_clusters
	`Varlist'				clustervars

	// Parameters
	`Real'					collinear_tol		// Tolerance used to determine regressors collinear with fixed effects
	`Boolean'				fast_regression		// Set to 1 to choose a quick but less accurate regression
	`Boolean'				report_constant

	// Names of partialled-out variables	
	`Varname'				depvar
	`Varlist'				fullindepvars       // indepvars including basevars
	`Varlist'				fullindepvars_bn    // as 'fullindepvars', but all factor variables have a "bn" (used by sumhdfe)
	`Varlist'				indepvars           // x1 x2 x3
	`RowVector'				indepvar_status		// 1: basevar (e.g. includes i2000.year in ib2000.year); 2: collinear with FEs; 3: collinear with partialled-out regressors; 0: Ok! Recall that first cell is the depvar!
	`Varlist'				varlist             // y x1 x2 x3 ; set by HDFE.partial_out()
	`Varname'				residuals_varname

	// Methods
	`Void'					new()
	`Void'					check_collinear_with_fe()
	`Void'					expand_results()
	`Void'					post()
}
	
// Set default values
`Void' Solution::new()
{
	collinear_tol = . // must be set by user based on HDFE.tolerance
	converged = 0 // perhaps redundant with stop code? remove? (used by MAP while SC is used by LSMR)
	stop_code = 0 // solver hasn't run yet
	fast_regression = 0
	iteration_count = 0
	cmd = "reghdfe"
	model = "ols"
	title = "Linear regression"
	residuals_varname = ""
	report_constant = . // better to set this explicitly
}
	

// Flag regressors collinear with fixed effects and trim solution.data accordingly
`Void' Solution::check_collinear_with_fe(`Integer' verbose)
{
	`RowVector'				is_collinear, ok_index, collinear_index
	`Boolean'				first_is_depvar
	`Integer'				num_collinear, i

	// Note: standardizing makes it hard to detect values that are perfectly collinear with the absvars
	// in which case they should be 0.00 but they end up as e.g. 1e-16
	// EG: reghdfe price ibn.foreign , absorb(foreign)
	// We use 'collinear_tol' for that
	
	assert_msg(inrange(collinear_tol, 0, 1e-1), "error partialling out; missing values found")
	assert_msg(K == cols(data))

	// Abort if there are no regressors (recall that first is depvar)
	if (K<=1) return

	// Find out which variables are collinear with the FEs (i.e. are zero after partialling-out)
	is_collinear = (diagonal(cross(data, data))' :/ norm2) :<= (collinear_tol)

	// We must keep the depvar even if it's collinear
	first_is_depvar = (depvar != "") & (depvar == varlist[1])
	if (first_is_depvar & is_collinear[1]) {
		assert_msg(indepvar_status[1] == 0)
		if (verbose > -1) printf("{txt}warning: dependent variable {res}%s{txt} is likely perfectly explained by the fixed effects (tol =%3.1e)\n", depvar, collinear_tol)
		is_collinear[1] = 0

		// We'll also zero-out the partialled-out dependent variable
		// Otherwise, computing e(tss_within) would lead to very low but non-zer values, 
		// and thus log() of that will not be missing, and e(ll) and e(ll0) would be incorrect
		data[., 1] = J(rows(data), 1, 0)
	}

	// Number of collinear variables (after keeeping depvar)
	num_collinear = sum(is_collinear)

	// Trim solution.data
	if (num_collinear) {
		data = select(data, !is_collinear)
		tss = select(tss, !is_collinear)
		means = select(means, !is_collinear)
		stdevs = select(stdevs, !is_collinear)
		norm2 = select(norm2, !is_collinear)
	}

	// Update solution.indepvar_status of collinear variables (set to '2')
	ok_index = selectindex(indepvar_status:==0)
	collinear_index = select(ok_index, is_collinear)
	indepvar_status[collinear_index] = J(1, num_collinear, 2)

	// Warn about collinear regressors
	for (i=1; i<=K; i++) {
		if (is_collinear[i] & verbose>-1) {
			printf("{txt}note: {res}%s{txt} is probably collinear with the fixed effects (all partialled-out values are close to zero; tol =%3.1e)\n", varlist[i], collinear_tol)
		}
	}

	K = cols(data)
}


`Void' Solution::expand_results(`String' bname,
					  `String' Vname,
					  `Integer' verbose)
{
	`RowVector'				ok_index
	`Integer'				full_K
	`Vector'				full_b
	`Matrix'				full_V

	// Add constant
	if (report_constant) {
		if (verbose > 0) printf("{txt}# Adding _cons to varlist\n")
		indepvar_status = indepvar_status, 0
		fullindepvars = fullindepvars,  "_cons"
		fullindepvars_bn = fullindepvars_bn,  "_cons"
		indepvars = indepvars, "_cons"
	}

	// Expand 'b' and 'V' to include base and omitted variables (e.g. collinear vars)
	if (verbose > 1) printf("{txt}# Expanding 'b' and 'V' to include base and omitted variables\n")
	full_K = cols(indepvar_status) - 1 // recall that first cell is the depvar
	assert(full_K >= K)
	full_b = J(full_K, 1, 0)
	full_V = J(full_K, full_K, 0)

	if (K | report_constant) {
		ok_index = selectindex(indepvar_status[2..cols(indepvar_status)]:==0)
		full_b[ok_index] = b
		full_V[ok_index, ok_index] = V
	}
	
	st_matrix(bname, full_b')
	st_matrix(Vname, full_V)
}


`Void' Solution::post()
{
	`String'        text
	`Integer'       i

	// Scalars

	st_numscalar("e(N)", N)
	st_numscalar("e(rank)", rank)
	st_numscalar("e(df_r)", df_r)
	st_numscalar("e(df_m)", df_m)
	st_numscalar("e(tss)", tss)
	st_numscalar("e(tss_within)", tss_within)
	st_numscalar("e(rss)", rss)
	st_numscalar("e(mss)", tss - rss)
	st_numscalar("e(rmse)", rmse)
	st_numscalar("e(F)", F)
	st_numscalar("e(ll)", ll)
	st_numscalar("e(ll_0)", ll_0)
	st_numscalar("e(r2)", r2)
	st_numscalar("e(r2_within)", r2_within)
	st_numscalar("e(r2_a)", r2_a)
	st_numscalar("e(r2_a_within)", r2_a_within)
	if (!missing(sumweights)) st_numscalar("e(sumweights)", sumweights)
	st_numscalar("e(ic)", iteration_count)
	st_numscalar("e(converged)", converged)
	st_numscalar("e(report_constant)", report_constant)

	if (!missing(N_clust)) {
		st_numscalar("e(N_clust)", N_clust)
		st_numscalar("e(N_clustervars)", num_clusters)
		for (i=1; i<=num_clusters; i++) {
			text = sprintf("e(N_clust%g)", i)
			st_numscalar(text, N_clust_list[i])
			text = sprintf("e(clustvar%g)", i)
			st_global(text, clustervars[i])
		}
		text = "Statistics robust to heteroskedasticity"
		st_global("e(title3)", text)
		st_global("e(clustvar)", invtokens(clustervars))
	}

	// Macros

	st_global("e(cmd)", cmd)
	st_global("e(cmdline)", cmdline)
	st_global("e(model)", model)
	st_global("e(predict)", "reghdfe_p")
	st_global("e(estat_cmd)", "reghdfe_estat")
	st_global("e(footnote)", "reghdfe_footnote")
	//st_global("e(marginsok)", "")
	st_global("e(marginsnotok)", "Residuals SCore")
	st_global("e(depvar)", depvar)
	st_global("e(indepvars)", invtokens(indepvars))

	assert(title != "")
	text = sprintf("HDFE %s", title)
	st_global("e(title)", text)
		
	if (residuals_varname != "") {
		st_global("e(resid)", residuals_varname)
	}

	// Stata uses e(vcetype) for the SE column headers
	// In the default option, leave it empty.
	// In the cluster and robust options, set it as "Robust"
	text = strproper(vcetype)
	if (text=="Cluster") text = "Robust"
	if (text=="Unadjusted") text = ""
	assert(anyof( ("", "Robust", "Jackknife", "Bootstrap") , text))
	if (text!="") st_global("e(vcetype)", text)

	text = vcetype
	if (text=="unadjusted") text = "ols"
	st_global("e(vce)", text)

	// Weights
	if (weight_type != "") {
		st_global("e(wexp)", "= " + weight_var)
		st_global("e(wtype)", weight_type)
	}

}

end
