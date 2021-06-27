// --------------------------------------------------------------------------
// Regression functions (least squares, compute VCE, etc.)
// --------------------------------------------------------------------------

mata:


`Void' function reghdfe_solve_ols(`FixedEffects' S, `Solution' sol, `String' vce_mode)
{
	`Matrix'				xx, inv_xx, W, inv_V, X
	`Vector'				xy, w, y
	`Integer'				tmp_N, used_df_r
	`Real'					stdev_y
	`RowVector'				stdev_x, is_collinear, idx

	// Hack: the first col of "sol.data" is actually y!

	if (S.verbose > 0) printf("{txt}# Solving least-squares regression of partialled-out variables\n\n")
	if (S.vcetype == "unadjusted" & S.weight_type=="pweight") S.vcetype = "robust"

	// vce_none: 		computes betas but not the VCE or extra info
	// vce_asymptotic:	DoF = N / (N-1)			note: where do I use this?
	// vce_small:		DoF = N / (N-K)
	assert_in(vce_mode, ("vce_none", "vce_small", "vce_asymptotic"))

	// Weight FAQ:
	// - fweight: obs. i represents w[i] duplicate obs. (there is no loss of info wrt to having the "full" dataset)
	// - aweight: obs. i represents w[i] distinct obs. that were mean-collapsed (so there is loss of info and hetero)
	//	 soln: normalize them so they sum to N (the true number of obs in our sample), and then treat them as fweight
	// - pweight: each obs. represents only one obs. from the pop, that was drawn from w[i] individuals
	//	          we want to make inference on the population, so if we interviewed 100% of the men and only 10% of women,
	//            then without weighting we would be over-representing men, which leads to a loss of efficiency +-+-
	// 			  it is the same as aweight + robust

	assert_msg(rows(sol.means) == 1, "nonconforming row(means)")
	assert_msg(cols(sol.means) == cols(sol.data), "nonconforming cols(means)")

	sol.K = cols(sol.data) - 1 // recall that first col is depvar 'y'
	sol.N = rows(sol.data) // This will change when using fweights

	// Regression weights
	if (rows(S.true_weights)) {
		// Custom case for IRLS (ppmlhdfe) where S.weights = mu * S.true_weights
		assert_msg(S.weight_type == "aweight")
		sol.N = sum(S.true_weights)
		w = S.weights * (sum(S.true_weights) / sum(S.weights))
	}
	else if (S.weight_type=="fweight") {
		sol.N = quadsum(S.weights)
		w = S.weights
	}
	else if (S.weight_type=="aweight" | S.weight_type=="pweight") {
		w = S.weights * (sol.N / quadsum(S.weights))
	}
	else {
		w = 1
	}

	assert_msg(!missing(S.weights_scale), "weights_scale must not be missing; did you ran .load_weights()?")
	sol.sumweights = S.weight_type != "" ? quadsum(S.weights) * S.weights_scale : sol.N

	// Build core matrices
	sol.K = cols(sol.data) - 1
	xx = quadcross(sol.data, w, sol.data) // at this point this is yx'yx not x'x
	sol.tss = sol.tss[1] // we initially save the TSS of all variables, not just depvar (at no cost, but helps with cache)
	sol.tss_within = xx[1,1]
	xy = sol.K ? xx[2..sol.K+1, 1] : J(0, 1, .)
	xx = sol.K ? xx[2..sol.K+1, 2..sol.K+1] : J(0, 0, .) // select all but the first row/column

	assert_msg( cols(sol.norm2) == sol.K+1		, "partial_out() was run with a different set of vars")
	assert_msg( cols(xx) == sol.K				, "x'x has the wrong number of columns")

	// Bread of the robust VCV matrix; compute this early so we can then update the list of collinear regressors
	inv_xx = reghdfe_rmcoll(xx, sol, is_collinear=., S.verbose)

	// Trim collinear variables (sol.norm2, sol.stdev, sol.means, xx, inv_xx, sol.K)
	// Also updates sol.indepvar_status, deletes sol.data, and creates 'y' and a trimmed-down 'X'
	reghdfe_trim_collinear(is_collinear, sol, xx, inv_xx, xy, y=., X=.)

	sol.df_m = sol.rank = sol.K
	sol.df_r = sol.N - S.df_a - sol.df_m // replaced when clustering

	// Compute betas
	if (!sol.K) {
		sol.b = J(0, 1, .)
	}
	else if (sol.fast_regression) {
		if (S.verbose > 0) printf("{txt}note: solving least-squares regression with a faster but less accurate method\n")
		// It's much faster to run qrsolve on xx and xy, but numerically inaccurate
		// See: http://www.stata.com/statalist/archive/2012-02/msg00956.html
		sol.b = qrsolve(xx, xy)
	}
	else if (S.has_weights) {
		sol.b = qrsolve(X :* sqrt(S.weights), y :* sqrt(S.weights))
	}
	else {
		sol.b = qrsolve(X, y)
	}

	// Compute residuals (simpler to do before adding back constant)
	if (vce_mode != "vce_none") sol.resid = y - X * sol.b

	// Add constant
	assert_boolean(sol.report_constant)
	if (sol.report_constant) {
		tmp_N = (S.weight_type=="aweight" | S.weight_type=="pweight") ? sol.N : sol.sumweights
		if (rows(S.true_weights)) tmp_N = sol.N
		reghdfe_extend_b_and_xx(sol, xx, inv_xx, tmp_N)
	}

	// Stop if no extra info needed (VCE, R2, RSS)
	if (vce_mode == "vce_none") {
		assert(!sol.is_standardized)
		return
	}

	sol.rss = quadcross(sol.resid, w, sol.resid) // do before reghdfe_vce_robust() modifies w

	// DO I NEED THIS???
	if (sol.report_constant) {
		X = sol.K ? X :+ sol.means[2..cols(sol.means)] : J(rows(X), 0, .)
	}

	// Compute VCE (update sol.V)
	assert_msg(anyof( ("unadjusted", "robust", "cluster") , S.vcetype), "invalid vcetype: " + S.vcetype)
	if (S.vcetype == "unadjusted") {
		reghdfe_vce_unadjusted(S, sol, inv_xx, vce_mode)
	}
	else if (S.vcetype == "robust") {
		reghdfe_vce_robust(S, sol, inv_xx, X, w, vce_mode)
	}
	else {
		reghdfe_vce_cluster(S, sol, inv_xx, X, w, vce_mode)
	}

	// Wald test: joint significance
	W = . // default when sol.K == 0
	if (sol.K) {
		idx = 1..sol.K
		inv_V = invsym(sol.V[idx, idx]) // this might not be of full rank but numerical inaccuracies hide it
		if (diag0cnt(inv_V)) {
			if (S.verbose > -1) printf("{txt}warning: missing F statistic; dropped variables due to collinearity or too few clusters\n")
			W = .
		}
		else {
			// We could probably do this with the simpler formula instead of Wald
			W = sol.b[idx]' * inv_V * sol.b[idx] / sol.df_m
			if (missing(W) & S.verbose > -1) printf("{txt}warning: missing F statistic\n")
		}
	}

	// V can be missing if all regressors are completely absorbed by the FEs
	if (missing(sol.V)) {
		if (S.verbose > 0) printf("{txt}   - VCE has missing values, setting it to zeroes (are your regressors all collinear?)\n")
		sol.V = J(rows(sol.V), rows(sol.V), 0)
	}

	// Undo standardization
	if (sol.is_standardized) {
		// Sanity checks
		assert(rows(sol.stdevs)==1)
		assert(sol.K == rows(sol.b) - sol.report_constant)
		assert(cols(sol.stdevs) - 1 == rows(sol.b) - sol.report_constant) // Subtract "y" on left; subtract "_cons" on right
		
		// Recover stdevs
		stdev_y = sol.stdevs[1]
		stdev_x = sol.K ? sol.stdevs[2..cols(sol.stdevs)] : J(1, 0, .)
		if (sol.report_constant) stdev_x = stdev_x, 1
		stdev_x = stdev_x :/ stdev_y
		
		// Transform output (note that S.tss is already ok)
		sol.rss = sol.rss * stdev_y ^ 2
		sol.tss_within = sol.tss_within * stdev_y ^ 2
		sol.resid = sol.resid * stdev_y
		sol.V = sol.V :/ (stdev_x' * stdev_x)
		sol.b = sol.b :/ stdev_x'
	}

	// Results
	used_df_r = sol.N - S.df_a - sol.df_m - S.df_a_nested
	sol.r2 = 1 - sol.rss / sol.tss
	sol.r2_a = 1 - (sol.rss / used_df_r) / (sol.tss / (sol.N - sol.has_intercept ) )
	sol.r2_within = 1 - sol.rss / sol.tss_within
	sol.r2_a_within = 1 - (sol.rss / used_df_r) / (sol.tss_within / (used_df_r + sol.rank))

	sol.ll = - 0.5 * sol.N * (1 + ln(2 * pi()) + ln(sol.rss / sol.N))
	sol.ll_0 = - 0.5 * sol.N * (1 + ln(2 * pi()) + ln(sol.tss_within / sol.N))

	sol.rmse = used_df_r ? sqrt(sol.rss / used_df_r) : sqrt(sol.rss)
	sol.F = W

	sol.vcetype = S.vcetype
	sol.weight_type = S.weight_type
	sol.weight_var = S.weight_var
	if (S.verbose > 0) printf("\n")
}


// --------------------------------------------------------------------------
// Remove collinear variables (based on ivreg2's s_rmcoll2)
// --------------------------------------------------------------------------
// This complements solution.check_collinear_with_fe()
`Matrix' reghdfe_rmcoll(`Matrix' xx, `Solution' sol, `RowVector' is_collinear, `Integer' verbose)
{
	`Integer'				num_dropped, i
	`Matrix'				inv_xx, alt_inv_xx
	`RowVector'				ok_index

	if (!sol.K) {
		is_collinear = J(1, 0, .)
		return(J(0, 0, .))
	}

	inv_xx = invsym(xx, 1..sol.K)
	num_dropped = diag0cnt(inv_xx)
	
	// Specifying the sweep order in invsym() can lead to incorrectly dropped regressors
	// (EG: with very VERY high weights)
	// We'll double check in this case and, if needed, fall back to the case without a preset sweep order
	if (num_dropped) {
		assert_msg(sol.K > 0 & sol.K == cols(inv_xx))
		alt_inv_xx = invsym(xx)
		if (num_dropped != diag0cnt(alt_inv_xx)) {
			inv_xx = alt_inv_xx
			num_dropped = diag0cnt(alt_inv_xx)
		}
	}

	// Flag collinear regressors
	is_collinear = !diagonal(inv_xx)'
	assert(cols(is_collinear) == sol.K)

	// Warn about collinear regressors
	ok_index = selectindex(sol.indepvar_status:==0)
	for (i=1; i<=sol.K; i++) {
		if (is_collinear[i] & verbose>-1) {
			printf("{txt}note: {res}%s{txt} omitted because of collinearity\n", sol.fullindepvars[ok_index[i]])
		}
	}
	
	return(inv_xx)
}


`Void' reghdfe_trim_collinear(`RowVector' is_collinear, `Solution' sol,
	                           `Matrix' xx, `Matrix' inv_xx, `Vector' xy,
	                           `Vector' y, `Matrix' X)
{
	`RowVector'				ok_index, collinear_index, extended_is_collinear
	`Integer'				num_collinear

	// Build 'y' and 'X' and discard sol.data
	y = sol.data[., 1]
	swap(X=., sol.data)
	X = select(X, (0, !is_collinear)) // Exclude 'y' and collinear regressors; will briefly use twice the memory taken by 'X'
	sol.K = cols(X)

	// Update xx, inv_xx, xy, and K
	ok_index = selectindex(!is_collinear)
	xx = xx[ok_index, ok_index]
	inv_xx = inv_xx[ok_index, ok_index]
	xy = xy[ok_index]

	assert_msg(cols(xx) == sol.K)
	assert_msg(cols(inv_xx) == sol.K)

	// Update other elements of sol
	extended_is_collinear = 0, is_collinear // includes 'y'
	sol.norm2 = select(sol.norm2, !extended_is_collinear)
	sol.stdevs = select(sol.stdevs, !extended_is_collinear)
	sol.means = select(sol.means, !extended_is_collinear)

	// Update solution.indepvar_status of collinear variables (set to '3')
	num_collinear = sum(extended_is_collinear)
	ok_index = selectindex(sol.indepvar_status:==0)
	collinear_index = select(ok_index, extended_is_collinear)
	if (num_collinear) sol.indepvar_status[collinear_index] = J(1, num_collinear, 3)
	assert_msg(sol.indepvar_status[1] == 0)  // ensure depvar was not touched
}


// --------------------------------------------------------------------------
// Use regression-through-mean and block partition formula to enlarge b and inv(XX)
// --------------------------------------------------------------------------
`Void' reghdfe_extend_b_and_xx(`Solution' sol, `Matrix' xx, `Matrix' inv_xx, `Integer' N)
{
	// How to add back _cons:
	// 1) To recover coefficient, apply "regression through means formula":
	//	  b0 = mean(y) - mean(x) * b1

	// 2) To recover variance ("full_inv_xx")
	//	  apply formula for inverse of partitioned symmetric matrix
	//    http://fourier.eng.hmc.edu/e161/lectures/gaussianprocess/node6.html
	//    http://www.cs.nthu.edu.tw/~jang/book/addenda/matinv/matinv/
	//
	//    Given A = [X'X X'1]		B = [B11 B21']		B = inv(A)
	//              [1'X 1'1]			[B21 B22 ]
	//
	//	  B11 is just inv(xx) (because of Frisch-Waugh)
	//	  B21 ("side") = means * B11
	//	  B22 ("corner") = 1 / sumweights * (1 - side * means')
	//
	//	- Note that means is NOT A12, but A12/N or A12 / (sum_weights)

	//  - Note: aw and pw (and unweighted) use normal weights,
	//	  but for fweights we expected sol.sumweights

	`RowVector'		means_x, side
	`Real'			b0, corner

	// Recalls that sol.means includes mean(y) in the first element
	means_x = sol.K ? sol.means[2..cols(sol.means)] : J(1, 0, .)

	// Update 'b'
	b0 = sol.means[1] - means_x * sol.b // mean(y) -  mean(x) * b1
	sol.b = sol.b \ b0

	// Update x'x
	corner = N
	side = N * means_x
	xx = (xx , side' \ side , corner)

	// Update inv(x'x)
	corner = (1 / N) + means_x * inv_xx * means_x'
	side = - means_x * inv_xx
	inv_xx = (inv_xx , side' \ side , corner)
}


`Void' reghdfe_vce_unadjusted(`FixedEffects' S, `Solution' sol, `Matrix' inv_xx, `String' vce_mode)
{
	`Integer'				dof_adj
	if (S.verbose > 0) {
		printf("{txt}   - Small-sample-adjustment: q = N / (N-df_m-df_a) = %g / (%g - %g - %g) = %g\n", sol.N, sol.N, sol.rank, S.df_a, sol.N / sol.df_r )
	}
	dof_adj = sol.N / sol.df_r
	if (vce_mode == "vce_asymptotic") dof_adj = sol.N / (sol.N-1) // 1.0
	sol.V = (sol.rss / sol.N) * dof_adj * inv_xx
}


// --------------------------------------------------------------------------
// Robust VCE
// --------------------------------------------------------------------------
// Advice: Delegate complicated regressions to -avar- and specialized routines
// BUGBUG: do we standardize X again? so V is well behaved?
// Notes:
// - robust is the same as cluster robust where cluster==_n
// - cluster just "collapses" X_i * e_i for each group, and builds M from that

`Void' reghdfe_vce_robust(`FixedEffects' S,
						   `Solution' sol,
						   `Matrix' D,
						   `Matrix' X,
						   `Variable' w,
						   `String' vce_mode)
{
	`Matrix'				M
	`Integer'				dof_adj

	if (S.verbose > 0) printf("{txt}# Estimating Robust Variance-Covariance Matrix of the Estimators (VCE)\n\n")
	if (S.verbose > 0) printf("{txt}   - VCE type: {res}%s{txt}\n", S.vcetype)
	if (S.verbose > 0) printf("{txt}   - Weight type: {res}%s{txt}\n", S.weight_type=="" ? "<none>" : S.weight_type)

	if (rows(S.true_weights)) {
		assert(S.weight_type=="aweight")
		w = (sol.resid :* w) :^ 2 :/ S.true_weights // resid^2 * aw^2 * fw
	}
	else if (S.weight_type=="") {
		w = sol.resid :^ 2
	}
	else if (S.weight_type=="fweight") {
		w = sol.resid :^ 2 :* w
	}
	else if (S.weight_type=="aweight" | S.weight_type=="pweight") {
		w = (sol.resid :* w) :^ 2
	}

	dof_adj = sol.N / (sol.N - S.df_a - sol.df_m)
	if (vce_mode == "vce_asymptotic") dof_adj = sol.N / (sol.N - 1) // 1.0
	if (S.verbose > 0 & vce_mode != "vce_asymptotic") {
		printf("{txt}   - Small-sample-adjustment: q = N / (N-df_m-df_a) = %g / (%g - %g - %g) = %g\n", sol.N, sol.N, sol.df_m, S.df_a, dof_adj )
	}
	M = sol.report_constant ? quadcross(X, 1, w, X, 1) : quadcross(X, w, X)
	sol.V = D * M * D * dof_adj
	w = J(0, 0, .) // ensure it's not used afterwards
}


`Void' reghdfe_vce_cluster(`FixedEffects' S,
						    `Solution' sol,
						    `Matrix' D,
						    `Matrix' X,
						    `Variable' w,
						    `String' vce_mode)
{
	`Matrix' 				M
	`Integer'				dof_adj, N_clust, df_r, nested_adj
	`Integer'				Q, q, g, sign, i, j
	`FactorPointers'		FPlist
	`FactorPointer'			FP
	`Varlist'				vars
	`String'				var, var_with_spaces
	`Boolean'				clustervar_is_absvar, required_fix
	`Matrix'				tuples
	`RowVector'				tuple
	`RowVector'				N_clust_list
	`Matrix'				joined_levels
	`Integer'				Msize

	w = sol.resid :* w
	Msize = sol.K + sol.report_constant

	vars = S.clustervars
	Q = cols(vars)
	if (S.verbose > 0) {
		printf("{txt}# Estimating Cluster Robust Variance-Covariance Matrix of the Estimators (VCE)\n\n")
		printf("{txt}   - VCE type: {res}%s{txt} (%g-way clustering)\n", S.vcetype, Q)
		printf("{txt}   - Cluster variables: {res}%s{txt}\n", invtokens(vars))
		printf("{txt}   - Weight type: {res}%s{txt}\n", S.weight_type=="" ? "<none>" : S.weight_type)
	}

	assert_msg(1 <= Q & Q <= 10)

	// Get or build factors associated with the clustervars
	FPlist = J(1, Q, NULL)
	N_clust_list = J(1, Q, .)
	for (q=1; q<=Q; q++) {
		var = vars[q]
		clustervar_is_absvar = 0
		for (g=1; g<=S.G; g++) {
			if (invtokens(S.factors[g].varlist, "#") == var) {
				clustervar_is_absvar = 1
				FP = &(S.factors[g])
				break
			}
		}
		var_with_spaces = subinstr(var, "#", " ")
		if (!clustervar_is_absvar) FP = &(factor(var_with_spaces, S.sample, ., "", ., ., ., 0))
		N_clust_list[q] = (*FP).num_levels
		if (S.verbose > 0) printf("{txt}   - {res}%s{txt} has {res}%g{txt} levels\n", var, N_clust_list[q])
		FPlist[q] = FP
	}

	// Build the meat part of the V matrix
	if (S.verbose > 0) printf("{txt}   - Computing the 'meat' of the VCE\n")
	M = J(Msize, Msize, 0)
	tuples = .
	for (q=1; q<=Q; q++) {
		tuples = reghdfe_choose_n_k(Q, q, tuples)
		sign = mod(q, 2) ? 1 : -1 // "+" with odd number of variables, "-" with even
		for (j=1; j<=rows(tuples); j++) {
			tuple = tuples[j, .]
			if (S.verbose > 0) printf("{txt}      - Level %g/%g; sublevel %g/%g; M = M %s ClusterVCE(%s)\n", q, Q, j, rows(tuples), sign > 0 ? "+" : "-" , invtokens(strofreal(tuple)))
			if (q==1) {
				assert(tuple==j)
				FP =  FPlist[j]
			}
			else if (q==2) {
				FP = &join_factors( *FPlist[tuple[1]] , *FPlist[tuple[2]] , ., ., 1)
			}
			else {
				joined_levels = (*FPlist[tuple[1]]).levels
				for (i=2; i<=cols(tuple); i++) {
					joined_levels = joined_levels, (*FPlist[tuple[i]]).levels
				}
				FP = &_factor(joined_levels, ., ., "", ., ., ., 0)
			}
			M = M + sign * reghdfe_vce_cluster_meat(FP, X, w, Msize, sol.report_constant)
		}
	}

	// Build VCE
	N_clust = min(N_clust_list)
	
	nested_adj = S.df_a_nested > 0 // minor adj. so we match xtreg when the absvar is nested within cluster
	// (when ..nested.., df_a is zero so we divide N-1 by something that can potentially be N (!))
	// so we either add the 1 back, or change the numerator (and the N_clust-1 factor!)
	// addendum: this just ensures we subtract the constant when we have nested FEs

	dof_adj = (sol.N - 1) / (sol.N - nested_adj - sol.df_m-S.df_a) * N_clust / (N_clust - 1) // adjust for more than 1 cluster
	if (vce_mode == "vce_asymptotic") dof_adj = N_clust / (N_clust - 1)  // 1.0
	if (S.verbose > 0) {
		printf("{txt}   - Small-sample-adjustment: q = (%g - 1) / (%g - %g) * %g / (%g - 1) = %g\n", sol.N, sol.N, sol.df_m+S.df_a+nested_adj, N_clust, N_clust, dof_adj)
	}
	sol.V = D * M * D * dof_adj
	if (Q > 1) {
		required_fix = reghdfe_fix_psd(sol.V, S.solution.report_constant)
		if (required_fix) printf("{txt}Warning: VCV matrix was non-positive semi-definite; adjustment from Cameron, Gelbach & Miller applied.\n")
	}

	// Store e()
	assert(!missing(sol.df_r))
	df_r = N_clust - 1
	if (sol.df_r > df_r) {
		sol.df_r = df_r
	}
	else if (S.verbose > 0) {
		printf("{txt}   - Unclustered df_r (N - df_m - df_a = %g) are {it:lower} than clustered df_r (N_clust-1 = %g)\n", sol.df_r, df_r)
		printf("{txt}     Thus, we set e(df_r) as the former.\n")
		printf("{txt}     This breaks consistency with areg but ensures internal consistency\n")
		printf("{txt}     between vce(robust) and vce(cluster _n)\n")
	}

	sol.N_clust = N_clust
	sol.N_clust_list = N_clust_list
	sol.num_clusters = Q
	sol.clustervars = S.clustervars
	if (S.verbose > 0) printf("\n")
}



`Matrix' reghdfe_vce_cluster_meat(`FactorPointer' FP,
                                  `Variables' X,
                                  `Variable' resid,
                                  `Integer' Msize,
                                  `Boolean' report_constant)
{
	`Integer'				i, N_clust
	`Variables'				X_sorted
	`Variable'				resid_sorted
	`Matrix'				X_tmp
	`Vector'				resid_tmp
	`RowVector'				Xe_tmp
	`Matrix'				M

	if (cols(X)==0 & !report_constant) return(J(0,0,0))

	N_clust = (*FP).num_levels
	(*FP).panelsetup()
	X_sorted = (*FP).sort(X)
	resid_sorted = (*FP).sort(resid)
	M = J(Msize, Msize, 0)

	if (cols(X)) {
		for (i=1; i<=N_clust; i++) {
			X_tmp = panelsubmatrix(X_sorted, i, (*FP).info)
			resid_tmp = panelsubmatrix(resid_sorted, i, (*FP).info)
			Xe_tmp = quadcross(1, 0, resid_tmp, X_tmp, report_constant) // Faster than colsum(e_tmp :* X_tmp)
			M = M + quadcross(Xe_tmp, Xe_tmp)
		}		
	}
	else {
		// Workaround for when there are no Xs except for _cons
		assert(report_constant)
		for (i=1; i<=N_clust; i++) {
			resid_tmp = panelsubmatrix(resid_sorted, i, (*FP).info)
			M = M + quadsum(resid_tmp) ^ 2
		}
	}

	return(M)
}


// Enumerate all combinations of K integers from N integers
// Kroneker approach based on njc's tuples.ado
`Matrix' reghdfe_choose_n_k(`Integer' n, `Integer' k, `Matrix' prev_ans)
{
	`RowVector' v
	`Integer' q
	`Matrix' candidate
	`Matrix' ans
	v = 1::n
	if (k==1) return(v)

	q = rows(prev_ans)
	assert(q==comb(n, k-1))
	assert(cols(prev_ans)==k-1)
	candidate = v # J(q, 1, 1)
	candidate = candidate , J(n, 1, prev_ans)	
	ans = select(candidate, candidate[., 1] :< candidate[., 2])
	return(ans)
}


// --------------------------------------------------------------------------
// Fix non-positive VCV
// --------------------------------------------------------------------------
// If the VCV matrix is not positive-semidefinite, use the fix from
// Cameron, Gelbach & Miller - Robust Inference with Multi-way Clustering (JBES 2011)
// 1) Use eigendecomposition V = U Lambda U' where U are the eigenvectors and Lambda = diag(eigenvalues)
// 2) Replace negative eigenvalues into zero and obtain FixedLambda
// 3) Recover FixedV = U * FixedLambda * U'
`Boolean' function reghdfe_fix_psd(`Matrix' V, `Boolean' report_constant) {
	`Matrix'				U
	`Vector'				lambda
	`Boolean' 				required_fix
	`Matrix'				V_backup
	`RowVector'				index

	if (!issymmetric(V)) _makesymmetric(V)
	if (!issymmetric(V)) exit(error(505))

	if (report_constant) {
		index = 1..cols(V)-1
		V_backup = V[index, index]
	}

	symeigensystem(V, U=., lambda=.)
	required_fix = min(lambda) < 0

	if (required_fix) {
		lambda = lambda :* (lambda :>= 0)
		V = quadcross(U', lambda, U') // V = U * diag(lambda) * U'
		
		// If V is positive-semidefinite, then submatrix V_backup also is
		if (report_constant) {
			(void) reghdfe_fix_psd(V_backup, 0)
			V[index, index] = V_backup
		}
	}

	return(required_fix)
}

end
