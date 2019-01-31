// Common functions ---------------------------------------------------------
mata:

// --------------------------------------------------------------------------
// BUGBUG: not sure if this is still used...
// --------------------------------------------------------------------------
`StringRowVector' clean_tokens(`String' vars)
{
	`StringRowVector' 		ans
	`Integer'				i
	ans = tokens(vars)
	for (i=1; i<=cols(ans); i++) {
		ans[i] = invtokens(tokens(ans[i]))
	}
	return(ans)
}


// --------------------------------------------------------------------------
// Workaround to st_data's odd behavior
// --------------------------------------------------------------------------
// This does three things:
// 1) Wrap up interactions in parens (up to four) to avoid Stata's quirk/bug
// 2) If issue persists, load columns one-by-one
// 1) Instead of returning it reuses existing matrices (might use less mem?)
//
// Example of the issue:
// 	sysuse auto, clear
//	mata: cols(st_data(., "1.rep78 2.rep78 3.rep78#1.foreign")) // expected 3, got 6
// Happens b/c st_data doesn't work variable by variable but expands the interactions
// We can partly fix it by surrounding interactions with parens
// But st_data() only supports up to 4 parens


`Void' _st_data_wrapper(`Variables' index, `StringRowVector' vars, `Variables' data, `Boolean' verbose)
{
	`RowVector'			is_interaction
	`StringRowVector'	fixed_vars
	`Integer'			i, k

	vars = tokens(invtokens(vars))

	// Add parenthesis only for Stata 11-14, as on Stata 15+ they are i) not needed and ii) corrupt output
	// For i) see "help set fvtrack"
	// For ii) see "test/stdata3.do" on Github
	if (st_numscalar("c(stata_version)") < 15) {
		is_interaction = strpos(vars, "#") :> 0
		is_interaction = is_interaction :& (runningsum(is_interaction) :<= 4) // Only up to 4 parenthesis supported
		fixed_vars = subinstr(strofreal(is_interaction), "0", "")
		fixed_vars = subinstr(fixed_vars, "1", "(") :+ vars :+ subinstr(fixed_vars, "1", ")")
	}
	else {
		fixed_vars = vars
	}

	// Override code above, to minimize any risk of incorrect data
	// Since this is an undocumented feature, it might or might not work on some older versions of Stata
	// (See also email from jpitblado@stata.com)
	fixed_vars = vars

	data = st_data(index, fixed_vars)
	k = cols(vars)

	if (cols(data) > k) {
	    if (verbose > 0) printf("{err}(some empty columns were added due to a bug/quirk in {bf:st_data()}; %g cols created instead of %g for {it:%s}; running slower workaround)\n", cols(data), k, invtokens(vars))
	    data = J(rows(data), 0, .)
	    for (i=1; i<=k; i++) {
	        data = data, st_data(index, vars[i])
	    }
	}
	assert(cols(data)==k)
}
	

// --------------------------------------------------------------------------
// Each col of A will have stdev of 1 unless stdev is quite close to 0
// --------------------------------------------------------------------------
`RowVector' function reghdfe_standardize(`Matrix' A)
{
	`RowVector'				stdevs, means
	`Integer'				K, N // i, 

	// We don't need to good accuracy for the stdevs, so we have a few alternatives:
	// Note: cross(1,A) is the same as colsum(A), but faster
	// Note: cross(A, A) is very fast, but we only need the main diagonals
	// [A: 1sec] stdevs = sqrt( (colsum(A:*A) - (cross(1, A) :^ 2 / N)) / (N-1) )
	// [B: .61s] stdevs = sqrt( (diagonal(cross(A, A))' - (cross(1, A) :^ 2 / N)) / (N-1) )
	// [C: .80s] stdevs = diagonal(sqrt(variance(A)))'
	// [D: .67s] means = cross(1, A) / N; stdevs =  sqrt(diagonal(crossdev(A, means, A, means))' / (N-1))

	assert_msg(!isfleeting(A), "input cannot be fleeting")
	N = rows(A)
	K = cols(A)

	stdevs = J(1, K, .)

	// (A) Very precise

	// (B) Precise
	// means = cross(1, A) / N
	// stdevs =  sqrt(diagonal(quadcrossdev(A, means, A, means))' / (N-1))

	// (C) 20% faster; don't use it if you care about accuracy
	stdevs = sqrt( (diagonal(cross(A, A))' - (cross(1, A) :^ 2 / N)) / (N-1) )
	assert_msg(!missing(stdevs), "stdevs are missing; is N==1?") // Shouldn't happen as we don't expect N==1
	stdevs = colmax(( stdevs \ J(1, K, 1e-3) ))
	A = A :/ stdevs

	// (D) Equilibrate matrix columns instead of standardize (i.e. just divide by column max)
	// _perhapsequilc(A, stdevs=.)
	// stdevs = 1 :/ stdevs
	// assert_msg(!missing(stdevs), "stdevs are missing; is N==1?")

	// (E) Don't do anything
	// stdevs = J(1, cols(A), 1)

	return(stdevs)
}


// --------------------------------------------------------------------------
// Divide two row vectors but adjust the denominator if it's too small
// --------------------------------------------------------------------------
`RowVector' safe_divide(`RowVector' numerator, `RowVector' denominator, | `Real' epsi) {
	 // If the numerator goes below machine precision, we lose accuracy
	 // If the denominator goes below machine precision, the division explodes
	 if (args()<3 | epsi==.) epsi = epsilon(1)
	return( numerator :/ colmax(denominator \ J(1,cols(denominator),epsi)) )
}


// If X is not square...
// `Matrix' R
// real colvector tau, p

// _hqrdp(A, tau, R, p=.)
// B = hqrdmultq1t(A, tau, B)
// rank = _solveupper(R, B, tol)
// B = B[invorder(p),.]
// +- +-

// invsym(makesymmetric(..))





// --------------------------------------------------------------------------
// Robust solver for Ax=b
// --------------------------------------------------------------------------
// Mata utility for sequential use of solvers
// Default is cholesky;
// if that fails, use QR;
// if overridden, use QR.
// Author: Schaffer, Mark E <M.E.Schaffer@hw.ac.uk>
// --------------------------------------------------------------------------
// Warning:
// cholqrsolve calls qrsolve which calls _qrsolve which calls ...
// Does all the indirection makes it too slow to use within a panel?
// --------------------------------------------------------------------------
`Matrix' function reghdfe_cholqrsolve(`Matrix' A,
                                      `Matrix' B,
                                    | `Boolean' useqr)
{
	`Matrix' C
	if (args()<3 | useqr==.) useqr = 0
	
	if (!useqr) {
		C = cholsolve(A, B)
		if (hasmissing(C)) useqr = 1
	}

	if (useqr) {
		C = qrsolve(A, B)
	}

	return(C)
}


// --------------------------------------------------------------------------
// OLS Regression
// --------------------------------------------------------------------------
`Void' function reghdfe_post_ols(`FixedEffects' S,
                                 `Variables' X,
                                 `String' bname,
                                 `String' Vname,
                                 `String' nname,
                                 `String' rname,
                                 `String' dfrname)
{
	`Integer'				N
	`Integer'				rank
	`Integer'				df_r
	`Vector'				b
	`Matrix'				V
	`Variable'				resid
	`Real'					eps
	`Integer'				i
	`RowVector'				kept
	`Vector'				not_basevar


	`Vector'				idx
	`Vector'				temp_b
	`Matrix'				temp_V
	`Integer'				k

	if (S.timeit) timer_on(90)
	reghdfe_solve_ols(S, X, b=., V=., N=., rank=., df_r=., resid=., kept=., "vce_small")
	assert(cols(X) - 1 == rows(b) - S.compute_constant) // The 1st column of X is actually Y
	assert((rows(b) == rows(V)) & (rows(b) == cols(V)))
	if (S.timeit) timer_off(90)

	// Add base vars
	if (S.compute_constant) {
		if (S.verbose > 1) printf("\n{txt}## Adding _cons to varlist\n")
		assert_msg(rows(S.not_basevar) == 1, "rows(S.not_basevar) == 1")
		S.not_basevar = S.not_basevar, 1
		S.fullindepvars = S.fullindepvars + " _cons"
		S.indepvars = S.indepvars + " _cons"
	}
	if (S.not_basevar != J(1, 0, .)) {
		if (S.verbose > 1) printf("\n{txt}## Adding base variables to varlist\n")
		k = cols(S.not_basevar)
		assert_msg(cols(S.not_basevar) == k, "cols(S.not_basevar) == k")
		idx = `selectindex'(S.not_basevar)
		swap(b, temp_b)
		swap(V, temp_V)
		b = J(k, 1, 0)
		V = J(k, k, 0)
		b[idx, 1] = temp_b
		V[idx, idx] = temp_V
	}

	st_matrix(bname, b')

	if (S.verbose > 1) printf("\n{txt}## Reporting omitted variables\n")
	// Add "o." prefix to omitted regressors
	eps = sqrt(epsilon(1))
	for (i=1; i<=rows(b); i++) {
		if (b[i]==0 & S.not_basevar[i] & S.verbose > -1) {
			printf("{txt}note: %s omitted because of collinearity\n", tokens(S.fullindepvars)[i])
			//stata(sprintf("_ms_put_omit %s", indepvars[i]))
			//indepvars[i] = st_global("s(ospec)")
			// This is now one in reghdfe.ado with -_ms_findomitted-
		}
	}

	st_matrix(Vname, V)
	st_numscalar(nname, N)
	st_numscalar(rname, rank)
	st_numscalar(dfrname, df_r)

	// Need to save resids if saving FEs, even if temporarily
	if (S.residuals == "" & S.save_any_fe) {
		S.residuals = "__temp_reghdfe_resid__"
	}

	if (S.residuals != "") {
		if (S.verbose > 0) printf("\n{txt}## Storing residuals in {res}%s{txt}\n\n", S.residuals)
		if (S.compact == 1) {
			S.residuals_vector = resid
		}
		else {
			S.save_variable(S.residuals, resid, "Residuals")
		}
	}
}


`Void' function reghdfe_solve_ols(`FixedEffects' S,
                                  `Variables' X,
                                  `Vector' b,
                                  `Matrix' V,
                                  `Integer' N,
                                  `Integer' rank,
                                  `Integer' df_r,
                                  `Vector' resid,
                                  `RowVector' kept,
                                  `String' vce_mode,
                                | `Variable' true_w)
{
	// Hack: the first col of X is actually y!
	`Integer'				K, KK, tmp_N
	`Matrix'				xx, inv_xx, W, inv_V, just_X
	`Vector' 				w
	`Integer'				used_df_r
	`Integer'				dof_adj
	
	`Boolean'				is_standardized
	`Real'					stdev_y
	`RowVector'				stdev_x

	if (true_w == . | args() < 11) true_w = J(0, 1, .)
	if (S.vcetype == "unadjusted" & S.weight_type=="pweight") S.vcetype = "robust"
	if (S.verbose > 0) printf("\n{txt}## Solving least-squares regression of partialled-out variables\n\n")
	assert_in(vce_mode, ("vce_none", "vce_small", "vce_asymptotic"))

	is_standardized = S.all_stdevs != J(1, 0, .)
	if (is_standardized) S.means = S.means :/ S.all_stdevs

	// Weight FAQ:
	// - fweight: obs. i represents w[i] duplicate obs. (there is no loss of info wrt to having the "full" dataset)
	// - aweight: obs. i represents w[i] distinct obs. that were mean-collapsed (so there is loss of info and hetero)
	//	 soln: normalize them so they sum to N (the true number of obs in our sample), and then treat them as fweight
	// - pweight: each obs. represents only one obs. from the pop, that was drawn from w[i] individuals
	//	          we want to make inference on the population, so if we interviewed 100% of the men and only 10% of women,
	//            then without weighting we would be over-representing men, which leads to a loss of efficiency +-+-
	// 			  it is the same as aweight + robust
	// We need to pick N and w
	N = rows(X) // Default; will change with fweights
	S.sumweights = S.weight_type != "" ? quadsum(S.weight) : N
	assert(rows(S.means) == 1)
	assert(cols(S.means) == cols(X))

	w = 1
	if (rows(true_w)) {
		// Custom case for IRLS (ppmlhdfe) where S.weight = mu * true_w
		assert_msg(S.weight_type == "aweight")
		N = sum(true_w)
		w = S.weight * sum(true_w) / sum(S.weight)
	}
	else if (S.weight_type=="fweight") {
		N = S.sumweights
		w = S.weight
	}
	else if (S.weight_type=="aweight" | S.weight_type=="pweight") {
		w = S.weight * (N / S.sumweights)
	}

	// Build core matrices
	if (S.timeit) timer_on(91)

	K = cols(X) - 1
	xx = quadcross(X, w, X)
	S.tss_within = xx[1,1]
	xx = K ? xx[| 2 , 2 \ K+1 , K+1 |] : J(0, 0, .)
	if (S.timeit) timer_off(91)

	// This matrix indicates what regressors are not collinear
	assert_msg(cols(S.kept)==K+1, "partial_out() was run with a different set of vars")

	// Bread of the robust VCV matrix
	// Compute this early so we can update the list of collinear regressors
	if (S.timeit) timer_on(95)
	assert_msg( cols(tokens(invtokens(S.indepvars)))==cols(xx) , "HDFE.indepvars is missing or has the wrong number of columns")
	inv_xx = reghdfe_rmcoll(tokens(invtokens(S.indepvars)), xx, kept) // this modifies -kept-

	// // Workaround for case with extremely high weights, where ivnsym loses precision and incorrectly excludes vars
	// if (S.has_weights) {
	// 	if (max(S.weight) > 1e5) {
	// 		kept = (1..K)
	// 	}
	// }

	S.df_m = rank = K - diag0cnt(inv_xx)
	KK = S.df_a + S.df_m
	S.df_r = N - KK // replaced when clustering
	if (S.timeit) timer_off(95)

	// Compute betas
	// - There are two main options
	//	 a) Use cholqrsolve on xx and xy. Faster but numerically inaccurate
	//      See: http://www.stata.com/statalist/archive/2012-02/msg00956.html
	//   b) Use qrsolve. More accurate but doesn't handle weights easily
	// - Ended up doing (b) with a hack for weights
	b = J(K, 1, 0)
	if (cols(kept)) {
		if (S.has_weights) {
			b[kept] = qrsolve(X[., 1:+kept] :* sqrt(S.weight), X[., 1] :* sqrt(S.weight))
		}
		else {
			b[kept] = qrsolve(X[., 1:+kept], X[., 1])
		}
	}

	if (S.timeit) timer_on(92)
	if (!isfleeting(resid) | vce_mode != "vce_none") resid = X * (1 \ -b) // y - X * b
	if (S.timeit) timer_off(92)

	if (S.compute_constant) {
		tmp_N = (S.weight_type=="aweight" | S.weight_type=="pweight") ? N : S.sumweights
		if (rows(true_w)) tmp_N = N
		reghdfe_extend_b_and_inv_xx(S.means, tmp_N, b, inv_xx)
	}

	// Stop if no VCE/R2/RSS needed
	if (vce_mode == "vce_none") {
		assert(!is_standardized)
		return
	}

	if (S.timeit) timer_on(93)
	if (S.vcetype != "unadjusted") {
		if (S.compute_constant) {
			if (isfleeting(X)) {
				// Save some memory... unsure if it helps
				swap(just_X, X)
				just_X = K ? just_X[., 2..K+1] :+ S.means[2..cols(S.means)] : J(rows(just_X), 0, .)
			}
			else {
				just_X = K ? X[., 2..K+1] :+ S.means[2..cols(S.means)] : J(rows(X), 0, .)
			}
		}
		else {
			just_X = K ? X[., 2..K+1] : J(rows(X), 0, .)
		}
	}
	if (S.timeit) timer_off(93)

	if (S.timeit) timer_on(94)
	S.rss = quadcross(resid, w, resid) // do before reghdfe_robust() modifies w
	if (S.timeit) timer_off(94)

	// Compute full VCE
	if (S.timeit) timer_on(96)
	assert_msg(anyof( ("unadjusted", "robust", "cluster") , S.vcetype), "invalid vcetype" + S.vcetype)
	if (S.vcetype == "unadjusted") {
		if (S.verbose > 0) {
			printf("{txt}   - Small-sample-adjustment: q = N / (N-df_m-df_a) = %g / (%g - %g - %g) = %g\n", N, N, rank, S.df_a, N / S.df_r )
		}
		dof_adj = N / S.df_r
		if (vce_mode == "vce_asymptotic") dof_adj = N / (N-1) // 1.0
		V = (S.rss / N) * dof_adj * inv_xx
	}
	else if (S.vcetype == "robust") {
		V = reghdfe_robust(S, just_X, inv_xx, resid, w, N, KK, vce_mode, true_w)
	}
	else {
		V = reghdfe_cluster(S, just_X, inv_xx, resid, w, N, KK, vce_mode)
	}
	if (S.timeit) timer_off(96)

	// Wald test: joint significance
	if (S.timeit) timer_on(97)
	inv_V = invsym(V[kept, kept]) // this might not be of full rank but numerical inaccuracies hide it
	if (diag0cnt(inv_V)) {
		if (S.verbose > -1) printf("{txt}warning: missing F statistic; dropped variables due to collinearity or too few clusters\n")
		W = .
	}
	else if (length(b[kept])==0) {
		W = .
	}
	else {
		// We could probably do this with the simpler formula instead of Wald
		W = b[kept]' * inv_V * b[kept] / S.df_m
		if (missing(W) & S.verbose > -1) printf("{txt}warning: missing F statistic\n")
	}
	if (S.timeit) timer_off(97)

	// V can be missing if b is completely absorbed by the FEs
	if (missing(V)) {
		if (S.verbose > 0) printf("{txt}   - VCE has missing values, setting it to zeroes (are your regressors all collinear?)\n")
		V = J(rows(V), rows(V), 0)
	}

	// Undo standardization
	if (is_standardized) {
		// Sanity checks
		assert(rows(S.all_stdevs)==1)
		assert(cols(S.all_stdevs) - 1 == rows(b) - S.compute_constant) // Subtract "y" on left; subtract "_cons" on right
		
		// Recover stdevs
		stdev_y = S.all_stdevs[1]
		stdev_x = K ? S.all_stdevs[2..cols(S.all_stdevs)] : J(1, 0, .)
		if (S.compute_constant) stdev_x = stdev_x, 1
		stdev_x = stdev_x :/ stdev_y
		
		// Transform output (note that S.tss is already ok)
		S.rss = S.rss * stdev_y ^ 2
		S.tss_within = S.tss_within * stdev_y ^ 2
		resid = resid * stdev_y
		V = V :/ (stdev_x' * stdev_x)
		b = b :/ stdev_x'
	}

	// Results
	S.title = "Linear regression"
	// S.model = "ols"
	used_df_r = N - KK - S.df_a_nested
	S.r2 = 1 - S.rss / S.tss
	S.r2_a = 1 - (S.rss / used_df_r) / (S.tss / (N - S.has_intercept ) )
	S.r2_within = 1 - S.rss / S.tss_within
	S.r2_a_within = 1 - (S.rss / used_df_r) / (S.tss_within / (used_df_r + rank))

	S.ll = - 0.5 * N * (1 + ln(2 * pi()) + ln(S.rss / N))
	S.ll_0 = - 0.5 * N * (1 + ln(2 * pi()) + ln(S.tss_within / N))

	S.rmse = sqrt(S.rss / used_df_r)
	if (used_df_r==0) S.rmse = sqrt(S.rss)
	S.F = W
	df_r = S.df_r // reghdfe_cluster might have updated it (this gets returned to the caller function)
}


// --------------------------------------------------------------------------
// Robust VCE
// --------------------------------------------------------------------------
// Advice: Delegate complicated regressions to -avar- and specialized routines
// BUGBUG: do we standardize X again? so V is well behaved?
// Notes:
// - robust is the same as cluster robust where cluster==_n
// - cluster just "collapses" X_i * e_i for each group, and builds M from that

`Matrix' reghdfe_robust(`FixedEffects' S,
                        `Variables' X,
						`Matrix' D,
						`Variable' resid,
						`Variable' w,
						`Integer' N,
						`Integer' K,
						`String' vce_mode,
					    `Variable' true_w)
{
	`Matrix'				M, V
	`Integer'				dof_adj

	if (S.verbose > 0) printf("\n{txt}## Estimating Robust Variance-Covariance Matrix of the Estimators (VCE)\n\n")
	if (S.verbose > 0) printf("{txt}   - VCE type: {res}%s{txt}\n", S.vcetype)
	if (S.verbose > 0) printf("{txt}   - Weight type: {res}%s{txt}\n", S.weight_type=="" ? "<none>" : S.weight_type)

	if (rows(true_w)) {
		assert(S.weight_type=="aweight")
		w = (resid :* w) :^ 2 :/ true_w // resid^2 * aw^2 * fw
	}
	else if (S.weight_type=="") {
		w = resid :^ 2
	}
	else if (S.weight_type=="fweight") {
		w = resid :^ 2 :* w
	}
	else if (S.weight_type=="aweight" | S.weight_type=="pweight") {
		w = (resid :* w) :^ 2
	}

	dof_adj = N / (N - K)
	if (vce_mode == "vce_asymptotic") dof_adj = N / (N-1) // 1.0
	M = S.compute_constant ? quadcross(X, 1, w, X, 1) : quadcross(X, w, X)
	if (S.verbose > 0) {
		printf("{txt}   - Small-sample-adjustment: q = N / (N-df_m-df_a) = %g / (%g - %g - %g) = %g\n", N, N, K-S.df_a, S.df_a, N / (N-K) )
	}
	V = D * M * D * dof_adj
	return(V)
}

`Matrix' reghdfe_cluster(`FixedEffects' S,
                        `Variables' X,
						`Matrix' D,
						`Variable' resid,
						`Variable' w,
						`Integer' N,
						`Integer' K,
						`String' vce_mode)
{
	`Matrix' 				M, V
	`Integer'				dof_adj, N_clust, df_r, nested_adj
	`Integer'				Q, q, g, sign, i, j
	pointer(`Factor') rowvector FPlist
	`FactorPointer'			FP
	`Varlist'				vars
	`String'				var, var_with_spaces
	`Boolean'				clustervar_is_absvar, required_fix
	`Matrix'				tuples
	`RowVector'				tuple
	`RowVector'				N_clust_list
	`Matrix'				joined_levels
	`Integer'				Msize

	w = resid :* w
	Msize = cols(X) + S.compute_constant

	vars = S.clustervars
	Q = cols(vars)
	if (S.verbose > 0) printf("\n{txt}## Estimating Cluster Robust Variance-Covariance Matrix of the Estimators (VCE)\n\n")
	if (S.verbose > 0) printf("{txt}   - VCE type: {res}%s{txt} (%g-way clustering)\n", S.vcetype, Q)
	if (S.verbose > 0) printf("{txt}   - Cluster variables: {res}%s{txt}\n", invtokens(vars))
	if (S.verbose > 0) printf("{txt}   - Weight type: {res}%s{txt}\n", S.weight_type=="" ? "<none>" : S.weight_type)
	assert_msg(0 < Q & Q < 10)

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
		sign = mod(q, 2) ? 1 : -1 // + with odd number of variables, - with even
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
			M = M + sign * reghdfe_vce_cluster_meat(FP, X, w, Msize, S.compute_constant)
		}
	}

	// Build VCE
	N_clust = min(N_clust_list)
	nested_adj = (S.df_a==0) // minor adj. so we match xtreg when the absvar is nested within cluster
	// (when ..nested.., df_a is zero so we divide N-1 by something that can potentially be N (!))
	// so we either add the 1 back, or change the numerator (and the N_clust-1 factor!)
	dof_adj = (N - 1) / (N - nested_adj - K) * N_clust / (N_clust - 1) // adjust for more than 1 cluster
	if (vce_mode == "vce_asymptotic") dof_adj = N_clust / (N_clust - 1)  // 1.0
	if (S.verbose > 0) {
		printf("{txt}   - Small-sample-adjustment: q = (%g - 1) / (%g - %g) * %g / (%g - 1) = %g\n", N, N, K+nested_adj, N_clust, N_clust, dof_adj)
	}
	V = D * M * D * dof_adj
	if (Q > 1) {
		required_fix = reghdfe_fix_psd(V)
		if (required_fix) printf("{txt}Warning: VCV matrix was non-positive semi-definite; adjustment from Cameron, Gelbach & Miller applied.\n")
	}

	// Store e()
	assert(!missing(S.df_r))
	df_r = N_clust - 1
	if (S.df_r > df_r) {
		S.df_r = df_r
	}
	else if (S.verbose > 0) {
		printf("{txt}   - Unclustered df_r (N - df_m - df_a = %g) are {it:lower} than clustered df_r (N_clust-1 = %g)\n", S.df_r, df_r)
		printf("{txt}     Thus, we set e(df_r) as the former.\n")
		printf("{txt}     This breaks consistency with areg but ensures internal consistency\n")
		printf("{txt}     between vce(robust) and vce(cluster _n)\n")
	}
	
	S.N_clust = N_clust
	S.N_clust_list = N_clust_list

	return(V)
}


`Matrix' reghdfe_vce_cluster_meat(`FactorPointer' FP,
                                  `Variables' X,
                                  `Variable' resid,
                                  `Integer' Msize,
                                  `Boolean' compute_constant)
{
	`Integer'				i, N_clust
	`Variables'				X_sorted
	`Variable'				resid_sorted
	`Matrix'				X_tmp
	`Vector'				resid_tmp
	`RowVector'				Xe_tmp
	`Matrix'				M

	if (cols(X)==0 & !compute_constant) return(J(0,0,0))

	N_clust = (*FP).num_levels
	(*FP).panelsetup()
	X_sorted = (*FP).sort(X)
	resid_sorted = (*FP).sort(resid)
	M = J(Msize, Msize, 0)

	if (cols(X)) {
		for (i=1; i<=N_clust; i++) {
			X_tmp = panelsubmatrix(X_sorted, i, (*FP).info)
			resid_tmp = panelsubmatrix(resid_sorted, i, (*FP).info)
			Xe_tmp = quadcross(1, 0, resid_tmp, X_tmp, compute_constant) // Faster than colsum(e_tmp :* X_tmp)
			M = M + quadcross(Xe_tmp, Xe_tmp)
		}		
	}
	else {
		// Workaround for when there are no Xs except for _cons
		assert(compute_constant)
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
`Boolean' function reghdfe_fix_psd(`Matrix' V) {
	`Matrix'				U
	`Matrix'				lambda
	`Boolean' 				required_fix

	if (!issymmetric(V)) _makesymmetric(V)
	if (!issymmetric(V)) exit(error(505))
	symeigensystem(V, U=., lambda=.)
	if (min(lambda)<0) {
		lambda = lambda :* (lambda :>= 0)
		// V = U * diag(lambda) * U'
		V = quadcross(U', lambda, U')
		required_fix = 1
	}
	else {
		required_fix = 0
	}
	return(required_fix)
}


// --------------------------------------------------------------------------
// Remove collinear variables
// --------------------------------------------------------------------------
// Based on ivreg2's s_rmcoll2
`Matrix' reghdfe_rmcoll(`Varlist' varnames,
                        `Matrix' xx,
                        `RowVector' kept)
{
	`Integer'				K, num_dropped
	`Matrix'				inv_xx, smat, alt_inv_xx
	`RowVector'				vl_drop, vl_keep

	assert(rows(xx)==cols(xx))
	K = cols(xx)
	inv_xx = K ? invsym(xx, 1..K) : J(0, 0, .)
	
	// Specifying the sweep order in invsym() can lead to incorrectly dropped regressors
	// (EG: with very VERY high weights)
	// We'll double check in this case
	num_dropped = diag0cnt(inv_xx)
	if (K & num_dropped) {
		alt_inv_xx = invsym(xx)
		if (num_dropped != diag0cnt(alt_inv_xx)) {
			inv_xx = alt_inv_xx
			num_dropped = diag0cnt(alt_inv_xx)
		}
	}

	st_numscalar("r(k_omitted)", num_dropped)
	smat = (diagonal(inv_xx) :== 0)'
	vl_drop = select(varnames, smat)
	vl_keep = select(varnames, !smat)
	if (cols(vl_keep)) st_global("r(varlist)", invtokens(vl_keep))
	if (cols(vl_drop)) st_global("r(omitted)", invtokens(vl_drop))
	kept = `selectindex'(!smat) // Return it, so we can exclude these variables from the joint Wald test
	return(inv_xx)
}


// --------------------------------------------------------------------------
// Use regression-through-mean and block partition formula to enlarge b and inv(XX)
// --------------------------------------------------------------------------
`Void' reghdfe_extend_b_and_inv_xx(
	`RowVector' 	means,
	`Integer'		N,
	`Vector'		b,
	`Matrix' 		inv_xx)
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
	//	  but for fweights we expected S.sumweights

	`RowVector'		means_x, side
	`Real'			corner

	means_x = cols(means) > 1 ? means[2..cols(means)] : J(1, 0, .)
	b = b \ means[1] - means_x * b // means * (1 \ -b)
	corner = (1 / N) + means_x * inv_xx * means_x'
	side = - means_x * inv_xx
	inv_xx = (inv_xx , side' \ side , corner)

}

end
