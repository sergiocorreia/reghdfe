// --------------------------------------------------------------------------
// Code that extends Factor() for standard FEs
// --------------------------------------------------------------------------

mata:

// Notes on preconditioning:
//- What we call "diagonal" is a type of right-preconditioning described here: https://web.stanford.edu/group/SOL/software/lsmr/
//- Instead of solving Ax=b, we solve a similar system where the columns of A have norm=2:
//	Call D = inv(diag(A'A)); then  Ax=b --> A D inv(D) x = b --> (AD) z = b.
//	After finding the solution "z", then we solve Dz=x and solve for x (which we might not even need to do if we just care about the resid Ax-b)
//	Note also that we can save the preconditioner as just a vector with size cols(A)
//	ADx --> first multiply Dx and then apply A
//	D'A'x --> compute A'x and then multipy by D'

// Main class ---------------------------------------------------------------

class FE_Factor extends Factor
{
	// Main properties
	`String' 				absvar				// expression: "firm#year", "i.firm##c.(x1 x2)", etc.
	`Varlist' 				ivars				// variables used for intercepts
	`Varlist' 				cvars				// variables used for slopes
	`Boolean'				has_intercept		// 1 for "firm" or "firm##c.year" ; 0 for "firm#c.year"
	`Integer' 				num_slopes			// number of slopes
	`String' 				target				// where the FE will be saved
	`Boolean'				save_fe				// 1 if we save the fixed effects (the "alphas")
	`Integer'				num_singletons		// Useful if we want to see which FE is causing more singletons (sumhdfe)

	// Properties used for weights
	`Boolean'				has_weights
	`Vector'				weights
	`Vector'				weighted_counts
	
	// Properties used for slope variables (cvars)
	`Matrix'				unsorted_x			// standardized slope variables "x1 x2.."
	`RowVector'				x_means				// means of slope variables
	`RowVector'				x_stdevs			// standard deviations of slope variables
	`Matrix'				inv_xx				// inv(x'x) for each fixed effect
	`Matrix'				xx_info				// Fake .info matrix used with block_diagonal preconditioner
	
	// Preconditioner info
	`String'				preconditioner // We need to remember the preconditioner used in order to undo it (when recovering the alphas)
	`Matrix'				preconditioner_intercept, preconditioner_slopes

	// Saved coefs
	`Vector'				tmp_alphas
	`Vector'				alphas

	`Boolean'				is_individual_fe
	`String'				individual_fe_function

	// Prevent compile errors (we must place them here and not in Individual_Factor class)
	//`Factor'				FG, FI
	//`BipartiteGraph'		bg
	`Vector'				group_p, group_inv_p

	// Methods
	`Void'					new()
	virtual `Vector'		drop_singletons()	// Adjust to dropping obs.
	virtual `Void'			init_diagonal_preconditioner()
	virtual `Void'			init_block_diag_preconditioner()
	virtual `Void'			mult()  			// b = Ax  , implemented as A.mult(x, b)
	virtual `Vector'		mult_transpose()  	// x = A'b , implemented as x = A.mult_transpose(b)
	virtual `Void'			undo_preconditioning()
	virtual `Matrix'		panelsum()
	virtual `Matrix'		panelmean()
	virtual `Void'			set_weights()

	// These methods don't do anything (are just pass through) but we declare them so "mata desc" does not crash
	virtual `Void'			drop_obs()
	virtual `Void'			panelsetup()

	virtual `Void'			cleanup_before_saving()
}


`Void' FE_Factor::new()
{
	is_individual_fe = 0
	num_singletons = 0
	has_weights = 0
	weights = 1
	weighted_counts = .
}


`Vector' FE_Factor::drop_singletons(| `Vector' fweight,
                                      `Boolean' zero_threshold)
{
	`Integer'				num_singletons
	`Vector'				mask, idx
	`Boolean'				has_fweight
	`Vector'				w_counts // don't name this "weighted_counts" b/c that's the name of a class property

	// - By default, this drops all singletons (obs where F.counts==1)
	// - If fweights are provided, we'll only drop those singletons with fweight of 1
	// - As a hack, if zero_threshold==1, we'll drop singletons AND all obs where 
	//   "w_counts" (actually depvar) is zero
	//   Also, we multiply by counts so we can track how many actual obs were dropped

	if (zero_threshold == .) zero_threshold = 0

	if (counts == J(0, 1, .)) {
		_error(123, "drop_singletons() requires the -counts- vector")
	}

	has_fweight = (args()>=1) & (fweight!=. & rows(fweight)>0) // either missing value or 0x1 vector are valid cases for no fweights

	if (has_fweight) {
		assert(rows(fweight)==num_obs)
		this.panelsetup()
		w_counts = this.panelsum(fweight, 1, 1)
		if (zero_threshold) {
			mask = (!w_counts :| (counts :== 1)) :* counts
		}
		else {
			mask = w_counts :== 1
		}
	}
	else {
		mask = (counts :== 1)
	}

	num_singletons = sum(mask)
	if (num_singletons == 0) return(J(0, 1, .))
	counts = counts - mask
	idx = selectindex(mask[levels, .])

	// Update and overwrite fweight
	if (has_fweight) {
		fweight = num_singletons == num_obs ? J(0, 1, .) : select(fweight, (!mask)[levels])
	}
	
	// Update contents of F based on just idx and the updated F.counts
	__inner_drop(idx)
	return(idx)
}


`Void' FE_Factor::panelsetup()
{
	super.panelsetup()
	assert(this.panel_is_setup==1)
}


`Matrix' FE_Factor::panelsum(`Matrix' x, | `Boolean' sort_input, `Boolean' ignore_weights)
{
	if (args()<2 | sort_input == .) sort_input = 0
	if (args()<3 | ignore_weights == .) ignore_weights = 0

	assert_boolean(sort_input)
	assert_boolean(ignore_weights)
	assert_boolean(this.has_weights)
	assert_msg(this.panel_is_setup == 1, "must call .panelsetup() before .panelsum()")
	assert_msg(this.has_weights==1 | this.weights==1)

	if (ignore_weights) {
		return(::panelsum(sort_input ? this.sort(x) : x,               this.info))
	}
	else {
		return(::panelsum(sort_input ? this.sort(x) : x, this.weights, this.info))
	}
}


`Variables' FE_Factor::panelmean(`Variables' y, | `Boolean' sort_input)
{
	`Variables' ans

	if (args()<2 | sort_input == .) sort_input = 0
	assert_boolean(this.has_weights)
	assert_msg(this.panel_is_setup == 1, "this.panel_setup() must be run beforehand")
	assert_msg(this.has_weights==1 | this.weights==1)

	ans = ::panelsum(sort_input ? this.sort(y) : y, this.weights, this.info)

	if (this.has_weights) {
		return(editmissing(ans :/ this.weighted_counts, 0))
	}
	else {
		return(ans :/ this.counts)
	}
}


`Void' FE_Factor::set_weights(`Vector' input_weights, | `Boolean' verbose)
{
	if (args()<2 | verbose == .) verbose = 0
	if (verbose > 0) printf("{txt}   - sorting weights for factor {res}%s{txt}\n", this.absvar)

	assert_msg(this.has_weights == 0)
	assert_msg(this.weights == 1)
	assert_msg(this.weighted_counts == .)
	assert_msg(this.num_obs == rows(input_weights))

	this.has_weights = 1
	this.weights = this.sort(input_weights)
	this.weighted_counts = this.panelsum(this.weights, 0, 1) // sort_input=0 ignore_weights=1
}


`Void' FE_Factor::drop_obs(`Vector' idx)
{
	super.drop_obs(idx)
}


`Void' FE_Factor::init_diagonal_preconditioner()
{
	`Matrix'				precond

	if (this.has_intercept) {
		precond = this.has_weights ? this.weighted_counts : this.counts
		this.preconditioner_intercept = sqrt(1 :/ precond)
	}

	if (this.num_slopes) {
		precond = this.panelsum(this.unsorted_x :^ 2, 1)
		this.preconditioner_slopes = 1 :/ sqrt(precond)
	}
}


`Void' FE_Factor::init_block_diag_preconditioner(`String' weight_type, `Real' weights_scale, `Vector' true_weights)
{
	`Matrix'				precond
	`Matrix'				sorted_x
	`Matrix'                x_means, inv_xx
	`Matrix'                tmp_inv_xx, tmp_x, tmp_xx
	`RowVector'             tmp_x_means
	`Vector'				tmp_w
	`Integer'				i, offset, L, K, FullK, N

	`Integer'				B11
	`RowVector'				B12

	`Matrix'				V
	`Vector'				d
	`Vector'				sorted_true_weights

	L = this.num_levels
	K = this.num_slopes
	FullK = this.has_intercept + this.num_slopes

	// No point in using block-diagonal with only one regressor
	if (FullK == 1) {
		this.init_diagonal_preconditioner()
		return
	}

	// Compute inv(xx) using block inverse formula
	// See: http://fourier.eng.hmc.edu/e161/lectures/gaussianprocess/node6.html
	// See code in : reghdfe_extend_b_and_xx()
	sorted_x = this.sort(this.unsorted_x)
	if (this.has_intercept) this.x_means = this.panelmean(sorted_x, 0) // We only need x_means if we have an intercept

	inv_xx = J(L*FullK, FullK, .)
	if (rows(true_weights)) sorted_true_weights = this.sort(true_weights)

	for (i=1; i<=L; i++) {
	
		tmp_x = panelsubmatrix(sorted_x, i, this.info)
		tmp_w = has_weights ? panelsubmatrix(this.weights, i, this.info) : 1
		offset = FullK * (i - 1)

		if (rows(true_weights)) {
			// Custom case for IRLS (ppmlhdfe) where S.weights = mu * S.true_weights
			assert_msg(weight_type == "aweight")
			N = quadsum(panelsubmatrix(sorted_true_weights, i, this.info))
		}
		else if (weight_type=="fweight") {
			N = quadsum(panelsubmatrix(this.weights, i, this.info))
		}
		else if (weight_type=="" | weight_type=="aweight" | weight_type=="pweight") {
			N = rows(tmp_x)
		}
		else {
			assert_msg(0, "preconditioner missing one case - iweight???")
		}

		if (1) {
			if (this.has_intercept)  {
			    tmp_x_means = K > 1 ? this.x_means[i, .] : this.x_means[i]
			    tmp_inv_xx = invsym(quadcrossdev(tmp_x, tmp_x_means, tmp_w, tmp_x, tmp_x_means))
			    // We just computed B22 (see link above), need to compute entire B matrix
			    B11 = 1 / N + tmp_x_means * tmp_inv_xx * tmp_x_means'
			    B12 = -tmp_x_means * tmp_inv_xx
			    tmp_inv_xx = (B11, B12 \ B12', tmp_inv_xx)
			}
			else {
			    tmp_inv_xx = invsym(quadcross(tmp_x, tmp_w, tmp_x))
			}

			// Square root of PSD matrix? Use eigenvalue VDV'
			// https://math.stackexchange.com/questions/349721/square-root-of-positive-definite-matrix
			// Note that if A=VDV' where D is diagonal, then V sqrt(D) V' is the square root (trivial to prove because V is orthogonal i.e. V'=inv(V))
			eigensystem(tmp_inv_xx, V=., d=.)
			tmp_inv_xx = makesymmetric(Re(V * diag(sqrt(d)) * V'))

			inv_xx[|offset + 1, 1 \ offset + FullK , . |] = tmp_inv_xx
		}
		else {
			// Doesn't work, probably an error in the formula or due to collinear vars
			tmp_xx = quadcross(tmp_x, tmp_w, tmp_x)
			if (this.has_intercept)  {
			    tmp_x_means = K > 1 ? this.x_means[i, .] : this.x_means[i]
			    tmp_xx = N, N*tmp_x_means \ N*tmp_x_means' , tmp_xx
			}

			//"ALT3"
			//tmp_xx = J(rows(tmp_x), 1, 1) , tmp_x
			//tmp_xx = quadcross(tmp_xx, tmp_xx)
			//invsym(tmp_xx)
			//tmp_xx

			eigensystem(tmp_xx, V=., d=.)
			inv_xx[|offset + 1, 1 \ offset + FullK , . |] = makesymmetric(Re(V * diag(1:/sqrt(d)) * V'))
		}
	}

	if (this.has_intercept) ::swap(x_means, this.x_means)
	::swap(inv_xx, this.inv_xx)
}


`Void' FE_Factor::mult(`Vector' coefs, `String' suggested_preconditioner, `Vector' ans)
{
	// TODO: OPTIMIZE (Q=1, num_slopes=1)
	// BUGBUG TODO: except for corner cases (pf->levels==1?) we should be able to replace [pf->levels, .] with [pf->levels], which is much faster

	`Integer'				NI, NS, NN, FullK
	`Boolean'				intercept_only
	`Matrix'				tmp

	NI = num_levels * has_intercept
	NS = num_levels * num_slopes
	NN = NI + NS
	FullK = has_intercept + num_slopes
	assert(num_levels * FullK == NN)

	intercept_only = this.has_intercept& !this.num_slopes
	this.preconditioner = (FullK == 1 & suggested_preconditioner == "block_diagonal") ? "diagonal" : suggested_preconditioner

	if (this.preconditioner == "none") {
		if (intercept_only) {
			ans = ans + (this.num_levels > 1 ? coefs[this.levels] : coefs[this.levels, .])
		}
		else if (!this.has_intercept) {
			ans = ans + rowsum( inv_vec(coefs, FullK)[this.levels, .] :* this.unsorted_x )
		}
		else {
			ans = ans + rowsum(inv_vec(coefs, FullK)[this.levels, .] :* (J(this.num_obs, 1, 1), this.unsorted_x))
		}
	}
	else if (this.preconditioner == "diagonal") {
		if (this.has_intercept) {
			ans = ans + (coefs[|1\NI|] :* this.preconditioner_intercept)[this.levels, .]
		}
		if (this.num_slopes) {
			ans = ans + rowsum( (inv_vec(coefs[|NI+1\NN|], this.num_slopes) :* this.preconditioner_slopes)[this.levels, .] :* this.unsorted_x )
		}
	}
	else if (this.preconditioner == "block_diagonal") {
		assert_msg(this.has_intercept + this.num_slopes > 1) // no point in using block-diagonal with only one regressor
		FullK = has_intercept + num_slopes
		
		// Goal: Given matrix A (N,LK), matrix D (LK,LK), and vector z (LK, 1) we want to compute vector b (N,1)
		// Interpretation: given a coefficient vector 'z' we compute the predicted values i.e.  Ax = ADz

		// 1) Multiply D*z - since D is a block-diagonal LK,LK matrix, and we only store "inv_x" as LK,K
		//    we'll have to do a trick:
		//    Suppose each block has a matrix DD and a coef vector zz. Then,
		//    DD*zz is equivalent to a) expanding zz, b) dot-product, c) rowsum

		// a) Reshape the coef vector into a (L, K) coef matrix
		tmp = inv_vec(coefs, FullK)
		// b) Premultiply by D (to show why this works, work out the math for multiplying a single block)
		tmp = tmp # J(FullK, 1, 1)
		tmp = this.inv_xx :* tmp
		tmp = rowsum(tmp)
		assert(rows(tmp) == num_levels * FullK)
		assert(cols(tmp) == 1)

		// 2) Multiply A*tmp
		// a) We first need to reshape and multiply by the cvars (which means the matrix needs to be N*K)
		tmp = inv_vec(tmp, FullK)
		assert(rows(tmp) == num_levels)
		assert(cols(tmp) == FullK)
		tmp = tmp[this.levels, .]
		assert(rows(tmp) == this.num_obs)
		if (this.has_intercept) {
			tmp = tmp :* (J(this.num_obs, 1, 1) , this.unsorted_x)
		}
		else {
			tmp = tmp :* this.unsorted_x
		}
		// b) Aggregate adding up the contribution of the intercept and each slope
		tmp = rowsum(tmp)

		ans = ans + tmp
	}
	else {
		_error(3498, sprintf("invalid preconditioner %s", this.preconditioner))
	}
}


`Vector' FE_Factor::mult_transpose(`Vector' y, `String' suggested_preconditioner)
{
	`Integer'				NI, NS, N, FullK
	`Matrix' 				alphas
	`Boolean'				intercept_only


	NI = num_levels * has_intercept
	NS = num_levels * num_slopes
	N = NI + NS
	FullK = has_intercept + num_slopes
	
	intercept_only = this.has_intercept& !this.num_slopes
	this.preconditioner = (FullK == 1 & suggested_preconditioner == "block_diagonal") ? "diagonal" : suggested_preconditioner
	
	assert(num_levels * FullK == N)

	// A'y: sum of y over each group (intercept); weighted sum with slopes and/or weights
	// note that weights are already presorted

	if (this.preconditioner == "none") {

		if (intercept_only) {
			alphas = this.panelsum(y, 1)
		}
		else if (!this.has_intercept) {
			alphas = fast_vec( this.panelsum(y :* this.unsorted_x, 1) )	
		}
		else {
			alphas = fast_vec( this.panelsum(y :* (J(this.num_obs, 1, 1), this.unsorted_x), 1) )
		}


	}
	else if (this.preconditioner == "diagonal") {
		alphas = J(N, 1, .)
		if (this.has_intercept) {
			alphas[|1\NI|] = this.panelsum(y, 1) :* this.preconditioner_intercept
		}
		if (this.num_slopes) {
			alphas[|NI+1\N|] = fast_vec( this.panelsum(y :* this.unsorted_x, 1) :* this.preconditioner_slopes )
		}
	}
	else if (this.preconditioner == "block_diagonal") {
		assert_msg(this.has_intercept + this.num_slopes > 1) // no point in using block-diagonal with only one regressor

		// Goal: given matrix 'A', block-diag matrix 'D', and data vector 'y', compute D'A'y=DA'y
		// TODO: fix the case with no intercept!

		// 1) Compute A'y (the alphas)
		if (this.has_intercept) {
			alphas = this.panelsum( (y, y :* this.unsorted_x) , 1)
		}
		else {
			alphas = this.panelsum(y :* this.unsorted_x, 1)
		}

		assert(rows(alphas) == num_levels)
		assert(cols(alphas) == FullK)

		// 2) Premultiply by 'D'
		alphas = alphas # J(FullK, 1, 1)
		assert(rows(alphas) == FullK * num_levels)
		assert(cols(alphas) == FullK)
		alphas = this.inv_xx :* alphas
		alphas = rowsum(alphas)
		assert(rows(alphas) == num_levels * FullK)
		assert(cols(alphas) == 1)
	}
	else {
		_error(3498, sprintf("invalid preconditioner %s", this.preconditioner))
	}
	return(alphas)
}



`Void' FE_Factor::undo_preconditioning(`Matrix' alphas)
{
	`Integer'				FullK, num_rows, num_cols
	`Boolean'				intercept_only
	
	assert_in(this.preconditioner, ("none", "diagonal", "block_diagonal"))

	FullK = has_intercept + num_slopes

	if (this.preconditioner == "none") {
		// pass (nothing to do)
	}
	else if (this.preconditioner == "diagonal") {
		if (this.has_intercept & !this.num_slopes) {
			alphas = alphas :* this.preconditioner_intercept
		}
		else if (!this.has_intercept & this.num_slopes) {
			alphas = alphas :* this.preconditioner_slopes
		}
		else {
			alphas = alphas :* (this.preconditioner_intercept , this.preconditioner_slopes ) 
		}
	}
	else if (this.preconditioner == "block_diagonal") {
		num_rows = rows(alphas)
		num_cols = cols(alphas)
		alphas = alphas # J(FullK, 1, 1)
		alphas = alphas :* this.inv_xx
		alphas = rowsum(alphas)
		alphas = inv_vec(alphas, num_cols)
		
		assert(num_rows == rows(alphas))
		assert(num_cols == cols(alphas))
	}

}


`Void' FE_Factor::cleanup_before_saving()
{
	super.cleanup_before_saving()
	//FG.cleanup_before_saving()
	//FI.cleanup_before_saving()
}


`FE_Factor' fe_factor(`Varlist' varnames,
					| `DataCol' touse, // either string varname or a numeric index
					  `Boolean' verbose,
					  `String' method,
					  `Boolean' sort_levels,
					  `Boolean' count_levels,
					  `Integer' hash_ratio,
					  `Boolean' save_keys)
{
	`Factor'			F
	`FE_Factor'			FE_F

	if (args()<2) touse = ""
	F = factor(varnames, touse, verbose, method, sort_levels, count_levels, hash_ratio, save_keys)
	F.swap(FE_F)
	return(FE_F)
}


`FE_Factor' _fe_factor(`DataFrame' data,
               | `Boolean' integers_only,
                 `Boolean' verbose,
                 `String' method,
                 `Boolean' sort_levels,
                 `Boolean' count_levels,
                 `Integer' hash_ratio,
                 `Boolean' save_keys,
                 `Varlist' vars, 			// hack
                 `DataCol' touse)		 	// hack
{
	`Factor'			F
	`FE_Factor'			FE_F
	F = _factor(data, integers_only, verbose, method, sort_levels, count_levels, hash_ratio, save_keys, vars, touse)
	F.swap(FE_F)
	return(FE_F)
}

`FE_Factor' fe_join_factors(`FE_Factor' F1,
                      `FE_Factor' F2, 
                    | `Boolean' count_levels,
                      `Boolean' save_keys,
                      `Boolean' levels_as_keys)
{
	//`Factor'			F
	//`FE_Factor'			FE_F
	//F = join_factors(F1, F2, count_levels, save_keys, levels_as_keys)
	//FE_F = FE_Factor()
	//F.swap(FE_F)
	//return(FE_F)
	return(join_factors(F1, F2, count_levels, save_keys, levels_as_keys))
}

end
