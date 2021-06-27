// --------------------------------------------------------------------------
// Mata Code: regressions with individual FEs and team-level outcomes
// --------------------------------------------------------------------------
// Project URL: ...


// Miscellanea --------------------------------------------------------------
	mata: mata set matastrict on
	mata: mata set mataoptimize on


// Include ftools -----------------------------------------------------------
    cap findfile "ftools.mata"
    if (_rc) {
        di as error "reghdfe requires the {bf:ftools} package, which is not installed"
        di as error `"    - install from {stata ssc install ftools:SSC}"'
        di as error `"    - install from {stata `"net install ftools, from("https://github.com/sergiocorreia/ftools/raw/master/src/")"':Github}"'
        exit 9
    }
    include "`r(fn)'"

// Type aliases -------------------------------------------------------------
	loc Factor				class Factor scalar
	loc OptionalFactor      class Factor
	loc Factors				class Factor rowvector
	loc FactorPointer   	pointer(`Factor') scalar
	loc FactorPointers   	pointer(`Factor') rowvector
	
	loc FE_Factor			class FE_Factor scalar
	loc Optional_FE_Factor	class FE_Factor
	loc FE_Factors			class FE_Factor rowvector
	loc FE_FactorPointer   	pointer(`FE_Factor') scalar
	loc FE_FactorPointers   pointer(`FE_Factor') rowvector
	
	loc Individual_Factor	class Individual_Factor scalar
	//loc Optional_FE_Factor	class FE_Factor
	//loc FE_Factors			class FE_Factor rowvector
	//loc FE_FactorPointer   	pointer(`FE_Factor') scalar
	//loc FE_FactorPointers   pointer(`FE_Factor') rowvector
	
	loc FixedEffects 		class FixedEffects scalar
	loc BipartiteGraph  	class BipartiteGraph scalar
	loc Solution			class Solution scalar

	loc DEBUG				"" // Set to empty string to run slow debugging code, set to "if (0)" for faster runtime

// Includes -----------------------------------------------------------------
	* Notes on include order:
	* - Include "common.mata" first
	* - Include class definitions ("FE.mata") before constructors ("FE_constructor.mata")
	* - Include class definitions before functions that use it ("FE.mata" before "Regression.mata", "DoF.mata", "LSMR.mata", etc)
	//loc includes Mock_Matrix Common Solution FE FE_Constructor Regression Bipartite DoF LSMR LSQR
	

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
// Bipartite Graphs ---------------------------------------------------------
// - For simplicity, assume the graph represent (firm, ceo) pairs
// - TODO: Check when we don't need all these objects anymore and clean them up!

mata:

class BipartiteGraph
{
	// Computed by init()
	`Boolean'				verbose
	`Integer'				N		// Num. obs
	`Integer'				N1		// Num. levels of FE 1
	`Integer'				N2		// Num. levels of FE 2
	`Integer'				N12		// N1 + N2
	`FE_FactorPointer'		PF1		// we need .has_weights so this needs to be of class FE_Factor instead of Factor
	`FE_FactorPointer'		PF2		// we need .has_weights so this needs to be of class FE_Factor instead of Factor
	`Factor'				F12
	`Factor'				F12_1
	`Factor'				F12_2

	// Computed by init_zigzag()
	`Vector' 				queue
	`Vector' 				stack
	`Vector' 				keys1_by_2
	`Vector' 				keys2_by_1
	`Integer'				num_subgraphs
	`Variable'				subgraph_id		// (optional)

	// Computed by compute_cores()
	`Vector' 				cores
	`Vector' 				drop_order

	// Computed after prune_1core()
	`Integer'               N_drop
	`Variable'              mask            // mask (0|1) of obs that are dropped after prunning of degree-1 edges
	`Boolean'               prune           // Whether to recursively prune degree-1 edges
	`Vector'                drop2idx
	`Matrix'                drop2info
	`Variable'              sorted_w
	`Boolean'				has_weights
	`Variable'              sorted_true_weight

	// Methods
	`Void'					init()
	`Real'					init_zigzag()
	`Void'					compute_cores()
	`Void'					prune_1core()
	`Variables'				expand_1core()
	`Variables'				partial_out()
	`Variables'				__partial_out_map()
	`Variables'				__partial_out_laplacian()

	`Void'					cleanup_before_saving()
}
	

`Void' BipartiteGraph::init(`FE_FactorPointer' PF1,
                            `FE_FactorPointer' PF2,
                            `Boolean' verbose)
{
	if (verbose>0) {
		printf("\n{txt}# Initializing bipartite graph\n")
		printf("    - FE #1: {res}%s{txt}\n", invtokens((*PF1).varlist))
		printf("    - FE #2: {res}%s{txt}\n", invtokens((*PF2).varlist))
	}
	this.verbose = verbose
	this.PF1 = PF1
	this.PF2 = PF2

	assert( (*PF1).num_obs == (*PF2).num_obs )

	N = (*PF1).num_obs
	N1 = (*PF1).num_levels
	N2 = (*PF2).num_levels

	N12 = N1 + N2
	(*PF1).panelsetup() // Just in case
	(*PF2).panelsetup() // Just in case

	// F12 must be created from F1.levels and F2.levels (not from the original keys)
	// This is set automatically by join_factors() with the correct flag:
	//			F12 = join_factors(F1, F2, ., ., 1)
	// But you can also run (slower)
	//			F12 = _fe_factor( (F1.levels, F2.levels) )
	//			asarray(F12.extra, "levels_as_keys", 1)
	if (verbose>0) printf("{txt}   - computing F12:  ")
	// join_factors(F1, (*PF2) [, count_levels, save_keys, levels_as_keys])
	F12 = fe_join_factors((*PF1), (*PF2), ., ., 1)
	if (verbose>0) printf("{txt} edges found: {res}%-10.0gc{txt}\n", F12.num_levels)
	F12.panelsetup()
	
	if (verbose>0) printf("{txt}   - computing F12_1:")
	// _factor(data [, integers_only, verbose, method, sort_levels, count_levels, hash_ratio, save_keys])
	F12_1 = _fe_factor(F12.keys[., 1], 1, 0, "", 1, 1, ., 0)
	if (verbose>0) printf("{txt} edges found: {res}%-10.0gc{txt}\n", F12_1.num_levels)
	F12_1.panelsetup()
	
	if (verbose>0) printf("{txt}   - computing F12_2:")
	F12_2 = _fe_factor(F12.keys[., 2], 1, 0, "", 1, 1, ., 0)
	if (verbose>0) printf("{txt} edges found: {res}%-10.0gc{txt}\n", F12_2.num_levels)
	F12_2.panelsetup()
}


// --------------------------------------------------------------------------
// init_zigzag()
// --------------------------------------------------------------------------
//   Construct -queue- and -stack- vectors that allow zigzag iteration
//   
//   queue: firm and CEOs that will be processed, in the required order
//   	   note: negative values indicate CEOs
//   
//   stack: for each firm/CEO, the list of nodes it connects to
//          note: stacks are zero-separated
//   
// --------------------------------------------------------------------------
//   As a byproduct, also computes the number of disjoint subgraphs.
//   See the algorithm from on Abowd, Creecy and Kramarz (WP 2002) p4. Sketch:
//
//   		g = 0
//   		While there are firms w/out a group:
//   		     g++
//   		     Assign the first firm w/out a group to group g
//   		     Repeat until no further changes:
//   		         Add all persons employed by a firm in g to g
//   		         Add all firms that employ persons in g to g
//   		 return(g)
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
`Real' BipartiteGraph::init_zigzag(| `Boolean' save_subgraphs)
{
	`Vector'				counter1
	`Vector'				counter2
	`Vector' 				done1
	`Vector' 				done2

	`Integer'				i_stack		// use to process the queue
	`Integer'				last_i		// use to fill out the queue
	`Integer'				start_j		// use to search for firms to start graph enumeration
	`Integer'				i_queue
	`Integer'				id 			// firm number if id>0; error if id=0; ceo number if id<0
	`Integer'				j 			// firm # (or viceversa)
	`Integer'				k 			// ceo # (or viceversa)
	`Integer'				c 			// temporary counter
	`Integer'				i 			// temporary iterator
	
	`Matrix'				matches // list of CEOs that matched with firm j (or viceversa)

	if (verbose>0) printf("\n{txt}# Initializing zigzag iterator for bipartite graph\n")
	assert(F12_1.panel_is_setup == 1)
	assert(F12_2.panel_is_setup == 1)
	//assert(asarray(F12.extra, "levels_as_keys") == 1)

	// If subgraph_id (mobility groups) is anything BUT zero, we will save them
	if (args()==0 | save_subgraphs==.) save_subgraphs = 0
	if (save_subgraphs) {
		subgraph_id = J(N2, 1, .)
	}

	queue = J(N12, 1, 0)
	stack = J(F12.num_levels + N12, 1, .) // there are N12 zeros
	counter1 = J(N1, 1, 0)
	counter2 = J(N2, 1, 0)

	keys1_by_2 = F12_2.sort(F12.keys[., 1])
	keys2_by_1 = F12_1.sort(F12.keys[., 2])
	done1 = J(N1, 1, 0) // if a firm is already on the queue
	done2 = J(N2, 1, 0) // if a CEO is already on the queue

	// Use -j- for only for firms and -k- only for CEOs
	// Use -i_queue- to iterate over the queue and -i_stack- over the stack
	// Use -last_i- to fill out the queue (so its the last filled value)
	// Use -i- to iterate arbitrary vectors
	// Use -id- to indicate a possible j or k (negative for k)
	// Use -start_j- to remember where to start searching for new subgraphs

	i_stack = 0
	last_i = 0
	start_j = 1
	num_subgraphs = 0

	for (i_queue=1; i_queue<=N12; i_queue++) {
		id = queue[i_queue] // >0 if firm ; <0 if CEO; ==0 if nothing yet
		j = k = . // just to avoid bugs
		
		// Pick starting point (useful if the graph is disjoint!)
		if (id == 0) {
			assert(last_i + 1 == i_queue)
			for (j=start_j; j<=N1; j++) {
				if (!done1[j]) {
					queue[i_queue] = id = j
					start_j = j + 1
					++last_i
					break
				}
			}
			// printf("{txt} - starting subgraph with firm %g\n", j)
			++num_subgraphs
			assert(id != 0) // Sanity check
		}

		if (id > 0) {
			// It's a firm
			j = id
			done1[j] = 1
			matches = panelsubmatrix(keys2_by_1, j, F12_1.info)
			for (i=1; i<=rows(matches); i++) {
				k = matches[i]
				c = counter2[k]
				counter2[k] = c + 1
				if (!done2[k]) {
					if (!c) {
						queue[++last_i] = -k
					}
				 	stack[++i_stack] = k
				}
			}
			stack[++i_stack] = 0
		}
		else {
			// It's a CEO
			k = -id
			done2[k] = 1
			matches = panelsubmatrix(keys1_by_2, k, F12_2.info)
			for (i=1; i<=rows(matches); i++) {
				j = matches[i]
				c = counter1[j]
				counter1[j] = c + 1
				if (!done1[j]) {
					if (!c) {
						queue[++last_i] = j
					}
					stack[++i_stack] = j
				}
			}
			stack[++i_stack] = 0
			if (save_subgraphs) subgraph_id[k] = num_subgraphs
		}
	}

	// Sanity checks
	assert(counter1 == F12_1.counts)
	assert(counter2 == F12_2.counts)
	assert(!anyof(queue, 0)) // queue can't have zeros at the end
	assert(allof(done1, 1))
	assert(allof(done2, 1))
	assert(!missing(queue))
	assert(!missing(stack))

	if (save_subgraphs) subgraph_id = subgraph_id[(*PF2).levels]
	
	if (verbose>0) printf("{txt}   - disjoint subgraphs found: {res}%g{txt}\n", num_subgraphs)
	return(num_subgraphs)
}


// --------------------------------------------------------------------------
// compute_cores()
// --------------------------------------------------------------------------
// 	  Computes vertex core numbers, which allows k-core pruning
//	  Algorithm used is listed here: https://arxiv.org/abs/cs/0310049
// --------------------------------------------------------------------------
// Note:
//    maybe use the k-cores for something useful? eg:
//    we might want to weight the core numbers by the strength (# of obs together)
//    https://arxiv.org/pdf/1611.02756.pdf --> # of butterflies in bipartite graph
//    this paper also has useful data sources for benchmarks
//    # of primary and secondary vertices, edges
// --------------------------------------------------------------------------

`Void' BipartiteGraph::compute_cores()
{
	`Factor'				Fbin
	`Boolean'				is_firm
	`Integer'				M, ND, j, jj
	`Integer'				i_v, i_u, i_w
	`Integer'				pv, pu, pw
	`Integer'				v, u, w
	`Integer'				dv, du
	`Vector'				bin, deg, pos, invpos, vert, neighbors

	if (verbose>0) printf("{txt}# Computing vertex core numbers\n\n")

	// v, u, w are vertices; <0 for CEOs and >0 for firms
	// vert is sorted by degree; deg is unsorted
	// pos[i] goes from sorted to unsorted, so:
	// 		vert[i] === original_vert[ pos[i] ]
	// invpos goes from unsorted to sorted, so:
	//		vert[invpos[j]] === original_vert[j]

	// i_u represents the pos. of u in the sorted tables
	// pu represents the pos. of u in the unsorted/original tables

	assert_msg(F12_1.panel_is_setup==1, "F12_1 not set up")
	assert_msg(F12_2.panel_is_setup==1, "F12_2 not set up")
	assert_msg(rows(queue)==N12, "Wrong number of rows in queue")
	assert_msg(rows(keys1_by_2)==F12.num_levels, "Wrong number of rows in keys1")
	assert_msg(rows(keys2_by_1)==F12.num_levels, "Wrong number of rows in keys2")

	deg = F12_1.counts \ F12_2.counts
	ND = max(deg) // number of degrees

	Fbin = _factor(deg, 1, 0)
	Fbin.panelsetup()

	bin = J(ND, 1, 0)
	bin[Fbin.keys] = Fbin.counts
	bin = rows(bin) > 1 ? runningsum(1 \ bin[1..ND-1]) : 1
	
	pos = Fbin.p
	invpos = invorder(Fbin.p)

	vert = Fbin.sort(( (1::N1) \ (-1::-N2) ))

	for (i_v=1; i_v<=N12; i_v++) {
		v = vert[i_v]
		is_firm = (v > 0)

		neighbors = is_firm ? panelsubmatrix(keys2_by_1, v, F12_1.info) : panelsubmatrix(keys1_by_2, -v, F12_2.info)
		M = rows(neighbors)
		
		for (j=1; j<=M; j++) {	
			pv = pos[i_v]
			jj = neighbors[j]
			pu = is_firm ? N1 + jj : jj // is_firm is *not* for the neighbor
			dv = deg[pv]
			du = deg[pu]
		
			if (dv < du) {
				i_w = bin[du]
				w = vert[i_w]
				u = is_firm ? -jj : jj // is_firm is *not* for the neighbor
				if (u != w) {
					pw = pos[i_w]
					i_u = invpos[pu]
					pos[i_u] = pw
					pos[i_w] = pu
					vert[i_u] = w
					vert[i_w] = u
					invpos[pu] = i_w
					invpos[pw] = i_u
				}
				bin[du] = bin[du] + 1
				deg[pu] = deg[pu] - 1
			}
		} // end for neighbor u (u ~ v)
	} // end for each node v

	if (verbose>0) {
		//printf("{txt}      Table: core numbers and vertex count\n")
		Fbin = _factor(deg, 1, 0)
		//printf("\n")
		mm_matlist(Fbin.counts, "%-8.0gc", 2, strofreal(Fbin.keys), "Freq.", "Core #")
		printf("\n")
	}

	// ((F1.keys \ F2.keys), (F12_1.keys \ -F12_2.keys))[selectindex(deg:==1), .]		
	
	// Store the values in the class
	swap(this.drop_order, vert)
	swap(this.cores, deg)
}

// --------------------------------------------------------------------------
// prune_1core()
// --------------------------------------------------------------------------
//    Prune edges with degree-1
//    That is, recursively remove CEOs that only worked at one firm,
//    and firms that only had one CEO in the sample, until every agent
//    in the dataset has at least two matches
// --------------------------------------------------------------------------
`Void' BipartiteGraph::prune_1core(| `Variable' weight)
{
	`Integer'               N_drop2, i, j, i1, i2, j1, j2, K_drop2
	`Vector'                drop1, drop2
	`Vector'                tmp_mask
	`Vector'                proj1, proj2
	`Variable'              w, tmp_weight

	has_weights = (args()>0 & rows(weight) > 1)
	if (has_weights) sorted_true_weight = (*PF1).sort(weight)
	tmp_weight = has_weights ? weight : J(N, 1, 1)

	N_drop = sum(cores :== 1)
	if (!N_drop) {
	    if (verbose>0) printf("{txt}   - no 1-core vertices found\n")
	    prune = 0
	    return
	}
	if (verbose>0) printf("{txt}   - 1-core vertices found: {res}%g{txt}\n", N_drop)

	drop_order = drop_order[1..N_drop]
	drop1 = selectindex(cores[1..N1] :== 1)
	cores = .
	drop1 = (1::N1)[drop1]
	drop2 = -select(drop_order, drop_order:<0)
	
	K_drop2 = rows(drop2)
	N_drop2 = K_drop2 ? sum((*PF2).info[drop2, 2] :- (*PF2).info[drop2, 1] :+ 1) : 0

	tmp_mask = J(N1, 1, 0)
	if (rows(drop1)) tmp_mask[drop1] = J(rows(drop1), 1, 1)
	mask = tmp_mask[(*PF1).levels, 1]
	tmp_mask = J(N2, 1, 0)
	if (K_drop2) tmp_mask[drop2] = J(K_drop2, 1, 1)
	mask = mask :| tmp_mask[(*PF2).levels, 1]
	tmp_mask = .

	drop2idx = J(N_drop2, 1, .)
	drop2info = J(N2, 2, .)

	j1 = 1
	for (i=1; i<=K_drop2; i++) {
	    j = drop2[i]
	    i1 = (*PF2).info[j, 1]
	    i2 = (*PF2).info[j, 2]

	    j2 = j1 + i2 - i1
	    drop2idx[j1::j2] = i1::i2
	    drop2info[j, .] = (j1, j2)
	    j1 = j2 + 1
	}

	if (!(*PF2).is_sorted) {
	    assert(((*PF2).p != J(0, 1, .)))
	    drop2idx = (*PF2).p[drop2idx, .]
	}

	if (!(*PF1).is_sorted) {
	    assert(((*PF1).inv_p != J(0, 1, .)))
	    drop2idx = invorder((*PF1).p)[drop2idx, .]
	}

	// To undo pruning, I need (*PF1).info[drop1, .] & drop2info & drop2idx

	// Set weights of pruned obs. to zero
	tmp_weight[selectindex(mask)] = J(sum(mask), 1, 0)

	// Update sorted weights for g=1,2
	PF1->set_weights(tmp_weight)
	PF2->set_weights(tmp_weight)
	tmp_weight = . // cleanup

	// Select obs where both FEs are degree-1 (and thus omitted)
	sorted_w = J(N, 1, 1)
	
	proj1 = PF1->panelmean(sorted_w)[PF1->levels, .]
	proj2 = PF2->panelmean(sorted_w)[PF2->levels, .]

	sorted_w = ((sorted_w - proj1) :!= 1) :| ((sorted_w - proj2) :!= 1)
	proj1 = proj2 = .
	sorted_w = (*PF1).sort(sorted_w)

	prune = 1
}

// --------------------------------------------------------------------------
// prune_1core()
// --------------------------------------------------------------------------
//    Prune edges with degree-1
//    That is, recursively remove CEOs that only worked at one firm,
//    and firms that only had one CEO in the sample, until every agent
//    in the dataset has at least two matches
// --------------------------------------------------------------------------
`Variables' BipartiteGraph::expand_1core(`Variables' y)
{
    `Boolean'               zero_weights
    `Variable'              sorted_y
    `Integer'               i, j, j1, j2, i2, k1, k2, nk
    `Matrix'                tmp_y
    `Vector'                tmp_w, tmp_idx, new_w
    `RowVector'             tmp_mean

    if (prune==0) return(y)
    if (verbose>0) printf("{txt}# Expanding 2-core into original dataset\n\n")
    assert(N_drop == rows(drop_order))

    sorted_y = (*PF1).sort(y)

    i2 = 0
    for (i=N_drop; i>=1; i--) {
        j = drop_order[i]
        if (j > 0) {
            j1 = (*PF1).info[j, 1]
            j2 = (*PF1).info[j, 2]

            tmp_y = sorted_y[| j1 , 1 \ j2 , . |] // panelsubmatrix(sorted_y, j, (*PF1).info)
            tmp_w = sorted_w[|j1, 1 \ j2, .|] // panelsubmatrix(sorted_w, j, (*PF1).info)
            new_w = has_weights ? sorted_true_weight[|j1, 1 \ j2, .|] : J(j2-j1+1, 1, 1)
            zero_weights = !sum(tmp_w)
            if (!zero_weights) {
                tmp_mean = mean(tmp_y, tmp_w)
                assert(!missing(tmp_mean)) // bugbug remove later
                sorted_y[| j1 , 1 \ j2 , . |] = tmp_y :- tmp_mean
            }
            sorted_w[| j1 , 1 \ j2 , 1 |] = new_w
        }
        else {
            ++i2
            j1 = drop2info[-j, 1]
            j2 = drop2info[-j, 2]
            tmp_idx = drop2idx[| j1 , 1 \ j2 , 1 |]
            tmp_y = sorted_y[tmp_idx, .]
            tmp_w = sorted_w[tmp_idx]
            zero_weights = !sum(tmp_w)
            if (zero_weights) {
                tmp_w = has_weights ? sorted_true_weight[tmp_idx] : J(j2-j1+1, 1, 1)
            }
            tmp_mean = mean(tmp_y, tmp_w)
            assert(!missing(tmp_mean)) // bugbug remove later
            nk = rows(tmp_idx)
            for (k1=1; k1<=nk; k1++) {
               k2 = tmp_idx[k1]
               sorted_y[k2, .] = sorted_y[k2, .] - tmp_mean
               sorted_w[k2] = has_weights ? sorted_true_weight[k2] : 1
            }
        }
    }

    if (verbose>0) printf("{txt}   - number of coefficients solved triangularly: {res}%s{txt}\n", strofreal(rows(drop_order)))
    return((*PF1).invsort(sorted_y))
}


`Variables' BipartiteGraph::partial_out(`Variables' y)
{

}


`Variables' BipartiteGraph::__partial_out_map(`Variables' y)
{

}


`Variables' BipartiteGraph::__partial_out_laplacian(`Variables' y)
{

}


`Void' BipartiteGraph::cleanup_before_saving()
{
	this.F12.cleanup_before_saving()
	this.F12_1.cleanup_before_saving()
	this.F12_2.cleanup_before_saving()

	// (*this.PF1).cleanup_before_saving()
	// (*this.PF2).cleanup_before_saving()

	// We need to set this to NULL
	// This potentially prevents us from using this after reloading, but then we don't need this for partialling out
	PF1 = NULL
	PF2 = NULL
}

end
// --------------------------------------------------------------------------
// Code that extends Factor() for individual FEs in a group setting
// --------------------------------------------------------------------------

mata:

class Individual_Factor extends FE_Factor
{
	// Main properties
	`String' 				aggregate_function // how do we add up the FEs (mean, sum)
	`FE_Factor'				FG, FI
	`BipartiteGraph'		bg
	`Integer'				verbose
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
	virtual `Void'			set_weights()

	virtual `Void'			drop_obs()  	// x = A'b , implemented as x = A.mult_transpose(b)
	virtual `Void'			panelsetup()
	`DataFrame'				sort()

	virtual `Void'			cleanup_before_saving()
}


`Void' Individual_Factor::new()
{
	is_individual_fe = 1
}


`Vector' Individual_Factor::drop_singletons(`Vector' sample, `Vector' indiv_sample) // , | `Vector' fweight, `Boolean' zero_threshold)
{
	// If an individual appears only once then his FE can perfectly explain the group outcome
	// Then, we should drop all the obs for this group, which can lead to more individuals appearing only once
	
	// When we do this recursively, it's equivalently to
	// a) setting all counts of the bg.F12_1 factor (based on FG) equal to 2
	// b) computing the core numbers of each node/vertex
	// To see why this works, look at the algorithm on p2 section 3 of https://arxiv.org/pdf/cs/0310049.pdf
	// This algo starts by setting core=degree and sorting by it (note that degree is bg.F12_1.counts)
	// It then looks at the nodes of the lowest degree, which can only be individuals with degree=1 (call them singleton individuals)
	// Then it reduces the cores of all the groups connected to those singleton individuals, and continues
	// Because we set counts=2 for all groups, we just need one singleton individual to make the group have core 1
	// Thus, all groups with cores=1 should be dropped as they are singleton observations

	if (verbose>2) printf("{txt}   - running individual_factor.drop_singletons()\n")
	`Vector'				mask, drop_mask, idx, indiv_idx
	`Integer'				N_drop

	// Recall that 'bg' is the bipartite graph with group and individual IDs
	assert(bg.F12_1.panel_is_setup==1)
	assert(allof(bg.F12_1.counts, 2))
	bg.compute_cores() // update 'cores'
	assert(rows(bg.cores) == FG.num_levels + FI.num_levels)

	// Simpler alternative to prune_1core(), for what we want
	drop_mask = bg.cores[1..FG.num_levels] :== 1 // we only care about the cores of groups, not individuals (which are at the end)
	drop_mask = drop_mask[FG.levels]
	assert_msg(rows(drop_mask) == rows(indiv_sample))
	N_drop = sum(drop_mask)
	if (verbose>2) printf("{txt}   - 1-core vertices found (i.e. singleton individuals): {res}%g{txt}\n", N_drop)

	// Recall that 'indiv_sample' have all the obs for each group (with the individual IDs), while 'sample' only has one obs for each group
	// We have 'drop_mask' for each obs in 'indiv_sample'; we want to match it to 'sample'

	// dropped obs are those that were in 'sample' and are no longer in 'indiv_sample' (all obs in 'indiv_sample' should also be in 'sample')
	if (N_drop) {

		// drop_mask has the obs within indiv_sample that we want to drop, but we want to map it to the obs in the data

		// 'idx' are the obs that are present in the old 'sample' but missing in the new 'indiv_sample'
		// we want a mask that is 0 except for the obs that existed in sample but are now dropped
		mask = J(st_nobs(), 1, 0)
		mask[sample] = 1::rows(sample)
		mask = mask :* create_mask(st_nobs(), 0, select(indiv_sample, drop_mask), 1)
		idx = select(mask, mask)

		// Drop observations in 'indiv_idx'
		
		//mask = J(st_nobs(), 1, 0)
		//mask[indiv_sample, .] = indiv_sample
		//mask = mask :* create_mask(st_nobs(), 0, select(indiv_sample, drop_mask), 1)
		//indiv_idx = select(mask, mask)

		indiv_idx = selectindex(drop_mask)

		if (this.FG.num_levels == rows(idx)) {
			//exit(error(2001))
			printf("{err}insufficient observations (dropped %f singletons due to individual fixed effects; check input to indiv() option)\n", rows(idx))
			exit(2001)
		}
		drop_obs(indiv_idx)
	}
	else {
		idx = J(0, 1, .)
	}

	return(idx)
}


`DataFrame' Individual_Factor::sort(`DataFrame' data)
{
	return(this.FI.sort(data))
}


`Void' Individual_Factor::drop_obs(`Vector' idx)
{
	this.FG.drop_obs(idx)
	this.FI.drop_obs(idx)

	// Update bipartite graph (because we updated FG and FI)
	bg = BipartiteGraph()
	bg.init(&FG, &FI, verbose-2)
	(void) bg.init_zigzag(0) // 1 => save subgraphs into bg.subgraph_id
	assert(rows(bg.F12_1.counts) == bg.F12_1.num_levels)
	if (verbose>2) printf("\n{txt}# Tricking Bipartite() into dropping singletons instead of doing 1-core prunning\n\n")
	bg.F12_1.counts = J(bg.F12_1.num_levels, 1, 2)
	this.num_levels = this.FI.num_levels
	this.num_obs = this.FG.num_levels

	this.group_p = this.group_inv_p = . // ensure they are empty
}


`Void' Individual_Factor::panelsetup()
{
	// BUGBUG haven't tested this fully
	this.FG.panelsetup()
	this.FI.panelsetup()
	assert(this.bg.F12_1.panel_is_setup==1)
	assert(this.bg.F12_2.panel_is_setup==1)
}


`Matrix' Individual_Factor::panelsum(`Matrix' x, | `Boolean' sort_input, `Boolean' ignore_weights)
{	
	return(this.FI.panelsum(x, sort_input, ignore_weights))
}


`Void' Individual_Factor::set_weights(`Vector' input_weights, | `Boolean' verbose)
{
	if (args()<2 | verbose == .) verbose = 0
	if (verbose > 0) printf("{txt}   - sorting weights for indiv. factor {res}%s{txt}\n", this.absvar)
	
	assert_msg(this.has_weights == 0)
	assert_msg(this.weights == 1)
	assert_msg(this.weighted_counts == .)
	assert_msg(this.FG.weights == 1)
	assert_msg(this.FI.weights == 1)
	assert_msg(this.FI.num_obs == rows(input_weights))

	this.has_weights = this.FG.has_weights = this.FI.has_weights = 1
	this.weights = this.FG.weights = this.FI.weights = this.FI.sort(input_weights)
	this.weighted_counts = this.FI.panelsum(this.weights, 0, 1) // sort_input=0 ignore_weights=1
}


`Void' Individual_Factor::init_diagonal_preconditioner()
{
	`Matrix'				precond
	`Matrix'				sorted_x

	assert(is_individual_fe == 1)

	if (this.has_weights) {
		assert(rows(this.weighted_counts) == this.FI.num_levels)
	}

	if (this.has_intercept) {
		precond = this.has_weights ? this.weighted_counts : this.FI.counts
		this.preconditioner_intercept = sqrt(1 :/ precond)
	}

	if (this.num_slopes) {
		precond = this.panelsum(this.unsorted_x :^ 2, 1)
		this.preconditioner_slopes = 1 :/ sqrt(precond)
	}
}


`Void' Individual_Factor::init_block_diag_preconditioner(`String' weight_type, `Real' weights_scale, `Vector' true_weights)
{
	assert(0)
	// BUGBUG
}


`Void' Individual_Factor::mult(`Vector' coefs, `String' suggested_preconditioner, `Vector' ans)
{
	// Add up the individual coefs for each group

	`Integer'				NI, NS, NN, FullK
	`Boolean'				intercept_only
	`Matrix'				tmp
	`Vector'				expanded_coefs, sum_coefs

	NI = num_levels * has_intercept
	NS = num_levels * num_slopes
	NN = NI + NS
	FullK = has_intercept + num_slopes

	assert(num_levels == FI.num_levels)
	assert(num_levels * FullK == NN)
	assert(aggregate_function == "sum" | aggregate_function == "mean")
	assert(rows(coefs) == NN)
	assert(rows(ans) == FG.num_levels)

	intercept_only = this.has_intercept & !this.num_slopes
	this.preconditioner = (FullK == 1 & suggested_preconditioner == "block_diagonal") ? "diagonal" : suggested_preconditioner

	// REMOVE THIS LATER
	assert(this.preconditioner == "none" | preconditioner == "diagonal")

	if (this.preconditioner == "none") {

		if (intercept_only) {
			expanded_coefs = this.FI.num_levels > 1 ? coefs[this.FI.levels] : coefs[this.FI.levels, .]
			sum_coefs = this.FG.panelsum(expanded_coefs, 1)
			if (aggregate_function == "mean") {
				sum_coefs = sum_coefs :/ this.FG.counts
			}
		}
		else {
			expanded_coefs = inv_vec(coefs, FullK)[this.FI.levels, .]
			sum_coefs = rowsum(expanded_coefs :* (J(this.FG.num_obs, this.has_intercept, 1), this.unsorted_x))
			sum_coefs = this.FG.panelsum(sum_coefs, 1)
			if (aggregate_function == "mean") {
				sum_coefs = sum_coefs :/ this.FG.counts
			}
		}

		//else if (!this.has_intercept) {
		//	expanded_coefs = inv_vec(coefs, FullK)[this.FI.levels, .]
		//	sum_coefs = rowsum(expanded_coefs :* this.unsorted_x)
		//	sum_coefs = this.FG.panelsum(sum_coefs, 1)
		//	if (aggregate_function == "mean") {
		//		sum_coefs = sum_coefs :/ this.FG.counts
		//	}
		//}
		//else {
		//	expanded_coefs = inv_vec(coefs, FullK)[this.FI.levels, .]
		//	sum_coefs = rowsum(expanded_coefs :* (J(this.num_obs, 1, 1), this.unsorted_x))
		//	sum_coefs = this.FG.panelsum(sum_coefs, 1)
		//	if (aggregate_function == "mean") {
		//		sum_coefs = sum_coefs :/ this.FG.counts
		//	}
		//}
	}
	else if (this.preconditioner == "diagonal") {
		sum_coefs = J(this.num_obs, 1, 0)

		if (this.has_intercept) {
			expanded_coefs = (coefs[|1\NI|] :* this.preconditioner_intercept)[this.FI.levels, .]
			sum_coefs = sum_coefs + this.FG.panelsum(expanded_coefs, 1)
		}

		if (this.num_slopes) {
			expanded_coefs = inv_vec(coefs[|NI+1\NN|], this.num_slopes) // transform vector to matrix of coefs
			expanded_coefs = (expanded_coefs :* this.preconditioner_slopes)[this.FI.levels, .] :* this.unsorted_x
			if (this.num_slopes > 1) expanded_coefs = rowsum(expanded_coefs)
			sum_coefs = sum_coefs + this.FG.panelsum(expanded_coefs, 1)
		}

		if (aggregate_function == "mean") {
			sum_coefs = sum_coefs :/ this.FG.counts
		}
	}
	else if (this.preconditioner == "block_diagonal") {
		assert(0)
		assert_msg(this.has_intercept + this.num_slopes > 1) // no point in using block-diagonal with only one regressor
		FullK = has_intercept + num_slopes
		
		// Goal: Given matrix A (N,LK), matrix D (LK,LK), and vector z (LK, 1) we want to compute vector b (N,1)
		// Interpretation: given a coefficient vector 'z' we compute the predicted values i.e.  Ax = ADz

		// 1) Multiply D*z - since D is a block-diagonal LK,LK matrix, and we only store "inv_x" as LK,K
		//    we'll have to do a trick:
		//    Suppose each block has a matrix DD and a coef vector zz. Then,
		//    DD*zz is equivalent to a) expanding zz, b) dot-product, c) rowsum

		//			// a) Reshape the coef vector into a (L, K) coef matrix
		//			tmp = inv_vec(coefs, FullK)
		//			// b) Premultiply by D (to show why this works, work out the math for multiplying a single block)
		//			tmp = tmp # J(FullK, 1, 1)
		//			tmp = this.inv_xx :* tmp
		//			tmp = rowsum(tmp)
		//			assert(rows(tmp) == num_levels * FullK)
		//			assert(cols(tmp) == 1)
		//			
		//			// 2) Multiply A*tmp
		//			// a) We first need to reshape and multiply by the cvars (which means the matrix needs to be N*K)
		//			tmp = inv_vec(tmp, FullK)
		//			assert(rows(tmp) == num_levels)
		//			assert(cols(tmp) == FullK)
		//			tmp = tmp[this.levels, .]
		//			assert(rows(tmp) == this.num_obs)
		//			tmp = tmp :* (J(this.num_obs, 1, 1) , this.unsorted_x)
		//			// b) Aggregate adding up the contribution of the intercept and each slope
		//			tmp = rowsum(tmp)
		//			
		//			ans = ans + tmp
	}
	else {
		_error(3498, sprintf("invalid preconditioner %s", this.preconditioner))
	}

	ans = ans + (this.FG.is_sorted ? sum_coefs : sum_coefs[group_inv_p])
}


`Vector' Individual_Factor::mult_transpose(`Vector' y, `String' suggested_preconditioner)
{
	// A'y returns the sum of individual FEs within each group +- +- +-
	// 'y' is a variable i.e. at the group level
	// A'y multiplies 'y' with each row of A i.e. for each individual reports the sum of 'y' over the groups it took part in
	// With an average function, then we need to compute the average because A is not 0/1

	`Integer'				NI, NS, N, FullK
	`Matrix' 				alphas
	`Boolean'				intercept_only
	`Vector'				sorted_y, expanded_y

	NI = num_levels * has_intercept
	NS = num_levels * num_slopes
	N = NI + NS
	FullK = has_intercept + num_slopes
	
	intercept_only = this.has_intercept & !this.num_slopes
	this.preconditioner = (FullK == 1 & suggested_preconditioner == "block_diagonal") ? "diagonal" : suggested_preconditioner
	
	assert_msg(num_levels * FullK == N, "error1")
	assert_msg(num_levels == this.FI.num_levels, sprintf("assertion failed: num_levels (%f) ≠ FI.num_levels (%f)", num_levels, FI.num_levels))
	assert_msg(FG.num_levels == rows(y), "error3", 3456)
	assert_msg(aggregate_function == "sum" | aggregate_function == "mean", "error4")

	// REMOVE THESE LATER
	assert_msg(this.preconditioner == "none" | this.preconditioner == "diagonal", "error5")

	// Sort input variable by group/team
	sorted_y = this.FG.is_sorted ? y : y[this.group_p]

	// A'y: sum of y over each group (intercept); weighted sum with slopes and/or weights
	// note that weights are already presorted

	if (aggregate_function == "sum") {
		expanded_y = sorted_y[this.FG.levels, .]
	}
	else if (aggregate_function == "mean") {
		// If aggregation=mean, we find out denom and divide
		expanded_y = (sorted_y :/ this.FG.counts)[this.FG.levels, .]
	}
	else {
		_error(3333, sprintf("aggregation function not implemented: %s\n", aggregate_function))
	}

	if (this.preconditioner == "none") {
		if (intercept_only) {
			alphas = this.FI.panelsum(expanded_y, 1)
		}
		else if (!this.has_intercept) {
			alphas = fast_vec(this.panelsum(expanded_y :* this.unsorted_x, 1))
		}
		else {
			alphas = fast_vec( this.panelsum(expanded_y :* (J(this.num_obs, 1, 1), this.unsorted_x), 1) )
		}
	}
	else if (this.preconditioner == "diagonal") {
		alphas = J(N, 1, .)

		if (this.has_intercept) {
			alphas[|1\NI|] = this.FI.panelsum(expanded_y, 1) :* this.preconditioner_intercept
		}
		if (this.num_slopes) {
			alphas[|NI+1\N|] = fast_vec( this.panelsum(expanded_y :* this.unsorted_x, 1) :* this.preconditioner_slopes )
		}
	}
	else if (this.preconditioner == "block_diagonal") {
		assert(0)
		// assert_msg(this.has_intercept + this.num_slopes > 1) // no point in using block-diagonal with only one regressor
		// 
		// // Goal: given matrix 'A', block-diag matrix 'D', and data vector 'y', compute D'A'y=DA'y
		// // TODO: fix the case with no intercept!
		// 
		// // 1) Compute A'y (the alphas)
		// alphas = panelsum(this.sort((y, y :* this.unsorted_x)), this.weights, this.info)
		// assert(rows(alphas) == num_levels)
		// assert(cols(alphas) == FullK)
		// 
		// // 2) Premultiply by 'D'
		// alphas = alphas # J(FullK, 1, 1)
		// assert(rows(alphas) == FullK * num_levels)
		// assert(cols(alphas) == FullK)
		// alphas = this.inv_xx :* alphas
		// alphas = rowsum(alphas)
		// assert(rows(alphas) == num_levels * FullK)
		// assert(cols(alphas) == 1)
	}
	else {
		_error(3498, sprintf("invalid preconditioner %s", this.preconditioner))
	}
	return(alphas)
}


`Void' Individual_Factor::undo_preconditioning(`Matrix' alphas)
{
	assert(0)
}



`Void' Individual_Factor::cleanup_before_saving()
{
	super.cleanup_before_saving()
	FG.cleanup_before_saving()
	FI.cleanup_before_saving()

	bg = BipartiteGraph() // we don't need this for partialling out; reset it to save space
	bg.cleanup_before_saving()
}

// indiv_factor(group_id, individual_id, this.sample, this.indiv_sample)
`Individual_Factor' indiv_factor(`Varname' group_id,
								 `Varname' individual_id,
								 `Vector' sample,
								 `Vector' indiv_sample,
								 `String' aggregate_function,
								 `Boolean' verbose)
{
	`Factor'				FG, FI
	`Individual_Factor'		F
	`BipartiteGraph'        bg

	assert(aggregate_function == "sum" | aggregate_function == "mean")

	// Note that fweights are meaningless here because data cannot be duplicated...

	if (verbose>2) printf("\n{txt}# Constructing factor for group id: {res}%s{txt}\n", group_id)
	FG = fe_factor(group_id, indiv_sample, ., "", ., 1, ., 0)
	//if (verbose>2) printf("{txt} groups found: {res}%-10.0gc{txt}\n", FG.num_levels)

	if (verbose>2) printf("\n{txt}# Constructing factor for individual id: {res}%s{txt}\n", group_id)
	FI = fe_factor(individual_id, indiv_sample, ., "", ., 1, ., 0)
	//if (verbose>2) printf("{txt} individuals found: {res}%-10.0gc{txt}\n", FI.num_levels)

	bg = BipartiteGraph()
	bg.init(&FG, &FI, verbose-2)
	(void) bg.init_zigzag(0) // 1 => save subgraphs into bg.subgraph_id

	assert(rows(bg.F12_1.counts) == bg.F12_1.num_levels)


	if (verbose>2) printf("\n{txt}# Tricking Bipartite() into dropping singletons instead of doing 1-core prunning\n\n")
	bg.F12_1.counts = J(bg.F12_1.num_levels, 1, 2)

	F = Individual_Factor()
	F.num_levels = FI.num_levels // Not sure why (Stata bug?) but in some cases this does not get propagated and stays missing
	F.num_obs = FG.num_levels
	F.varlist = F.absvar = FG.varlist
	F.aggregate_function = aggregate_function
	F.verbose = verbose
	F.group_p = order(st_data(sample, group_id), 1)
	F.group_inv_p = invorder(F.group_p)
	swap(F.FG, FG)
	swap(F.FI, FI)
	swap(F.bg, bg)

	return(F)
}


// `Individual_Factor' _indiv_factor(`DataFrame' data,
//                | `Boolean' integers_only,
//                  `Boolean' verbose,
//                  `String' method,
//                  `Boolean' sort_levels,
//                  `Boolean' count_levels,
//                  `Integer' hash_ratio,
//                  `Boolean' save_keys,
//                  `Varlist' vars, 			// hack
//                  `DataCol' touse)		 	// hack
// {
// 	assert(0)
// 	`Factor'			F
// 	`Individual_Factor'			FE_F
// 	//F = _factor(data, integers_only, verbose, method, sort_levels, count_levels, hash_ratio, save_keys, vars, touse)
// 	FE_F = Individual_Factor()
// 	//F.swap(FE_F)
// 	return(FE_F)
// }


// `Individual_Factor' fe_join_factors(`Individual_Factor' F1,
//                       `Individual_Factor' F2, 
//                     | `Boolean' count_levels,
//                       `Boolean' save_keys,
//                       `Boolean' levels_as_keys)
// {
// 	`Factor'			F
// 	`Individual_Factor'			FE_F
// 	F = join_factors(F1, F2, count_levels, save_keys, levels_as_keys)
// 	FE_F = Individual_Factor()
// 	F.swap(FE_F)
// 	return(FE_F)
// }

end
// --------------------------------------------------------------------------
// Helper functions
// --------------------------------------------------------------------------

mata:

`Matrix' inv_vec(`Vector' x, `Integer' num_cols)
{
	// Invert the vec() function
	//return(num_cols == 1 ? x : transposeonly(rowshape(x, num_cols)))
	return(num_cols == 1 ? x : colshape(x, num_cols))
}


`Vector' fast_vec(`Matrix' x)
{
	//return(cols(x) == 1 ? x :  vec(x))
	return(cols(x) == 1 ? x :  colshape(x, 1))
}


`Matrix' givens_rotation(`Real' a, `Real' b, `Real' r)
{
	// See http://www.netlib.org/lapack/lawnspdf/lawn150.pdf for technical details and an improved alternative
	`Real' c, s
	if (b) {
		r = hypot(a, b) // sqrt(a ^ 2 + b ^ 2)
		c = a / r
		s = b / r
	}
	else {
		// If a and b are zero, then r=0 and c=s=.
		// To avoid these issues, we follow the "Stable Calculation" of https://en.wikipedia.org/wiki/Givens_rotation
		r = a
		c = 1
		s = 0
	}
	return(c, -s \ s, c)
}


`Real' hypot(`Real' x, `Real' y)
{
	// Naive algorithm: sqrt(x ^ 2 + y ^ 2)
	
	// This algorithm: corrected unfused from page 11 of https://arxiv.org/pdf/1904.09481.pdf
	
	// See also:
	// - https://en.wikipedia.org/wiki/Hypot
	// - https://arxiv.org/abs/1904.09481
	// - https://walkingrandomly.com/?p=6633

	`Real' ax, ay, scale, h, delta
	ax = abs(x)
	ay = abs(y)

	// Ensure ax >= ay
	if (ax < ay) swap(ax, ay)

	// Operands vary widely (ay is much smaller than ay)
	if (ay <= ax * sqrt(epsilon(0.5))) return(ax)

	// Operands do not vary widely
	scale = epsilon(sqrt(smallestdouble())) // rescaling constant
	if (ax > sqrt(maxdouble()/2)) {
		ax = ax * scale
		ay = ay * scale
		scale = 1 / scale
	}
	else if (ay < sqrt(smallestdouble())) {
		ax = ax / scale
		ay = ay / scale
	}
	else {
		scale = 1
	}
	h = sqrt(ax^2+ay^2)

	// Single branch
	delta = h - ax
	h = h - (delta*(2*(ax-ay)) + (2*delta-ay)*ay + delta^2) / (2*h)
	return(h*scale)
}


`Void' assign(`Vector' x, `Real' x1, `Real' x2)
{
	// "assign(x, a, b)" is equivalent to "a, b = x"
	assert(rows(x) == 2)
	assert_msg(!hasmissing(x), "input x has missing values") // BUGBUG remove once we are done debugging
	x1 = x[1]
	x2 = x[2]
}


`Vector' normalize(`Vector' x, `Vector' weights, `Real' norm, `String' msg)
{
	`Vector' normalized_x
	assert(weights!=.) // If we don't want to pass any weights just set them to "1", not to missing
	norm = vector_norm(x, weights, msg)

	if (norm < epsilon(1)) {
		norm = 0
		return(x)
	}
	else {
		return(x / norm)
	}
}


`Real' vector_norm(`vector' x, `Vector' weights, `String' msg)
{
	// Compute 2-norm (much faster than calling -norm-)

	`Real'		ans
	if (weights != 1) assert_msg(rows(x) == rows(weights), "weights non-conforming size")
	ans = sqrt(quadcross(x, weights, x))
	assert_msg(!missing(ans), msg) // TODO: remove in production
	return(ans)
	
	// Note: matrix_norm() computed as not used; formula for Frobenius norm was:
	// sqrt(sum(x :^ 2)) == sqrt(trace(cross(x, x)))
}


`Variables' panelmean(`Variables' y, `FE_Factor' f)
{
	assert(0) // replaced by f.panelmean() which feels more natural
	assert_boolean(f.has_weights)
	assert_msg(f.panel_is_setup, "F.panel_setup() must be run beforehand")
	if (f.has_weights) {
		return(editmissing(f.panelsum(y) :/ f.weighted_counts, 0))
	}
	else {
		return(f.panelsum(y) :/ f.counts)
	}
}


`Matrix' precompute_inv_xx(`FE_Factor' f)
{
	`Integer'               i, L, K, offset
	`Vector'             	tmp_x, tmp_w
	`Matrix'                inv_xx, tmp_inv_xx, sorted_x
	`RowVector'             tmp_xmeans

	// note: f.x and f.weights must be already sorted by the factor f
	assert_boolean(f.has_weights)
	assert_boolean(f.has_intercept)

	L = f.num_levels
	K = cols(f.unsorted_x)
	inv_xx = J(L * K, K, .)
	sorted_x = f.sort(f.unsorted_x)

	for (i=1; i<=L; i++) {
		tmp_x = panelsubmatrix(sorted_x, i, f.info)
		tmp_w = f.has_weights ? panelsubmatrix(f.weights, i, f.info) : 1
		if (f.has_intercept) {
			tmp_xmeans = K > 1 ? f.x_means[i, .] : f.x_means[i]
			tmp_inv_xx = invsym(quadcrossdev(tmp_x, tmp_xmeans, tmp_w, tmp_x, tmp_xmeans))
		}
		else {
			tmp_inv_xx = invsym(quadcross(tmp_x, tmp_w, tmp_x))
		}
		offset = K * (i - 1)
		inv_xx[|offset + 1, 1 \ offset + K , . |] = tmp_inv_xx
	}
	return(inv_xx)
}


`Variables' panelsolve_invsym(`Variables' y, `FE_Factor' f, | `Matrix' alphas)
{
	`Integer'               i, L, K, offset
	`Boolean'               save_alphas
	`Matrix'             	xbd, tmp_xbd, tmp_x, tmp_y, tmp_xy, tmp_inv_xx, sorted_x
	`Vector'             	tmp_w, b
	`RowVector'             tmp_xmeans, tmp_ymeans

	// note: y, f.x, and f.weights must be already sorted by the factor f
	assert_boolean(f.has_weights)
	assert_boolean(f.has_intercept)
	assert_msg(rows(f.unsorted_x) == rows(y))

	save_alphas = args()>=3 & alphas!=J(0,0,.)
	if (save_alphas) assert(cols(y)==1)

	L = f.num_levels
	K = cols(f.unsorted_x)
	xbd = J(rows(y), cols(y), .)
	sorted_x = f.sort(f.unsorted_x)

	for (i=1; i<=L; i++) {
		tmp_y = panelsubmatrix(y, i, f.info)
		tmp_x = panelsubmatrix(sorted_x, i, f.info)
		tmp_w = f.has_weights ? panelsubmatrix(f.weights, i, f.info) : 1
		offset = K * (i - 1)
		tmp_inv_xx = f.inv_xx[|offset + 1, 1 \ offset + K , . |]

		if (f.has_intercept) {
			tmp_ymeans = mean(tmp_y, tmp_w)
			tmp_xmeans = K > 1 ? f.x_means[i, .] : f.x_means[i]
			tmp_xy = quadcrossdev(tmp_x, tmp_xmeans, tmp_w, tmp_y, tmp_ymeans)
			b = tmp_inv_xx * tmp_xy
			tmp_xbd = (tmp_x :- tmp_xmeans) * b :+ tmp_ymeans
			if (save_alphas) alphas[i, .] = tmp_ymeans - tmp_xmeans * b, b'
		}
		else {
			tmp_xy = quadcross(tmp_x, tmp_w, tmp_y)
			b = tmp_inv_xx * tmp_xy
			tmp_xbd = tmp_x * b
			if (save_alphas) alphas[i, .] = b'
		}
		xbd[|f.info[i,1], 1 \ f.info[i,2], .|] = tmp_xbd
	}
	return(f.invsort(xbd))
}


`RowVector' function compute_stdevs(`Matrix' A)
{
	`RowVector'				stdevs
	`Integer'				K, N
	// Note: each col of A will have stdev of 1 unless stdev is quite close to 0

	// We don't need good accuracy for the stdevs, so we have a few alternatives:
	// Note: cross(1,A) is the same as colsum(A), but faster
	// Note: cross(A, A) is very fast, but we only need the main diagonals
	// [A: 1sec] stdevs = sqrt( (colsum(A:*A) - (cross(1, A) :^ 2 / N)) / (N-1) )
	// [B: .61s] stdevs = sqrt( (diagonal(cross(A, A))' - (cross(1, A) :^ 2 / N)) / (N-1) )
	// [C: .80s] stdevs = diagonal(sqrt(variance(A)))'
	// [D: .67s] means = cross(1, A) / N; stdevs =  sqrt(diagonal(crossdev(A, means, A, means))' / (N-1))

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

	// (D) Equilibrate matrix columns instead of standardize (i.e. just divide by column max)
	// _perhapsequilc(A, stdevs=.)
	// stdevs = 1 :/ stdevs
	// assert_msg(!missing(stdevs), "stdevs are missing; is N==1?")

	// (E) Don't do anything
	// stdevs = J(1, cols(A), 1)

	return(stdevs)
}


`RowVector' function standardize_data(`Matrix' A)
{
	`RowVector'				stdevs
	assert_msg(!isfleeting(A), "input cannot be fleeting")
	stdevs = compute_stdevs(A)
	A = A :/ stdevs
	return(stdevs)
}


`Variables' st_data_pool(`Vector'  sample,		// observations to load
					     `Varlist' vars,		// variables to load
					     `Integer' poolsize,	// how many variables to load each time
					   | `Boolean' compact,		// whether to trim the dataset to save memory or not
					     `Varlist' keepvars,	// what extra vars to keep if we trim the dataset (clustervars, timevar, panelvar)
					     `Boolean' verbose)
{
	`Integer'               i, j, k
	`Varlist'				temp_keepvars
	`Variables'				data

	if (args()<4 | compact == .) compact = 0
	if (args()<6 | verbose == .) verbose = 0

	k = cols(vars)
	assert_msg(poolsize > 0, "poolsize must be a positive integer")

	if (k <= poolsize) {
		data = st_data(sample, vars)
		if (compact) compact_dataset(keepvars)
	}
	else {
		data = J(rows(sample), 0, .)
		for (i=1; i<=k; i=i+poolsize) {
			j = i + poolsize - 1
			if (j>k) j = k
			data = data, st_data(sample, vars[i..j])
			if (compact) {
				temp_keepvars = j == k ? keepvars :  keepvars , vars[j+1..k]
				compact_dataset(temp_keepvars)
			}
		}
	}

	assert_msg(k == cols(data), "could not load data into Mata correctly; please contact author")
	return(data)
}


`Void' compact_dataset(`Varlist' keepvars)
{
	keepvars = tokens(invtokens(keepvars)) // deal with empty strings
	if (cols(keepvars)) {
		stata(sprintf("fvrevar %s, list", invtokens(keepvars)))
		stata(sprintf("novarabbrev keep %s", st_global("r(varlist)")))
	}
	else {
		stata("clear")
	}
}


`Void' add_undocumented_options(`String' object, `String' options, `Integer' verbose)
{
	`StringRowVector'		tokens
	`Integer'				num_tokens, i
	`String'				key, val, cmd, msg

	if (options == "") return
	tokens = tokens(options, " ()")
	
	msg = sprintf("options {bf:%s} not allowed", options)
	assert_msg(!mod(cols(tokens), 4), msg, 198, 0) // 198: option not allowed; 0:traceback off

	num_tokens = trunc(cols(tokens) / 4)
	for (i=0; i<num_tokens; i++) {
		key = tokens[4*i+1]
		assert_msg(tokens[4*i+2] == "(")
		val = tokens[4*i+3]
		assert_msg(tokens[4*i+4] == ")")
		if (verbose > 1) printf("{txt}    %s.%s = %s\n", object, key, val)
		cmd = sprintf("cap mata: %s.%s = %s", object, key, val)
		stata(cmd)
		assert_msg(!st_numscalar("c(rc)"), sprintf("option {bf:%s} not allowed", key), 198, 0)
	}
}


// This short-but-difficult function receives a list of dropped rows (idx) at the group level
// And returns the list of dropped rows (indiv_idx) at the individual level
`Vector' get_indiv_idx(`Vector' sample, `Vector' indiv_sample, `Varname' group_id, `Vector' idx)
{
	`Factor'				F
	`Vector'				mask, dropped_levels, indiv_idx

	// varnames, touse, verbose, method, SORT_LEVELS, count_levels, hash_ratio, save_keys
	F = factor(group_id, indiv_sample, ., "", 1, 1, ., 0) // not completely sure that we need sort_levels==1
	
	// Map observations to group identifiers
	mask = J(st_nobs(), 1, 0)
	mask[indiv_sample] = F.levels

	// List what groups we are dropping
	dropped_levels = mask[sample[idx]]

	// Create mask indicating whether we drop each group or not
	mask = create_mask(F.num_levels, 0, dropped_levels, 1)

	// Select rows within indiv_sample that correspond to the groups we are dropping
	indiv_idx = selectindex(mask[F.levels])

	return(indiv_idx)
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

end
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
// --------------------------------------------------------------------------
// FixedEffects main class
// --------------------------------------------------------------------------

mata:

class FixedEffects
{
	// TODO: rename N to num_obs and M to num_levels (or sum_num_levels?)

	// Factors
	`Integer'				G                   // # FEs
	`Integer'				N                   // # obs
	`Integer'				M                   // # levels of FEs (sum over all FEs)
	`FE_Factors'			factors				// vector of Factor() objects
	`Vector'				sample				// regression sample (index vector)
	`String'				tousevar			// indicator variable that matches the regression sample
	`Integer'				num_singletons		// Number of detected singletons
	`Varlist'               absvars
	`Varlist'               extended_absvars	// used to label variables (estimates of the FEs)

	// Specific to individual FEs
	`Varname'				group_id
	`Varname'				individual_id		// We have an individual FE if (individual_id != "")
	`Varname'				indiv_tousevar
	`String'				function_individual
	`Vector'				indiv_sample
	`Variable'				indiv_weights		// unsorted weight
	`Boolean'				has_indiv_fe
	
	// Weight-specific
	`Boolean'				has_weights
	`Variable'				weights             // unsorted weight
	`String'				weight_var          // Weighting variable
	`String'				weight_type         // Weight type (pw, fw, etc)
	`Variable'				true_weights        // custom case for IRLS (ppmlhdfe) where weights = mu * true_weights
	`Real'					weights_scale		// if there are very low weights (e.g. ppmlhdfe) we will rescale weights to ensure no weights are below a given epsilon; afterwards multiply by this to obtain the correct sum of weights

	// Saving FEs
	`Boolean'				save_any_fe
	`Boolean'				save_all_fe
	`Boolean'				storing_alphas		// set to 1 after coef estimation, where we are backing-out alphas

	// Constant-related
	`Boolean'				has_intercept		// 1 if at least one FE has intercepts (i.e. not slope-only)

	// Multiprocessing
	`Integer'				parallel_maxproc	// Up to how many parallel processes will we launch (0=disabled)
	`Boolean'				parallel_force		// Force the use of parallel even for one worker
	`String'				parallel_opts		// Options passed to parallel_map
	`String'				parallel_dir		// Path where we will store temporary objects passed to the workers

	// Misc
	`Boolean'				initialized			// 0 until we run .init()
	`Integer'				verbose
	`Boolean'				timeit

	// Parallel
	`Integer'				parallel_numproc	// How many worker processes are we running
	`RowVector'				parallel_poolsizes	// Pool size of each block


	// Solver settings
	`String'				technique			//
	`Integer'				maxiter
	`Integer'				miniter
	`Real'					tolerance
	`String'				transform           // Kaczmarz Cimmino Symmetric_kaczmarz (k c s)
	`String'				acceleration        // Acceleration method. None/No/Empty is none\
	`String'				preconditioner
	`Boolean'               abort               // Raise error if convergence failed?

	// MAP Solver settings
	`Integer'				min_ok				// minimum number of "converged" iterations until we declare convergence
	`Integer'               accel_start         // Iteration where we start to accelerate // set it at 6? 2?3?
	`Integer'               accel_freq          // Specific to Aitken's acceleration

	// Memory usage
	`Integer'				poolsize
	`Boolean'				compact				//

	// Solver results
	`Solution'				solution			// Stata does not allow "s = HDFE.partial_out()" so we have to keep it here

	
	//`Real'				finite_condition
	//`Real'				compute_rre         // Relative residual error: || e_k - e || / || e ||
	//`Real'				rre_depvar_norm
	//`Vector'				rre_varname
	//`Vector'				rre_true_residual
	`String'				panelvar
	`String'				timevar

	// Absorbed degrees-of-freedom computations
	`StringRowVector'		dofadjustments
	`Integer'				G_extended          // Number of intercepts plus slopes
	`Vector'				doflist_K			// For each extended absvar, number of coefs
	`Vector'				doflist_M			// For each extended absvar, number of redundant coefs
	`Vector'				doflist_M_is_exact  // For each extended absvar, whether redundant coef computation is exact
	`Vector'				doflist_M_is_nested // For each extended absvar, whether coefs are redundant due to cluster nesting

	`Integer'				df_a_nested 		// Redundant coefs due to being nested; used for: r2_a r2_a_within rmse
	`Integer'				df_a_redundant      // e(mobility)
	`Integer'				df_a_initial		// df_a assuming there are no redundant coefs at all
	`Integer'				df_a                // df_a_inital - df_a_redundant

	// VCE and cluster variables (we use clustervars when computing DoF so this belongs to HDFE object instead of solution object)
	`String'				vcetype
	`Integer'				num_clusters
	`Varlist'				clustervars
	`Varlist'				base_clustervars

	`Boolean'				drop_singletons
	`String'				absorb              // contents of absorb()

	// Methods
	`Void'					new()
	`Void'					init()

	`Void'					partial_out()
	`Void'					save_touse()
	`Void'					save_variable()
	`Void'					post_footnote()
	`Void'					store_alphas()
	`Void'					store_alphas_map()
	`Void'					store_alphas_lsmr()

	// Matrix methods; used by LSQR/LSMR
	`Integer'				num_rows()
	`Integer'				num_cols()
	`Vector'				mult()
	`Vector'				mult_transpose()

	// Methods used to load and reload weights (and update cvars, etc. when weights change)
	`Void'					load_weights() // calls update_preconditioner, etc.
	`Void'					update_preconditioner()

	// Methods used by alterating projections (MAP)
	`Variables'				project_one_fe()
}


`Void' FixedEffects::new()
{
	// Variables not yet initialized
	initialized = 0
	absvars = tousevar = ""
	num_singletons = 0
	sample = indiv_sample = J(0, 1, .)

	// Default values
	technique = "lsmr"
	preconditioner = "block_diagonal"
	drop_singletons = 1
	weight_type = weight_var = ""
	weights = indiv_weights = 1 // set to 1 so cross(x, S.weights, y) == cross(x, y)
	true_weights = J(0, 1, .)
	weights_scale = . // must be a real after initializing
	abort = 1
	verbose = 0
	timeit = 0
	has_indiv_fe = 0

	// MAP Solver settings
	min_ok = 1
	accel_start = 6
	accel_freq = 3

	// Individual FEs
	group_id = ""
	individual_id = ""
	indiv_tousevar = ""

	// Multiprocessing
	parallel_maxproc = 0
	parallel_force = 0
	parallel_opts = "" // commented out else fgetmatrix() crashes
	parallel_dir = "" // commented out else fgetmatrix() crashes

	storing_alphas = 0
}


`Void' FixedEffects::init()
{
	`StringRowVector'		ivars, cvars, targets
	`RowVector'				intercepts, num_slopes
	`Vector'                idx, indiv_idx
	`Integer'				num_singletons_i
	`Integer'				i, g, gg, j, remaining, rc
	`Boolean'				already_found_individual
	`Boolean'				is_fw_or_iw, zero_threshold

	// To initialize, we use the following attributes:
	// 	.absvars			: contents of absorb()
	// 	.tousevar			: variable indicating the sample
	// 	.weight_type		: fweight, aweight, pweight, iweight
	// 	.weight_var			: weighting variable
	// 	.drop_singletons	: drop singletons if =1
	// 	.preconditioner		: type of preconditioner used
	// 	.verbose			: how much extra info to show (-1, 0, 1, 2, 3)

	assert_msg(initialized == 0, "HDFE already initialized")

	has_indiv_fe = this.individual_id != ""

	// Construct touse if needed
	if (this.tousevar == "") {
		this.tousevar = st_tempname()
		st_store(., st_addvar("byte", this.tousevar, 1), J(st_nobs(), 1, 1))
	}

	// Parse contents of absorb()
	rc = _stata(`"ms_parse_absvars "' + this.absvars)
	assert_msg(!rc, "error in option absorb()", rc, 0)
	assert_msg(st_global("s(options)") == "", sprintf("invalid suboptions in absorb(): %s", st_global("s(options)")), 198, 0)
	if (verbose > 0) stata(`"sreturn list"')

	// Update sample with variables used in the fixed effects
	stata(sprintf("markout %s %s, strok", this.tousevar, st_global("s(basevars)")))
	if (has_indiv_fe) {
		stata(sprintf("markout %s %s, strok", this.indiv_tousevar, st_global("s(basevars)")))
	}

	// Tokenize strings
	this.absvars = tokens(st_global("s(absvars)"))
	this.has_intercept = strtoreal(st_global("s(has_intercept)"))
	ivars = tokens(st_global("s(ivars)"))
	cvars = tokens(st_global("s(cvars)"))
	targets = strtrim(tokens(st_global("s(targets)")))
	intercepts = strtoreal(tokens(st_global("s(intercepts)")))
	num_slopes = strtoreal(tokens(st_global("s(num_slopes)")))

	// Quick sanity check
	this.G = strtoreal(st_global("s(G)"))
	assert(this.G == cols(this.absvars))

	// Fill out object
	this.factors = FE_Factor(G)
	this.extended_absvars = tokens(st_global("s(extended_absvars)"))
	this.save_any_fe = strtoreal(st_global("s(save_any_fe)"))
	this.save_all_fe = strtoreal(st_global("s(save_all_fe)"))

	// Load sample and weight variables
	this.sample = selectindex(st_data(., tousevar))
	if (this.weight_var != "") this.weights = st_data(this.sample, this.weight_var) // just to use it in f.drop_singletons()

	if (has_indiv_fe) {
		assert(this.weight_type != "fweight") // because there are no duplicated obs, it makes no sense to have fw and indiv FEs with groups
		this.indiv_sample = selectindex(st_data(., indiv_tousevar))
		if (this.weight_var != "") this.indiv_weights = st_data(this.indiv_sample, this.weight_var) // just to use it in f.drop_singletons() // BUGBUG do I really need/use this?
		this.group_id = group_id
		this.individual_id = individual_id
	}

	if (verbose > 0) printf("\n{txt}{title:Loading fixed effects information:}\n")

	if (G == 1 & ivars[1] == "_cons") {
		// Special case without any fixed effects
		this.N = rows(this.sample)
		this.factors[1].num_obs = this.N
		this.factors[1].counts = this.N
		this.factors[1].num_levels = 1
		this.factors[1].levels = J(this.N, 1, 1)
		this.factors[1].is_sorted = 1
		this.factors[1].method = "none"
		this.factors[1].is_individual_fe = 0
		if (verbose > 0) {
			logging_singletons("start", absvars)
			logging_singletons("iter1", absvars, this.factors[1], intercepts[1], num_slopes[1], rows(this.sample), i, 1)
			logging_singletons("iter2", absvars, this.factors[1])
			printf(" %10.0g   {c |}", 0)
			printf("\n")
			logging_singletons("end", absvars)
		}
	}
	else {
		// (1) create the factors and remove singletons
		if (verbose > 0) logging_singletons("start", absvars)
		remaining = this.G
		i = 0
		already_found_individual = 0
		while (remaining) {
			++i
			g = 1 + mod(i-1, this.G)

			if (verbose > 0) logging_singletons("iter1", absvars, this.factors[g], intercepts[g], num_slopes[g], rows(this.sample), i, g)

			if (rows(this.sample) < 2) {
				if (verbose > 0) printf("\n")
				exit(error(2001))
			}

			if (i <= this.G) {
				if (ivars[g] == this.individual_id) {
					// Careful that Mata has a few weird bugs around class inheritance and virtual methods
					already_found_individual = 1
					//this.factors[g] = Individual_Factor()
					this.factors[g] = indiv_factor(group_id, individual_id, this.sample, this.indiv_sample, function_individual, verbose)
					assert_msg(!missing(factors[g].num_levels), "Stata/Mata bug: Individual_Factor.num_levels is missing")
					assert(this.factors[g].is_individual_fe == 1)
				}
				else {
					// We don't need to save keys (or sort levels but that might change estimates of FEs)
					this.factors[g] = fe_factor(ivars[g], this.sample, ., "", ., 1, ., 0)
					assert(this.factors[g].is_individual_fe == 0)
				}
			}

			if (verbose > 0) logging_singletons("iter2", absvars, this.factors[g])

			if (drop_singletons) {
				is_fw_or_iw = (this.weight_type == "fweight") | (this.weight_type == "iweight")
				zero_threshold = (this.weight_type == "iweight") ? 1 : 0 // used by ppmlhdfe
				if (this.factors[g].is_individual_fe) {
					idx = factors[g].drop_singletons(this.sample, this.indiv_sample)
				}
				else {
					idx = this.factors[g].drop_singletons(is_fw_or_iw ? this.weights : J(0, 1, .), zero_threshold)
				}
				num_singletons_i = rows(idx)
				this.factors[g].num_singletons = this.factors[g].num_singletons + num_singletons_i
				this.num_singletons = this.num_singletons + num_singletons_i
				if (verbose > 0) {
					printf(" %10.0g   {c |}", num_singletons_i)
					displayflush()
				}

				if (num_singletons_i == 0) {
					--remaining
				}
				else {
					remaining = this.G - 1

					// BUGBUG we don't do this if the FE we just updated is the indiv one
					if (has_indiv_fe) {
						indiv_idx = get_indiv_idx(this.sample, this.indiv_sample, this.group_id, idx)
						this.indiv_sample[indiv_idx] = J(rows(indiv_idx), 1, 0)
						this.indiv_sample = select(this.indiv_sample, this.indiv_sample)
					}
					
					// sample[idx] = . // not allowed in Mata; instead, make 0 and then select()
					this.sample[idx] = J(rows(idx), 1, 0)

					// factor.drop_singletons() above will trim this.weights for fw and iw but not pw or aw, so we need to do it here
					if (this.weight_type == "aweight" | this.weight_type == "pweight") this.weights = select(this.weights, this.sample)

					this.sample = select(this.sample, this.sample)

					for (j=i-1; j>=max((1, i-remaining)); j--) {
						gg = 1 + mod(j-1, this.G)
						this.factors[gg].drop_obs(this.factors[gg].is_individual_fe ? indiv_idx : idx)
						if (verbose > 0) printf("{res} .")
					}
				}
			}
			else {
				if (verbose > 0) printf("      n/a     {c |}")
				--remaining
			}
			if (verbose > 0) printf("\n")
		}
		if (verbose > 0) logging_singletons("end", absvars)
	}

	// Save updated e(sample); needed b/c singletons reduce sample; required to parse factor variables that we want to partial out
	if (drop_singletons & num_singletons) this.save_touse()

	if (has_indiv_fe) {
		assert_msg(already_found_individual, "individual id not found on absorb()", 100, 0)
	}
	
	if (drop_singletons & this.num_singletons>0 & verbose>-1 | this.factors[1].num_obs<2) {
		if (this.weight_type=="iweight") {
			// PPML-specific
			printf(`"{txt}(dropped %s observations that are either {browse "http://scorreia.com/research/singletons.pdf":singletons} or {browse "http://scorreia.com/research/separation.pdf":separated} by a fixed effect)\n"', strofreal(this.num_singletons))
		}
		else {
			printf(`"{txt}(dropped %s {browse "http://scorreia.com/research/singletons.pdf":singleton observations})\n"', strofreal(this.num_singletons))
		}
	}

	if (this.factors[1].num_obs < 2) {
		exit(error(2001))
	}

	this.N = factors[1].num_obs // store number of obs.
	//assert_msg(this.N == factors[this.G].num_obs, "N check #1")
	assert_msg(this.N > 1, "N check #2")
	assert_msg(!missing(this.N), "N check #3")

	// (2) run *.panelsetup() after the sample is defined
	// (3) load cvars
	if (verbose > 0) printf("\n{txt}# Initializing panelsetup() and loading slope variables for each FE\n\n")
	for (g=1; g<=this.G; g++) {
		assert_in(this.factors[g].is_individual_fe, (0, 1))

		this.factors[g].absvar = absvars[g]
		this.factors[g].ivars = tokens(ivars[g])
		this.factors[g].cvars = tokens(cvars[g])
		this.factors[g].has_intercept = intercepts[g]
		this.factors[g].num_slopes = num_slopes[g]
		this.factors[g].target = targets[g]
		this.factors[g].save_fe = (targets[g] != "")

		if (verbose > 0) printf("{txt}   - Fixed effects: {res}%s{txt}\n", this.factors[g].absvar)
		if (verbose > 0) printf("{txt}     panelsetup()\n")

		this.factors[g].panelsetup()

		if (this.factors[g].num_slopes) {
			// Load, standardize, sort by factor, and store
			// Don't precompute aux objects (xmeans, inv_xx) as they depend on the weights
			// and will be computed on step (5)
			if (verbose > 0) printf("{txt}     loading cvars({res}%s{txt})\n", strtrim(invtokens(cvars)))
			//this.factors[g].x = this.factors[g].sort(st_data(this.sample, this.factors[g].cvars))
			if (this.factors[g].is_individual_fe) {
				this.factors[g].unsorted_x = st_data(this.indiv_sample, this.factors[g].cvars)
			}
			else {
				this.factors[g].unsorted_x = st_data(this.sample, this.factors[g].cvars)
			}
			this.factors[g].x_stdevs = standardize_data(this.factors[g].unsorted_x)
		}
	}

	// (5) load weight
	this.load_weights(1) // update this.has_weights, this.factors, etc.

	// (6) Update sort order of indiv FEs, if needed
	if (has_indiv_fe) {
		for (g=1; g<=this.G; g++) {
			if (this.factors[g].is_individual_fe==0) {
				this.factors[g].group_p = order(st_data(this.sample, this.group_id), 1)
				this.factors[g].group_inv_p = invorder(this.factors[g].group_p)
			}
		}
	}


	// Mark the FixedEffects() object as correctly initialized
	this.initialized = 1
}


`Integer' FixedEffects::num_rows()
{
	assert(!missing(N))
	return(N)
}


`Integer' FixedEffects::num_cols()
{	
	assert_msg(!missing(M))
	return(M)
}



`Vector' FixedEffects::mult(`Vector' x)
{
	// Return A*x (given a vector -x- of alphas, this is the sum of the alphas for each obs.)
	`Integer'				g, i, n
	`Vector'				ans

	assert_msg(M == rows(x), "Matrix A and vector x not conformable")
	`DEBUG' assert_msg(!hasmissing(x), ".mult() received missing values")
	assert_msg(cols(x) == 1)
	ans = J(this.N, 1, 0) // sum starts as empty, then we add to it for each FE

	for (g=i=1; g<=G; g++) {
		n = factors[g].num_levels * (factors[g].has_intercept + factors[g].num_slopes)
		factors[g].mult(x[|i \ i + n -1|], preconditioner, ans)
		i = i + n
	}
	`DEBUG' assert_msg(!hasmissing(ans), ".mult() returned missing values")
	`DEBUG' assert_msg(cols(ans) == 1)
	return(ans)
}


`Vector' FixedEffects::mult_transpose(`Vector' y)
{
	// Return A'x (given a variable -x-, this is the sum of -x- within each group)
	// BUGBUG TODO: is there a preconditioner that makes the intercept+slope case equivalent to inv(x'x)x'y for cols(x)=2 ???

	`Integer'				g, i, n
	`Vector' 				alphas

	assert_msg(this.N == rows(y), "Matrix A and vector y not conformable")
	assert_msg(cols(y) == 1)
	alphas = J(this.M, 1, .)

	for (g=i=1; g<=this.G; g++) {
		n = factors[g].num_levels * (factors[g].has_intercept + factors[g].num_slopes)
		alphas[| i \ i+n-1 |] = factors[g].mult_transpose(y, preconditioner)
		i = i + n
	}
	`DEBUG' assert_msg(cols(alphas) == 1)
	return(alphas)
}


// This adds/removes weights or changes their type
`Void' FixedEffects::load_weights(`Boolean' verbose)
{
	`Integer'				g

	this.has_weights = (this.weight_type != "" & this.weight_var != "")
	if (this.verbose > 0 & verbose > 0 & this.has_weights) printf("\n{txt}# Loading weights [{res}%s{txt}={res}%s{txt}]\n\n", this.weight_type, this.weight_var)

	// If needed, load weight from dataset (we'll need to update weights if we later find singletons)
	// If there are singletons that we drop, we need to ensure we reload/update weights

	if (this.has_weights) {
		if (this.weights == 1) this.weights = st_data(this.sample, this.weight_var)
		assert_msg(rows(this.weights) == rows(this.sample))
		assert_msg(!hasmissing(this.weights), "weights can't be missing")
	}
	else {
		this.weights = 1
	}

	// Rescale weights so there are no weights below 1, to avoid numerical issues
	// This is not always possible b/c if weights_scale is to low then we get other types of problems
	// We pick 1.0X-017 ~ 1e7 because it's the precision of floats
	// see: https://journals.sagepub.com/doi/pdf/10.1177/1536867X0800800207
	weights_scale = 1
	if (this.weight_type != "fweight") {
		weights_scale = clip(colmin(this.weights), 1.0X-017, 1)
		if (weights_scale != 1) this.weights = this.weights :/ weights_scale
	}

	M = 0 // Count number of fixed effects

	for (g=1; g<=this.G; g++) {
		M = M + factors[g].num_levels * (factors[g].has_intercept + factors[g].num_slopes)

		if (this.has_weights) {
			factors[g].set_weights(factors[g].is_individual_fe ? this.indiv_weights : this.weights, this.verbose > 0 & verbose > 0)
		}
	}

	this.update_preconditioner()
}


`Void' FixedEffects::update_preconditioner()
{
	`FE_Factor'				F
	`Integer'				g

	if (this.weight_type == "iweight") {
		// (this is meaningless with iweights) --> why?
		// pass
	}
	else if (this.technique == "map") {

		for (g=1; g<=G; g++) {
			assert(factors[g].panel_is_setup==1) // bugbug remove
			// Update mean(z; w) and inv(z'z; w) where z is a slope variable and w is the weight
			if (factors[g].num_slopes) {
				if (verbose > 0) printf("{txt}   - preprocessing cvars of factor {res}%s{txt}\n", factors[g].absvar)
				if (factors[g].has_intercept) {
					factors[g].x_means = factors[g].panelmean(factors[g].unsorted_x, 1)
				}
				factors[g].inv_xx = precompute_inv_xx(factors[g])
			}
		}
	}
	else if (preconditioner == "none") {
		// pass
	}
	else if (preconditioner == "diagonal") {
		for (g=1; g<=G; g++) {
			factors[g].init_diagonal_preconditioner()
		}
	}
	else if (preconditioner == "block_diagonal") {
		for (g=1; g<=this.G; g++) {
			// Use the simpler diagonal preconditioner when we have only one coef. per FE
			if (factors[g].has_intercept + factors[g].num_slopes > 1) {
				factors[g].init_block_diag_preconditioner(this.weight_type, this.weights_scale, true_weights)
			}
			else {
				factors[g].init_diagonal_preconditioner()
			}
		}
	}
	else {
		assert(0)
	}
}


`Void' FixedEffects::save_touse(| `Varname' touse, `Boolean' replace)
{
	`Integer'				idx
	`Vector'				mask

	// Set default arguments
	if (args()<1 | touse=="") {
		assert(tousevar != "")
		touse = tousevar
	}

	// Note that args()==0 implies replace==1 (else how would you find the name)
	if (args()==0) replace = 1
	if (args()==1 | replace==.) replace = 0

	if (verbose > 0) printf("\n{txt}# %s e(sample)\n", replace ? "Updating" : "Saving")

	// Compute dummy vector
	mask = create_mask(st_nobs(), 0, sample, 1)

	// Save vector as variable
	if (replace) {
		st_store(., touse, mask)
	}
	else {
		idx = st_addvar("byte", touse, 1)
		st_store(., idx, mask)
	}


	// Do the same for indiv touse
	
	if (has_indiv_fe) {
		// Compute dummy vector
		mask = create_mask(st_nobs(), 0, indiv_sample, 1)

		// Save vector as variable
		if (replace) {
			st_store(., indiv_tousevar, mask)
		}
		else {
			idx = st_addvar("byte", indiv_sample, 1)
			st_store(., idx, mask)
		}
	}
}


`Void' FixedEffects::partial_out(`Varlist' varnames,
								|`Boolean' save_tss,
								 `Boolean' standardize_inputs)
{
	// Features:
	// 1) Allow both varlists and matrices as inputs
	// 2) solution.check(depvar_pos) is better than first_is_depvar
	// 3) can compute TSS of whatever column we want (and save it in solution.tss[i])

	// Notes:
	// 1) don't need standardize_inputs==2 as we will have solution.undo_standardize()
	// 2) first_is_depvar is problematic with cache... 
	// 3) Let's just save tss of all variables, it's easier with cache and simpler overall

	if (this.timeit) timer_on(40)

	`Boolean'				input_is_varlist
	`Boolean'				solver_failed
	`Integer'				i, g
	`Solution'				temp_sol
	`Matrix'				data
	`String'				msg
	`Matrix'				tss
	`Varlist'				keepvars

	// Used by MAP:
	`FunctionP'				fun_transform
	`FunctionP'				fun_accel

	assert_msg(initialized == 1, "must run HDFE.init() first")

	// Defaults
	if (args()<2 | save_tss==.) save_tss = 1
	if (args()<3 | standardize_inputs==.) standardize_inputs = 1

	solution = Solution() // clean up any previous results
	solution.varlist = tokens(invtokens(varnames)) // trick to allow string scalars and string vectors
	solution.K = cols(solution.varlist)
	// TODO: when implementing _partial.out, set .varlist ="variable #" :+ strofreal(1..cols(data))

	if (verbose > 0) {
		msg = poolsize > solution.K ? " in a single block" : sprintf(" in blocks of %g", poolsize)
		printf("\n{txt}# Loading and partialling out %g variables%s\n\n", solution.K, msg)
	}

	// Sanity checks
	assert_in(technique, ("map", "lsmr", "lsqr", "gt"), sprintf(`"invalid algorithm "%s""', technique))
	assert_in(preconditioner, ("none", "diagonal", "block_diagonal"), sprintf(`"invalid LSMR/LSQR preconditioner "%s""', preconditioner))
	if (G==1) acceleration = "none" // Shortcut for trivial case (1 FE)
	if (transform=="kaczmarz" & acceleration=="conjugate_gradient") printf("{err}(WARNING: convergence is {bf:unlikely} with transform=kaczmarz and accel=CG)\n")

	// This is from cleanup_for_parallel() in the Parallel .mata file
	// However, it helps to save space (specially for indiv FEs where the bipartite graph might take space)
	// So we delete this now, as good practice
	for (g=1; g<=this.G; g++) {
		this.factors[g].cleanup_before_saving()
	}

	// Load data
	if (verbose > 0) printf("{txt}   - Loading data into Mata: {res}%s{txt}\n", invtokens(varnames))
	if (compact) keepvars = this.base_clustervars , this.timevar, this.panelvar
	if (this.timeit) timer_on(41)
	data = st_data_pool(sample, solution.varlist, poolsize, compact, keepvars, verbose)
	if (this.timeit) timer_off(41)

	// Sanity check
	assert(weights == 1 | rows(weights) == rows(data))

	if (this.timeit) timer_on(42)

	// Compute 2-norm of each var, to see if we need to drop as regressors (BUGBUG: use weights?)
	solution.norm2 = diagonal(cross(data, data))'

	// Compute and save means of each var
	solution.means = mean(data, weights)
	

	// Compute TSS
	if (save_tss) {
		if (verbose > 0) printf("{txt}   - Computing total sum of squares\n")
		if (has_intercept) {
			tss = crossdev(data, solution.means, weights, data, solution.means) // Sum of w[i] * (y[i]-y_mean) ^ 2
		}
		else {
			tss = cross(data, weights, data) // Sum of w[i] * y[i] ^ 2
		}
		tss = diagonal(tss)' // Go from K*K to 1*K
		if (weight_type=="aweight" | weight_type=="pweight") tss = tss :* (rows(data) / sum(weights))
		solution.tss = tss
	}

	
	// Standardize variables
	if (standardize_inputs) {
		if (verbose > 0) printf("{txt}   - Standardizing variables\n")
		solution.stdevs = standardize_data(data)
		solution.norm2 = solution.norm2 :/ solution.stdevs :^ 2
		solution.means = solution.means :/ solution.stdevs
	}

	// Attach more info to solution
	solution.is_standardized = standardize_inputs
	solution.has_intercept = has_intercept

	if (this.timeit) timer_off(42)
	
	// Partial out
	if (verbose > 1) printf("\n")

	if (parallel_maxproc > 0) {
		cleanup_for_parallel(this)
		if (this.timeit) timer_on(49)
		save_before_parallel(parallel_dir, this, data) // updates this.parallel_poolsize
		if (this.timeit) timer_off(49)
		if (this.timeit) timer_off(40)
		exit()
	}

	if (technique == "map") {
		// Load transform pointer
		if (transform=="cimmino") fun_transform = &transform_cimmino()
		if (transform=="kaczmarz") fun_transform = &transform_kaczmarz()
		if (transform=="symmetric_kaczmarz") fun_transform = &transform_sym_kaczmarz()
		if (transform=="random_kaczmarz") fun_transform = &transform_rand_kaczmarz() // experimental
		if (transform=="unrolled_sym_kaczmarz") fun_transform = &transform_unrolled_sym_kaczmarz() // experimental
		if (transform=="commutative_kaczmarz") fun_transform = &transform_commutative_kaczmarz() // experimental
		assert_msg(fun_transform != NULL, "NULL transform function")

		// Pointer to acceleration routine
		if (acceleration=="test") fun_accel = &accelerate_test()
		if (acceleration=="none") fun_accel = &accelerate_none()
		if (acceleration=="conjugate_gradient") fun_accel = &accelerate_cg()
		if (acceleration=="steepest_descent") fun_accel = &accelerate_sd()
		if (acceleration=="aitken") fun_accel = &accelerate_aitken()
		if (acceleration=="hybrid") fun_accel = &accelerate_hybrid()
		assert_msg(fun_accel != NULL, "NULL accelerate function")

		if (verbose>0) printf("{txt}   - Running solver (acceleration={res}%s{txt}, transform={res}%s{txt} tol={res}%-1.0e{txt})\n", acceleration, transform, tolerance)
		if (verbose==1) printf("{txt}   - Iterating:")
		if (verbose>1) printf("{txt}      ")

		if (this.timeit) timer_on(46)
		assert(this.solution.K == cols(data)) // we can't have this assert in -map_solver- as it will fail with parallel()
		map_solver(this, data, poolsize, fun_accel, fun_transform) // , maxiter, tolerance, verbose)
		if (this.timeit) timer_off(46)

	}
	else if (technique == "lsmr" | technique == "lsqr") {
		solution.stop_code = 0
		solution.iteration_count = 0
		solution.alphas = J(this.num_cols(), 0, .)  // 'this' is like python's self

		if (verbose > 0) printf("{txt}   - Partialling out (%s) in a pools up to size %g\n", strupper(technique), poolsize)
		for (i=1; i<=cols(data); i++) {
			if (this.timeit) timer_on(46)
			if (technique == "lsmr") {
				temp_sol = lsmr(this, data[., i], miniter, maxiter, tolerance, tolerance, ., verbose)
			}
			else {
				temp_sol = lsqr(this, data[., i], miniter, maxiter, tolerance, tolerance, ., verbose)
			}
			if (this.timeit) timer_off(46)
			assert(temp_sol.stop_code > 0)
			solver_failed = (temp_sol.stop_code >= 10)
			if (solver_failed) printf("{err}convergence not achieved in %s iterations (stop code=%g); try increasing maxiter() or decreasing tol().\n", strtrim(sprintf("%8.0fc", temp_sol.iteration_count)), temp_sol.stop_code)
			if (solver_failed & temp_sol.stop_code == 13 & this.abort==0) solver_failed = 0 // Don't exit if we set abort=0 and reached maxiter
			if (solver_failed) exit(430)

			solution.stop_code = max(( solution.stop_code , temp_sol.stop_code )) // higher number is generally worse
			solution.converged = solution.stop_code < 10
			solution.iteration_count = max(( solution.iteration_count , temp_sol.iteration_count ))
			if (this.timeit) timer_on(47)
			solution.alphas = solution.alphas , temp_sol.alphas
			data[., i] = temp_sol.data
			if (this.timeit) timer_off(47)
		}
		swap(solution.data, data)
	}
	else {
		_assert_abort(90, "ALGORITHM NOT CURRENTLY IMPLEMENTED", 1)
	}
	data = . // data should be empty now but let's make that obvious (not solution.data though!) 

	if (verbose == 0) printf(`"{txt}({browse "http://scorreia.com/research/hdfe.pdf":MWFE estimator} converged in %s iteration%s)\n"', strofreal(solution.iteration_count), solution.iteration_count > 1 ? "s" : "s")
	if (verbose > 0) printf("\n") // add space

	assert_msg(!hasmissing(solution.data), "error partialling out; missing values found")
	if (this.timeit) timer_off(40)
}


`Void' FixedEffects::save_variable(`Varname' varname,
								  `Variables' data,
								| `Varlist' varlabel)
{
	// We keep this as a method instead of standalone because it uses the .sample vector
	`RowVector'				idx
	`Integer'				i
	idx = st_addvar("double", tokens(varname))
	st_store(this.sample, idx, data)
	if (args()>=3 & varlabel!="") {
		for (i=1; i<=cols(data); i++) {
			st_varlabel(idx[i], varlabel[i])
		}
	}
}


`Void' FixedEffects::post_footnote()
{
	`Matrix'                table
	`StringVector'          rowstripe
	`StringRowVector'       colstripe
	`String'                text

	assert(!missing(G))
	st_numscalar("e(N_hdfe)", G)
	st_numscalar("e(N_hdfe_extended)", G_extended)
	st_numscalar("e(df_a)", df_a)
	st_numscalar("e(df_a_initial)", df_a_initial)
	st_numscalar("e(df_a_redundant)", df_a_redundant)
	st_numscalar("e(df_a_nested)", df_a_nested)
	st_global("e(dofmethod)", invtokens(dofadjustments))

	if (absvars == "") {
		absvars = extended_absvars = "_cons"
	}

	text = sprintf("Absorbing %g HDFE %s", G, plural(G, "group"))
	st_global("e(title2)", text)

	st_global("e(absvars)", invtokens(absvars))
	text = invtokens(extended_absvars)
	text = subinstr(text, "1.", "")
	st_global("e(extended_absvars)", text)

	// Absorbed degrees-of-freedom table
	table = (doflist_K \ doflist_M \ (doflist_K-doflist_M) \ !doflist_M_is_exact \ doflist_M_is_nested)'
	rowstripe = extended_absvars'
	rowstripe = J(rows(table), 1, "") , extended_absvars' // add equation col
	colstripe = "Categories" \ "Redundant" \ "Num Coefs" \ "Inexact?" \ "Nested?" // colstripe cannot have dots on Stata 12 or earlier
	colstripe = J(cols(table), 1, "") , colstripe // add equation col
	st_matrix("e(dof_table)", table)
	st_matrixrowstripe("e(dof_table)", rowstripe)
	st_matrixcolstripe("e(dof_table)", colstripe)

	st_numscalar("e(drop_singletons)", drop_singletons)
	st_numscalar("e(num_singletons)", num_singletons)
	
	// ivreghdfe doesn't write e(N) at this stage so we won't be able to show e(N_full)
	if (length(st_numscalar("e(N)"))) st_numscalar("e(N_full)", st_numscalar("e(N)") + num_singletons)


	// // Prune specific
	// if (prune==1) {
	// 	st_numscalar("e(pruned)", 1)
	// 	st_numscalar("e(num_pruned)", num_pruned)
	// }
	// 
	// if (!missing(finite_condition)) st_numscalar("e(finite_condition)", finite_condition)
}


`Variables' FixedEffects::project_one_fe(`Variables' y, `Integer' g)
{
	`Boolean'				store_these_alphas
	`Matrix'				alphas, proj_y

	// Cons+K+W, Cons+K, K+W, K, Cons+W, Cons = 6 variants

	// FOR NOW AT least (not working with Indiv FEs)
	assert(individual_id == "")

	store_these_alphas = storing_alphas & factors[g].save_fe
	if (store_these_alphas) assert(cols(y)==1)

	if (factors[g].num_slopes==0) {
		if (store_these_alphas) {
			factors[g].tmp_alphas = alphas = factors[g].panelmean(y, 1)
			return(alphas[factors[g].levels, .])
		}
		else {
			if (cols(y)==1 & factors[g].num_levels > 1) {
				// For speed purposes this is the most important line with technique(map):
				return(factors[g].panelmean(y, 1)[factors[g].levels])
			}
			else {
				return(factors[g].panelmean(y, 1)[factors[g].levels, .])
			}
		}
	}
	else {
		// This includes both cases, with and w/out intercept (## and #)
		if (store_these_alphas) {
			alphas = J(factors[g].num_levels, factors[g].has_intercept + factors[g].num_slopes, .)
			proj_y = panelsolve_invsym(factors[g].sort(y), factors[g], alphas)
			factors[g].tmp_alphas = alphas
			return(proj_y)
		}
		else {
			return(panelsolve_invsym(factors[g].sort(y), factors[g]))
		}
	}
}


`Void' FixedEffects::store_alphas(`Anything' d_varname)
{
	`Variable'              d

	if (verbose > 0) printf("\n{txt}## Storing estimated fixed effects (technique=%s)\n\n", technique)

	// 1) Load variable 'd'; it can be either the data or the variable name
	if (eltype(d_varname) == "string") {
		if (verbose > 0) printf("{txt}   - Loading d = e(depvar) - xb - e(resid)\n")
		d = st_data(sample, d_varname)
	}
	else {
		d = d_varname
	}
	assert(!missing(d))

	// 2) Obtain the alphas
	if (technique == "map") {
		store_alphas_map(d)
	}
	else {
		// LSMR and LSQR
		store_alphas_lsmr(d)
	}
}


`Void' FixedEffects::store_alphas_lsmr(`Variable' d)
{
	`Integer'				g, row_index, col_index, i
	`Integer'				k, n
	`Boolean'				intercept_only
	`Solution'				sol
	`Matrix'				alphas
	`StringRowVector'       varlabel
	`RowVector'             tmp_stdev
	`Real'					intercept_mean

	if (technique == "lsmr") {
		sol = lsmr(this, d, miniter, maxiter, tolerance, tolerance, ., verbose) // we want "sol.alphas"
	}
	else {
		sol = lsqr(this, d, miniter, maxiter, tolerance, tolerance, ., verbose) // we want "sol.alphas"
	}

	if (verbose > 0) printf("{txt}   - SSR of d wrt FEs: %g\n", quadcross(sol.data, sol.data))

	for (g=row_index=col_index=1; g<=G; g++) {

		k = factors[g].has_intercept + factors[g].num_slopes
		n = factors[g].num_levels * k
		intercept_only = factors[g].has_intercept& !factors[g].num_slopes

		// Get alphas corresponding to fixed effect 'g'
		alphas = sol.alphas[|row_index \ row_index + n -1|]
		row_index = row_index + n

		// 'alphas' is a matrix if we have intercept and slopes
		// TODO: maybe change the Factor_FE.mult_transpose and .mult functions so this code can be simplified
		if (factors[g].preconditioner == "none" | factors[g].preconditioner == "block_diagonal") {
			alphas = inv_vec(alphas, k)
		}
		else {
			if (!intercept_only) {
				if (factors[g].has_intercept) {
					alphas = alphas[|1\factors[g].num_levels|] , inv_vec(alphas[|factors[g].num_levels+1\.|] , k-1)
				}
				else {
					alphas = inv_vec(alphas , k)
				}
			}

		}

		factors[g].undo_preconditioning(alphas)

		if (verbose > 0) printf("{txt}   - [FE %g] Creating varlabels\n", g)
		// We'll label the FEs for ease of use
		varlabel = J(1, k, "")
		for (i=1; i<=k; i++) {
			varlabel[i] = sprintf("[FE] %s", extended_absvars[col_index])
			col_index++
		}

		if (factors[g].num_slopes) {
			if (verbose > 0) printf("{txt}   - [FE %g] Recovering unstandardized variables\n", g)
			tmp_stdev = factors[g].x_stdevs
			if (factors[g].has_intercept) tmp_stdev = 1, tmp_stdev

			// We need to *divide* the coefs by the stdev, not multiply!
			alphas = alphas :/ tmp_stdev
		}


		// The following is just for consistency, to standardize that the avg of FEs is zero
		// (applies to intercept-only FEs)
		// All cases automatically do it except LSMR/LSQR with no preconditioning
		if (factors[g].preconditioner=="none" & factors[g].has_intercept & !factors[g].num_slopes) {
			assert_boolean(factors[g].has_weights)
			if (factors[g].has_weights) {
				intercept_mean = mean(alphas[., 1], factors[g].weighted_counts)
			}
			else {
				intercept_mean = mean(alphas[., 1], factors[g].counts)

			}
			alphas[., 1] = alphas[., 1] :- intercept_mean
		}


		if (factors[g].save_fe) {
			if (verbose > 0) printf("{txt}   - [FE %g] Storing alphas in dataset\n", g)
			save_variable(factors[g].target, alphas[factors[g].levels, .], varlabel)
		}
	}
}


`Void' FixedEffects::store_alphas_map(`Variable' d)
{
	`Integer'               g, i, j
	`StringRowVector'       varlabel
	`RowVector'             tmp_stdev

	// Create empty alphas
	if (verbose > 0) printf("{txt}   - Initializing alphas\n")
	for (g=j=1; g<=G; g++) {
		if (!factors[g].save_fe) continue
		factors[g].alphas = factors[g].tmp_alphas = J(factors[g].num_levels, factors[g].has_intercept + factors[g].num_slopes, 0)
	}

	// Fill out alphas
	if (verbose > 0) printf("{txt}   - Computing alphas\n")
	this.storing_alphas = 1
	this.solution.converged = 0 // otherwise, accelerate_sd() will raise an error
	d = accelerate_sd(this, d, &transform_kaczmarz())
	this.storing_alphas = 0

	if (verbose > 0) printf("{txt}   - SSR of d wrt FEs: %g\n", quadcross(d,d))

	// Store alphas in dataset
	if (verbose > 0) printf("{txt}   - Creating varlabels\n")
	for (g=j=1; g<=G; g++) {
		if (!factors[g].save_fe) {
			j = j + factors[g].has_intercept + factors[g].num_slopes
			continue
		}
		varlabel = J(1, factors[g].has_intercept + factors[g].num_slopes, "")
		for (i=1; i<=cols(varlabel); i++) {
			varlabel[i] = sprintf("[FE] %s", extended_absvars[j])
			j++
		}

		if (factors[g].num_slopes) {
			if (verbose > 0) printf("{txt}   - Recovering unstandardized variables\n")
			tmp_stdev = factors[g].x_stdevs
			if (factors[g].has_intercept) tmp_stdev = 1, tmp_stdev

			// We need to *divide* the coefs by the stdev, not multiply!
			factors[g].alphas = factors[g].alphas :/ tmp_stdev
		}

		if (verbose > 0) printf("{txt}   - Storing alphas in dataset\n")
		save_variable(factors[g].target, factors[g].alphas[factors[g].levels, .], varlabel)
		factors[g].alphas = .
		assert(factors[g].tmp_alphas == J(0, 0, .))
	}
}


`Void' logging_singletons(`String' step, | `StringRowVector' absvars, `Optional_FE_Factor' f, `Boolean' has_intercept, `Integer' num_slopes, `Integer' num_obs, `Integer' i, `Integer' g)
{
	`Integer' 	spaces, G
	`String'	name
	spaces = max((0, max(strlen(absvars))-4))
	G = cols(absvars)

	if (step == "start") {
		printf("\n{txt}# Initializing Mata object for %g fixed %s\n\n", G, plural(G, "effect"))
		printf("{txt}   {c TLC}{hline 4}{c TT}{hline 3}{c TT}{hline 1}%s{hline 6}{c TT}{hline 6}{c TT}{hline 9}{c TT}{hline 11}{c TT}{hline 12}{c TT}{hline 9}{c TT}{hline 8}{c TT}{hline 14}{c TRC}\n", "{hline 1}" * spaces)
		printf("{txt}   {c |}  i {c |} g {c |} %s Name {c |} Int? {c |} #Slopes {c |}    Obs.   {c |}   Levels   {c |} Sorted? {c |} Indiv? {c |} #Drop Singl. {c |}\n", " " * spaces)
		printf("{txt}   {c LT}{hline 4}{c +}{hline 3}{c +}{hline 1}%s{hline 6}{c +}{hline 6}{c +}{hline 9}{c +}{hline 11}{c +}{hline 12}{c +}{hline 9}{c +}{hline 8}{c +}{hline 14}{c RT}\n", "{hline 1}" * spaces)
	}
	else if (step == "iter1"){
		name = absvars[g]
		if (name == "") name = "_cons"
		printf("{txt}   {c |} %2.0f {c |} %1.0f {c |} {res}%s{txt} {c |} ", i, g, (spaces+5-strlen(name)) * " " + name)
		printf("{txt}{%s}%3s{txt}  {c |}    %1.0f    {c |}", has_intercept ? "txt" : "err", has_intercept ? "Yes" : "No", num_slopes)
		printf("{res}%10.0g{txt} {c |}", num_obs)
	}
	else if (step == "iter2") {
		printf(" {res}%10.0g{txt} {c |} %7s {c |} %6s {c |}", f.num_levels, f.is_sorted ? "Yes" : "No", f.is_individual_fe ? "Yes" : "No")
	}
	else if (step == "end") {
		printf("{txt}   {c BLC}{hline 4}{c BT}{hline 3}{c BT}{hline 1}%s{hline 6}{c BT}{hline 6}{c BT}{hline 9}{c BT}{hline 11}{c BT}{hline 12}{c BT}{hline 9}{c BT}{hline 8}{c BT}{hline 14}{c BRC}\n", "{hline 1}" * spaces)
	}
	else {
		assert(0)
	}

	displayflush()
}

end
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
// --------------------------------------------------------------------------
// Estimate degrees-of-freedom (DoF) lost due to fixed effects
// --------------------------------------------------------------------------

mata:

`Void' estimate_dof(`FixedEffects' S, `StringRowVector' dofadjustments, `Varname' groupvar)
{
	`Integer'               g, n, i, j
	`Boolean'				skip_next
	`RowVector'				intercept_index
	`Integer'				indiv_fe_i, indiv_fe_g

	
	if (S.verbose > 0) printf("\n{txt}# Estimating degrees-of-freedom absorbed by the fixed effects\n\n")

	// Accessible properties:
	// S.factors[g] 
	//	.has_intercept
	//	.num_slopes
	//  .is_individual_fe

	// Here the focus is not on 'G' but on 'G_extended' which counts slopes and intercepts separately
	// EG: "i.turn i.trunk##c.(price gear)" has G=2 but G_extended=4: i.turn i.trunk i.trunk#c.price i.trunk#c.gear

	S.df_a_nested = 0
	
	// Compute G_extended and initialize table with number of coefficients per extended absvar
	S.doflist_K = J(1, 0, .)
	S.G_extended = 0
	for (g=1; g<=S.G; g++) {
		n = S.factors[g].has_intercept + S.factors[g].num_slopes
		S.G_extended = S.G_extended + n
		S.doflist_K = S.doflist_K , J(1, n, S.factors[g].num_levels)
	}

	// Initialize table with number of redundant coefficients for each extended absvar
	S.doflist_M = S.doflist_M_is_exact = S.doflist_M_is_nested = J(1, S.G_extended, 0)
	
	// intercept index: for each extended absvar, either 0 or 'g' (when the 'g' has intercept)
	intercept_index = J(1, S.G_extended, 0)
	for (g=i=1; g<=S.G; g++) {
		if (S.factors[g].has_intercept & !S.factors[g].is_individual_fe) intercept_index[i] = g
		if (S.factors[g].is_individual_fe) {
			assert(indiv_fe_i == .) // There can only be ONE indiv FE
			indiv_fe_i = i
			indiv_fe_g = g
		}
		n = S.factors[g].has_intercept + S.factors[g].num_slopes
		i = i + n
	}

	// Look for absvars that are clustervars or are nested within clustervars
	assert(0 <= S.num_clusters & S.num_clusters <= 10)
	if (S.num_clusters > 0 & anyof(dofadjustments, "clusters")) {
		dof_update_nested(S, intercept_index)
	}

	// (Intercept-Only) Every intercept but the first has at least one redundant coef.
	// Note that this excludes the ones nested within clusters
	// If we DO have FEs nested within clusters, we should also include the first intercept
	n = cols(selectindex(intercept_index))
	if ((n > 1) | (n > 0 & S.df_a_nested > 0)) {
		if (S.verbose > 0) printf("{txt}   - there is at least one redundant coef. for every set of FE intercepts after the first one\n")
		skip_next = S.df_a_nested == 0
		for (i=1; i<=S.G_extended; i++) {
			g = intercept_index[i]
			if (!g) continue
			if (skip_next) {
				skip_next = 0
				continue
			}
			S.doflist_M[i] = 1
		}
	}

	// First intercept absvar always has an exactly computed DoF
	for (i=1; i<=S.G_extended; i++) {
		g = intercept_index[i]
		if (!g) continue
		S.doflist_M_is_exact[i] = 1
		break
	}

	// Mobility group algorithm
	if (anyof(dofadjustments, "firstpair") | anyof(dofadjustments, "pairwise")) {
		// TODO: add group3hdfe
		dof_update_mobility_group(S, intercept_index, dofadjustments, groupvar)
	}

	if (anyof(dofadjustments, "continuous")) {
		dof_update_cvars(S)
	}

	// Algorithm for individual FEs
	if (S.individual_id != "") {
		dof_update_individual_fe(S, indiv_fe_i, indiv_fe_g)
	}

	// Store results (besides doflist_..., etc.)
	S.df_a_initial = sum(S.doflist_K)
	S.df_a_redundant = sum(S.doflist_M)
	S.df_a = S.df_a_initial - S.df_a_redundant
	S.dofadjustments = dofadjustments
	if (S.verbose > 0) printf("\n")
}


`Void' dof_update_nested(`FixedEffects' S,
						 `RowVector' intercept_index)
{
	`Integer'				i, g, i_cluster
	`String'                absvar, clustervar, vars
	`DataCol'               cluster_data
	`Factor'                F

	// (1) (Intercept-Only) Look for absvars that are clustervars
	for (i=1; i<=S.G_extended; i++) {
		g = intercept_index[i]
		if (!g) continue
		absvar = invtokens(S.factors[g].ivars, "#")
		if (anyof(S.clustervars, absvar)) {
			S.doflist_M[i] = S.doflist_K[i]
			S.df_a_nested = S.df_a_nested + S.doflist_M[i]
			S.doflist_M_is_exact[i] = S.doflist_M_is_nested[i] = 1
			intercept_index[i] = 0
			if (S.verbose > 0) printf("{txt}   - categorical variable {res}%s{txt} is also a cluster variable, so it doesn't reduce DoF\n", absvar)
		}
	}

	// (2) (Intercept-Only) Look for absvars that are nested within a clustervar
	for (i_cluster=1; i_cluster<= S.num_clusters; i_cluster++) {
		cluster_data = .
		if (!any(intercept_index)) break // no more absvars to process
		for (i=1; i<=S.G_extended; i++) {
			g = intercept_index[i]
			if (!g | S.doflist_M_is_exact[i]) continue // nothing to do
			absvar = invtokens(S.factors[g].ivars, "#")
			clustervar = S.clustervars[i_cluster]

			// Load clustervar if needed
			if (cluster_data == .) {
				if (strpos(clustervar, "#")) {
					vars = subinstr(clustervar, "#", " ", .)
					F = factor(vars, S.sample, ., "", 0, 0, ., 0)
					cluster_data = F.levels
					F = Factor() // clear
				}
				else {
					cluster_data = __fload_data(clustervar, S.sample, 0)
				}
			}

			if (S.factors[g].nested_within(cluster_data)) {
				S.doflist_M[i] = S.doflist_K[i]
				S.doflist_M_is_exact[i] = S.doflist_M_is_nested[i] = 1
				S.df_a_nested = S.df_a_nested + S.doflist_M[i]
				intercept_index[i] = 0
				if (S.verbose > 0) printf("{txt}   - categorical variable {res}%s{txt} is nested within cluster variable {res}%s{txt}, so it doesn't reduce DoF\n", absvar, clustervar)
			}
		}
	}
	cluster_data = . // save memory
}


`Void' dof_update_mobility_group(`FixedEffects' S,
								 `RowVector' intercept_index,
								 `StringRowVector' dofadjustments,
								 `Varname' groupvar)
{
	`BipartiteGraph'        bg
	`Integer'               bg_verbose                  // verbose level when calling BipartiteGraph()
	`Integer'               pair_count
	`Integer'               i, j
	`Integer'               g, h
	`Boolean'				save_subgraph
	`String'				grouplabel
	`Integer'               m                           // Mobility groups between a specific pair of FEs

	bg = BipartiteGraph()
	bg_verbose = max(( S.verbose - 1 , 0 ))
	pair_count = 0

	for (i=1; i<=S.G_extended; i++) {
		g = intercept_index[i]
		if (!g) continue

		for (j=i+1; j<=S.G_extended; j++) {
			h = intercept_index[j]
			if (!h) continue

			bg.init(&(S.factors[g]), &(S.factors[h]), bg_verbose)
			++pair_count
			save_subgraph = (pair_count == 1) & (groupvar != "")
			m = bg.init_zigzag(save_subgraph)
			if (S.verbose > 0) printf("{txt}   - mobility groups between FE intercepts #%f and #%f: {res}%f\n", g, h, m)
			if (save_subgraph) {
				if (S.verbose > 0) printf("{txt}   - Saving identifier for the first mobility group: {res}%s\n", groupvar)
				st_store(S.sample, st_addvar("long", groupvar), bg.subgraph_id)
				grouplabel = sprintf("Mobility group between %s and %s", invtokens(S.factors[g].ivars, "#"), invtokens(S.factors[h].ivars, "#"))
				st_varlabel(groupvar, grouplabel)
			}
			S.doflist_M[j] = max(( S.doflist_M[j] , m ))
			if (pair_count==1) S.doflist_M_is_exact[j] = 1
			if (pair_count & anyof(dofadjustments, "firstpair")) break
		}
		if (pair_count & anyof(dofadjustments, "firstpair")) break
	}
}


`Void' dof_update_cvars(`FixedEffects' S)
{
	`Integer'               i, j, g, k
	`Boolean'               has_int
	`Matrix'                tmp, sorted_x
	`RowVector'             zeros, results

	for (i=g=1; g<=S.G; g++) {
		// If model has intercept, redundant cvars are those that are CONSTANT
		// Without intercept, a cvar has to be zero within a FE for it to be redundant
		// Since S.fes[g].x are already demeaned IF they have intercept, we don't have to worry about the two cases
		has_int = S.factors[g].has_intercept
		if (has_int) i++
		k = S.factors[g].num_slopes
		if (!k) continue
		if (S.factors[g].is_individual_fe) continue // currently not adjusting DoFs with slopes of individual FEs
		assert(k == cols(S.factors[g].unsorted_x))
		results = J(1, k, 0)
		// float(1.1) - 1 == 2.384e-08 , so let's pick something bigger, 1e-6
		zeros = J(1, k, 1e-6)

		sorted_x = S.factors[g].sort(S.factors[g].unsorted_x)

		// BUGBUG: This part is probably VERY SLOW !!!
		// This can be speed up by moving the -if- outside the -for-
		// Maybe we can just call a function with fcollapse ???
		for (j = 1; j <= S.factors[g].num_levels; j++) {
			tmp = colminmax(panelsubmatrix(sorted_x, j, S.factors[g].info)) // TODO: this line might be SLOW!!!
			if (has_int) {
				results = results + ((tmp[2, .] - tmp[1, .]) :<= zeros)
			}
			else {
				results = results + (colsum(abs(tmp)) :<= zeros)
			}
		}

		if (sum(results)) {
			if (has_int  & (S.verbose > 0)) printf("{txt}   - the slopes in the FE #%f are constant for {res}%f{txt} levels, which don't reduce DoF\n", g, sum(results))
			if (!has_int & (S.verbose > 0)) printf("{txt}   - the slopes in the FE #%f are zero for {res}%f{txt} levels, which don't reduce DoF\n", g, sum(results))
			S.doflist_M[i..i+S.factors[g].num_slopes-1] = results
		}
		i = i + S.factors[g].num_slopes
	}
}

`Void' dof_update_individual_fe(`FixedEffects' S,
								`Integer' i,
								`Integer' g)
{
	assert(S.doflist_M[i] == 0)
	// TODO: add algorithm(s) here
	// For instance, two authors that worked at the SAME papers are redundant
	// There's probably a nice algorithm or data structure that allows us to see which indivs have the same groups
	// More generally, if indiv A always worked with B or C, and B,C only worked with A, then A is redundant
}

end
	
mata:

// --------------------------------------------------------------------------
// LSMR: least squares solver of Ax=b
// --------------------------------------------------------------------------

// Reference: 	http://web.stanford.edu/group/SOL/software/lsmr/
// - Julia:		https://github.com/JuliaMath/IterativeSolvers.jl/blob/master/src/lsmr.jl
// - Python:		https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.lsmr.html
// - Matlab:	 	https://github.com/timtylin/lsmr-SLIM/blob/master/lsmr.m
// - Fotran:		https://web.stanford.edu/group/SOL/software/lsmr/f90/lsmrf90.zip
// (note that most implementations are extremely similar to each other, up to comments)

// Usage:
// x = lsmr(A, B | , converged=?, maxiter=?, atol, btol, conlim, verbose)

// TODO:
// - Add preconditioning [PARTLY]
// - Allow A.weights [DONE]
// - Allow matrix B [DONE] -> REMOVED
// - Allow initial values for -x-
//	 This would be VERY (!!) useful when computing regressions with multiple LHSs ..
//	 like "reghdfe y x1, ... savecache(xyz)" and then reghdfe y x1 ... "usecache()".
//	 We just feed the last Xs as initial values. They don't even have to be the same varibles!
//	 So if e.g. we now include x2, we can always reuse what we had for y and x1.
// - Allow partial reorthogonalization
// - Standardize inputs???
// - Output: we want several objects:
//		- r=b-Ax
//		- in general we don't care about x, but maybe sometimes?? not sure
//		- standard errors of some xs?
//		- estimate of the ||A|| which the lsqr paper says might be useful for debugging
//		- estimate of the condition number of A (for which of the partial outs?)
//		- ...
//	  Problem: If we compute "r" maybe then compute ||r|| directly instead of relying on the product of sines formula
// - General: add -fast- option that replaces QR by the faster alternative (at the end, not part of lsmr)

// Concepts:

// residual: r = b - Ax
// relative residual: norm(r) / norm(b)  --> useful measure of how are we doing (reported by Matlab)


`Solution' lsmr(`FixedEffects' A,
			    `Matrix' b,
				`Integer' miniter,
			  | `Integer' maxiter,
			    `Real' atol,
			    `Real' btol,
			    `Real' conlim, // stop iterations if estimate of cond(A) exceeds conlim; intended to regularize ill-conditioned systems
			    `Integer' verbose)
{

	// Type definitions
	`Integer' 	m, n, iter
	`Vector' 	x
	`Solution'	sol

	// Type definitions (LSMR)
	`Real'		α, β, ρ, θ, ζ, α_bar, ρ_bar, θ_bar, ζ_bar
	`Real'		last_ρ, last_ρ_bar
	`Vector' 	u, v, h, h_bar
	`Matrix'	P, P_bar, P_tilde // Givens rotations (2x2 matrices defined by c and s)

	// Type definitions (residual norm and stopping criteria)
	`Real' 		β_hat, β_umlaut, ρ_dot, θ_tilde, β_dot, β_tilde, τ_tilde, τ_dot, ρ_tilde
	`Real' 		last_ζ, last_θ_tilde
	`Real' 		norm_A 		// Estimate of Frobenius norm ||A|| = sqrt(trace(A'A)); shouldn't exceed sqrt(n) if columns of A have been scaled to length 1
	`Real' 		norm_r 		// Estimate of ||r||
	`Real' 		norm_At_r 	// Estimate of ||A'r||
	`Real' 		norm_x 		// Estimate of ||x||
	`Real' 		cond_A 		// Estimate of cond(A)
	`Real' 		norm_b 		// Actual ||b||
	`Real' 		max_ρ_bar, min_ρ_bar, norm_A2
	`Real' 		ρ_temp
	`Real'		rel_res, rtol, rel_neq

	`String'	msg

	assert(A.weights != .)

	// Initialize constants and other scalars
	m = A.num_rows()
	n = A.num_cols()
	sol = Solution()

	// Default parameters
	if (miniter == .) miniter = 0
	if (maxiter == .) maxiter = 1e5  // Julia uses "m", Matlab uses "min(m, 20)" (recall that in exact arithmetic LSQR takes at most "m" iterations)
	if (atol == .) atol = 1e-6
	if (btol == .) btol = 1e-6
	if (conlim == .) conlim = 1e8
	if (verbose == .) verbose = 0

	// Tolerances cannot exceed roundoff error (machine precision)
	if (atol < epsilon(1)) atol = epsilon(1)
	if (btol < epsilon(1)) btol = epsilon(1)
	if (conlim == 0 | conlim > 1 / epsilon(1)) conlim = 1 / epsilon(1)

	// Sanity checks
	assert_msg(m == rows(b), "Matrix A and vector b not conformable")
	assert_msg(inrange(miniter, 0, 1e10), "miniter outside valid range: [0, 1e+10]")
	assert_msg(inrange(maxiter, 1, 1e10), "maxiter outside valid range: [1, 1e+10]")
	assert_msg(miniter <= maxiter, "miniter should be below maxiter")
	assert_msg(inrange(atol, epsilon(1), 1e-1), "atol outside valid range: [1e-16, 0.1]")
	assert_msg(inrange(btol, epsilon(1), 1e-1), "btol outside valid range: [1e-16, 0.1]")
	assert_msg(inrange(conlim, 1, 1e+16), "conlim outside valid range: [1, 1e+16]")
	assert_in(verbose, (-1, 0, 1, 2, 3, 4), "-verbose- must be an integer between -1 and 4")

	if (verbose>1) logging_lsmr("start")

	// Initial LSMR iteration
	u = normalize(b, A.weights, β=., "Failed βu = b")
	// BUGBUG: if β=0 then b=0 and we stop as any -x- solves this system
	v = normalize(A.mult_transpose(u), 1, α=., "Failed αv = A'u")
	// BUGBUG: if α=0 then normal equations already hold and we have found the solution
	
	// Initialize variables for LSMR
	h = v
	h_bar = x = J(n, 1, 0)
	α_bar = α
	ζ_bar = α * β
	ρ = ρ_bar = 1
	P_bar = I(2)

	// Initialize variables for estimation of residual norm ||r||
	β_umlaut = β
	β_dot = τ_tilde = θ_tilde = ζ = 0
	ρ_dot = 1
	ρ_tilde = 1 // This is not in the paper, but on iter=1 ρ_tilde=hypot(ρ_dot, 0)=1
	norm_b = vector_norm(b, A.weights, "Failed ||b||")

	if (norm_b == 0) {
		if (verbose>1) logging_lsmr("end")
		if (verbose>0) printf("{txt}   - LSMR converged in 0 iterations (trivial solution: 'b' is zero)\n")
		sol.stop_code = 1 // missing values
		sol.alphas = J(n, 1, 0) // could also just do "sol.alphas = x"
		sol.data = b
		sol.iteration_count = 0
		return(sol)
	}

	// Initialize variables for estimation of ||A|| and cond(A)
	norm_A = cond_A = norm_x = -1 // BUGBUG do we need to pre-initialize all these???? NOT SURE!
	norm_A2 = α ^ 2
	max_ρ_bar = 0
	min_ρ_bar= 1e100

	// Iterate Golub-Kahan
	for (iter=1; iter<=maxiter; iter++) {

		// Keep track of past values for Golub-Kahan
		last_ρ = ρ
		last_ρ_bar = ρ_bar

		// Keep track of past values for residual norm
		last_ζ = ζ
		last_θ_tilde = θ_tilde

		// 1) LSMR Algorithm:

		// Bidiagonalization
		assert_msg(!hasmissing(v), sprintf("iter !%g", iter)) // bugbug todo remove?
		u = normalize(A.mult(v) - α * u, A.weights, β=., "Failed βu = Av - αu")
		// BUGBUG We can stop if beta=0, as normal eqns hold!!!
		// if reorthogonalize: Store v in v_stack
		v = normalize(A.mult_transpose(u) - β * v, 1, α=., "Failed αv = A'u - βv")
		// if (reorthogonalize) reorthogonalize(v, V, iter)
		// BUGBUG We can stop if alpha=0, as normal eqns hold!!!

		// Construct rotation P from (α_bar, β)
		P = givens_rotation(α_bar, β, ρ=.)

		// Apply rotation P: (α, 0)  -->  (α_bar, θ)
		assign(P * (α \ 0), α_bar=., θ=.)

		// Apply rotation P_bar: (c_bar * ρ, θ)  -->  (ρ_temp, θ_bar)
		assign(P_bar * (ρ \ 0), ρ_temp=., θ_bar=.) // need to do this before overwriting P_bar below

		// Construct rotation P_bar from (ρ_temp, θ)
		P_bar = givens_rotation(ρ_temp, θ, ρ_bar=.)
		
		// Apply rotation P_bar: (0, ζ_bar) -->  (ζ_bar, ζ)
		assign(P_bar * (0 \ ζ_bar), ζ_bar=., ζ=.)

		// Update h, h_bar, x
		h_bar = h - (θ_bar * ρ) / (last_ρ * last_ρ_bar) * h_bar
		x = x + ζ / (ρ * ρ_bar) * h_bar
		h = v - (θ / ρ) * h

		// 2) Compute residual norm ||r||:

		// Apply rotation P: (0, β_umlaut)  -->  (β_umlaut, β_hat)
		assign(P * (0 \ β_umlaut), β_umlaut=., β_hat=.)

		if (iter > 1) {
			// Construct rotation P_tilde from (ρ_dot, θ_bar)
			P_tilde = givens_rotation(ρ_dot, θ_bar, ρ_tilde=.)
			
			// Apply rotation P_tilde: (ρ_bar, 0)  -->  (ρ_dot, θ_tilde)
			assign(P_tilde * (ρ_bar \ 0), ρ_dot=., θ_tilde=.)

			// Apply rotation P_tilde: (β_hat, β_dot)  -->  (β_dot, β_tilde)
			assign(P_tilde * (β_hat \ β_dot), β_dot=., β_tilde=.) // Note that "β_tilde" is never used
		}

		// Update t_tilde by forward substitution
		τ_tilde = (last_ζ - last_θ_tilde * τ_tilde) / ρ_tilde
		τ_dot = (ζ - θ_tilde * τ_tilde) / ρ_dot

		// 3) Estimate norms:

		// Estimate ||r||
		norm_r = hypot(β_dot - τ_dot, β_umlaut)
		
		// Estimate ||A||
		norm_A2 = norm_A2 + β ^ 2
		norm_A = sqrt(norm_A2)
		norm_A2 = norm_A2 + α ^ 2

		// Estimate cond(A)
		max_ρ_bar = max((max_ρ_bar, last_ρ_bar))
		if (iter > 1) min_ρ_bar = min((min_ρ_bar, last_ρ_bar))
		cond_A = max((max_ρ_bar, ρ_temp)) / min((min_ρ_bar, ρ_temp))
		
		// Estimate ||error|| = ||A'r||
		norm_At_r = abs(ζ_bar)

		// Estimate ||x||
		norm_x = vector_norm(x, 1, "Failed ||x||")

		// 4) Test for convergence
		rel_res = norm_r / norm_b
		rtol = btol + atol * norm_A * norm_x / norm_b
		rel_neq = norm_At_r / norm_A / norm_x
		if (verbose>1) logging_lsmr("iter", iter, rel_res, rtol, rel_neq, atol, 1/cond_A, 1/conlim, norm_A)

		if (iter < miniter) {
			continue
		}

		if (rel_res == . | rel_neq == .) {
			msg = sprintf("{err}@ LSMR stopped; missing values in rel_res or rel_neq\n")
			sol.stop_code = 12 // missing values
			break
		}
		else if (rel_res <= rtol) {
			msg = sprintf("{txt}   - LSMR converged in %g iterations (criterion: relative residual)\n", iter)
			sol.stop_code = 2 // consistent systems (with exact solutions)
			break
		}
		else if (rel_neq <= atol) {
			msg = sprintf("{txt}   - LSMR converged in %g iterations (criterion: normal equation)\n", iter)
			sol.stop_code = 3 // inconsistent systems (least squares)
			break
		}
		else if (cond_A >= conlim) {
			msg = sprintf("{err}@ LSMR stopped; A is ill-conditioned\n")
			sol.stop_code = 11 // ill-conditioned
			break
		}
	}

	if (sol.stop_code == 0) {
		assert_msg(iter == maxiter+1)
		iter = maxiter
		msg = sprintf("{err}@ LSMR stopped; maximum number of iterations reached\n")
		sol.stop_code = 13 // max-iter reached
	}

	if (verbose>1) logging_lsmr("end")
	if (verbose>0 | (verbose>-1 & sol.stop_code>=10)) printf(msg)
	
	sol.data = b - A.mult(x) // BUGBUG: this will use A LOT of space for a sec: 1) "x" 2) A.mult(x) 3) b 4) the substraction; sadly we don't have in-place substract
	swap(sol.alphas, x) // BUGBUG: we need to adjust the alphas to undo any preconditioning ; see https://web.stanford.edu/group/SOL/software/lsmr/
	sol.iteration_count = iter
	return(sol)
}


`Void' logging_lsmr(`String' step, | `Integer' iter, `Real' rel_res, `Real' rtol, `Real' rel_neq, `Real' atol, `Real' invcond, `Real' ctol, `Real' norm_A)
{
	`Integer' 	col
	`String' 	table_row, color1, color2, color3

	// See "help smcl##ascii"
	if (step == "start") {
		printf("\n{txt}## Solving linear system via LSMR:\n")
		printf(" {txt}{c TLC}{hline 5}{c TT}{hline 22}{c TT}{hline 22}{c TT}{hline 22}{c TT}{hline 10}{c TRC}\n")
		printf(" {txt}{c |}{space 5}{c |}        Rule 1        {c |}        Rule 2        {c |}        Rule 3        {c |}{space 10}{c |}\n")
		printf(" {txt}{c |}  i  {c LT}{hline 12}{c TT}{hline 9}{c +}{hline 12}{c TT}{hline 9}{c +}{hline 12}{c TT}{hline 9}{c RT}   ||A||  {c |}\n")
		printf(" {txt}{c |}{space 5}{c |}  rel.res.  {c |}   rtol  {c |}  rel.neq.  {c |}   atol  {c |}  1/cond(A) {c |}   ctol  {c |}{space 10}{c |}\n")
		printf(" {txt}{c LT}{hline 5}{c +}{hline 12}{c +}{hline 9}{c +}{hline 12}{c +}{hline 9}{c +}{hline 12}{c +}{hline 9}{c +}{hline 10}{c RT}\n")
	}
	else if (step == "end") {
		printf(" {txt}{c BLC}{hline 5}{c BT}{hline 12}{c BT}{hline 9}{c BT}{hline 12}{c BT}{hline 9}{c BT}{hline 12}{c BT}{hline 9}{c BT}{hline 10}{c BRC}\n")
	}
	else {
		col = 0
		color1 = (rel_res <= rtol) ? "{res}" : "{txt}"
		color2 = (rel_neq <= atol) ? "{res}" : "{txt}"
		color3 = (invcond <= ctol) ? "{res}" : "{txt}"
		table_row = sprintf("{txt} {c |} %3.0f {col %3.0f}{c |}", iter, col = col + 8)
		
		table_row = table_row + sprintf("%s%11.0e{txt}{col %3.0f}{c |}", color1, rel_res, col = col + 13)
		table_row = table_row + sprintf("%6.0e{col %3.0f}{c |}", rtol, col = col + 10)
		//table_row = table_row + sprintf("%8.1f{col %3.0f}{c |}", 100 * rtol / rel_res, col = col + 10)
		
		table_row = table_row + sprintf("%s%11.0e{txt}{col %3.0f}{c |}", color2, rel_neq, col = col + 13)
		table_row = table_row + sprintf("%6.0e{col %3.0f}{c |}", atol, col = col + 10)
		//table_row = table_row + sprintf("%8.1f{col %3.0f}{c |}", 100 * atol / rel_neq, col = col + 10)
		
		table_row = table_row + sprintf("%s%11.0e{txt}{col %3.0f}{c |}", color3, invcond, col = col + 13)
		table_row = table_row + sprintf("%6.0e{col %3.0f}{c |}", ctol, col = col + 10)
		//table_row = table_row + sprintf("%8.1f{col %3.0f}{c |}", 100 * ctol / invcond, col = col + 10)
		
		table_row = table_row + sprintf("%9.1f {col %3.0f}{c |}", norm_A, col = col + 11)
		printf(table_row + "\n")
	}

	displayflush()
}

end
mata:

// --------------------------------------------------------------------------
// LSQR: least squares solver of Ax=b
// --------------------------------------------------------------------------

// - Method derived from: https://web.stanford.edu/group/SOL/software/lsqr/lsqr-toms82a.pdf
// - Code uses scaffolding created by SAC in LSMR.mata
// - Most code adapted from Julia LSQR code 
// - Julia code: https://github.com/JuliaMath/IterativeSolvers.jl/blob/master/src/lsqr.jl
// - Note: If atol & btol are 1e-9, residual norm should be accurate to about 9 digits


`Solution' lsqr(`FixedEffects' A, 
				`Matrix' b,
				`Integer' miniter,
			  | `Integer' maxiter,
				`Real' atol,
				`Real' btol,
				`Real' ctol,
				`Integer' verbose)
{


	// Essential type definitions
	`Integer' 	m, n, iter
	`Vector' 	x
	`Solution'	sol

	// Type definitions for LSQR
	`Real'		α, β, ρ, ρ1, θ, ρ_bar, φ, φ_bar, φ_bar1, ζ, ζ_bar
	`Real'		rhs, δ, γ_bar, γ, cs2, sn2, extra_var
	`Vector' 	u, v, w
	`Matrix'	P, P_bar

	// Type definitions for residual norm and stopping criteria
	`Real' 		norm_r, norm_A, norm_b, norm_x, norm_r1, norm_r2, norm_dd, norm_xx, cond_A
	`Real' 		norm_Ar, res1, res2, r1sq, test1, test2 
	`Real' 		ρ_temp, w_ρ, cs1, sn1, ψ, τ, t1, t2, test3, rtol 

	`String'	msg
	
	assert(A.weights != .)

	// Initialize constants and other scalars
	m = A.num_rows()
	n = A.num_cols()
	sol = Solution()

	// Default parameters
	if (miniter == .)	miniter = 0
	if (maxiter == .) 	maxiter = 1e5
	if (atol == .) 		atol = 1e-6
	if (btol == .) 		btol = 1e-6
	if (ctol == .) 		ctol = 1e-6
	if (verbose == .) 	verbose = 0

	// Tolerances cannot exceed roundoff error (machine precision)
	if (atol < epsilon(1)) atol = epsilon(1)
	if (btol < epsilon(1)) btol = epsilon(1)
	if (ctol < epsilon(1)) ctol = epsilon(1)

	// Sanity checks
	assert_msg(m == rows(b), "Matrix A and vector b not conformable")
	assert_msg(inrange(miniter, 0, 1e10), "miniter outside valid range: [0, 1e+10]")
	assert_msg(inrange(maxiter, 1, 1e10), "maxiter outside valid range: [1, 1e+10]")
	assert_msg(miniter <= maxiter, "miniter should be below maxiter")
	assert_msg(inrange(atol, 1e-16, 1e-1), "atol outside valid range: [1e-16, 0.1]")
	assert_msg(inrange(btol, 1e-16, 1e-1), "btol outside valid range: [1e-16, 0.1]")
	assert_msg(inrange(ctol, 1e-16, 1e-1), "ctol outside valid range: [1e-16, 0.1]")
	assert_in(verbose, (-1, 0, 1, 2, 3, 4), "-verbose- must be an integer between -1 and 4")
	if (verbose>1) logging_lsqr("start")
	// if (verbose > 0) printf("\n{txt}## Solving linear system via LSQR:\n")

	// Initial LSQR iteration
	u = normalize(b, A.weights, β=., "Failed βu = b") // βu = b 
	v = normalize(A.mult_transpose(u), A.weights, α=., "Failed αv = A'u") // αv = A'u 

	// New section -- if norm(b) is zero, just exit 
	norm_b = vector_norm(b, A.weights, "Failed ||b||")
	if (norm_b == 0) {
		if (verbose>1) logging_lsqr("end")
		if (verbose>0) printf("{txt}   - LSQR converged in 0 iterations (trivial solution: 'b' is zero)\n")
		sol.stop_code = 1 // missing values
		sol.alphas = J(n, 1, 0) // could also just do "sol.alphas = x"
		sol.data = b
		sol.iteration_count = 0
		return(sol)
	}
	
	// Initialize vars for LSQR
	x 		= J(n, 1, 0)
	w 		= v 
	norm_A 	= cond_A = norm_dd = res2 = norm_x = norm_xx = sn2 = ζ = 0
	cs2 	= -1
	norm_Ar = α * β
	ρ_bar 	= α
	φ_bar	= norm_b = norm_r = norm_r1 = norm_r2 = β
	
	// Start the loop
	for (iter=1; iter<=maxiter; iter++) {

		// Lanczos process: generate vector v and scalars α, β 

		// Step 3: Bidiagonalization --> copied code from Sergio 
		u = normalize(A.mult(v) - α * u, A.weights, β=., "Failed βu = Av - αu") // βu = Av - αu
		v = normalize(A.mult_transpose(u) - β * v, A.weights, α=., "Failed αv = A'u - βv") // αv = A'u - βv 
		
		// Step 4: Construct & apply next orthogonal transformation
		// The following section is the calculation of x we need 
 
		φ_bar1  = φ_bar 
		ρ 		= hypot(ρ_bar, β)

		// givens rotations function returns (c, -s \ s, c)
		P = givens_rotation(ρ_bar, β, ρ=.)

		// apply rotations
		// multiple rotations here because I need some extra values (compared to LSMR)
		assign(P * (0 \ -α),     θ = ., ρ_bar = .)

		assign(P * (φ_bar1 \ 0), φ = ., φ_bar = .) // we need φ for updating t1, τ 
		assign(P * (φ \ 0)     , τ = ., extra_var = .) // need τ in the calculation of norm_Ar 


		// Calculate some vars following Givens Rotation 
		t1 		= φ / ρ
		t2 		=-θ / ρ

		// update w, x 
		x 		= t1 * w + x 
		w 		= t2 * w + v 

		// Now we move to constructing norms that are required for stopping criteria 
		// Use a plane rotation on the right to elimatine the 
		// super-diagonal element (θ) of the upper-bidiaganol matrix to estimate norm(x)

		// Want to use Givens Rotations instead of keeping cs2, sn2 
		γ_bar = -ρ	
		P_bar 	= givens_rotation(γ_bar, θ, γ = .)


		assign(P_bar * (0 \ -ρ), δ = ., γ_bar = .)


		rhs 	= φ - (δ * ζ)
		ζ_bar	= rhs / γ_bar

		norm_x 	= sqrt(norm_xx + (ζ_bar^2))
		γ 		= hypot(γ_bar, θ)

		ζ 		= rhs / γ
		norm_xx = norm_xx + ζ^2 

		// Estimate cond(A_bar), ||r_bar||, and ||A_bar'r_bar||
		w_ρ     = w * (1/ρ)
		norm_dd = norm_dd + norm(w_ρ)

		norm_A 	= sqrt(norm_A^2 + α^2 + β^2) // can't use hypot() as written bc of 3 args 
		
		cond_A 	= norm_A * sqrt(norm_dd)
		res1 	= φ_bar ^ 2	
		norm_Ar = α * abs(τ)
		norm_r 	= sqrt(res1)

		// Distinguish between norms 
		r1sq 	= norm_r^2
		norm_r1 = sqrt(r1sq)
		norm_r2 = norm_r

		// Use these norms to estimate other quantities
		// These quantities --> 0 near a solution 
		test1	= norm_r / norm_b
		test2 	= norm_Ar / (norm_A * norm_r)
		test3 	= 1 / cond_A
		t1 		= test1 / (1 + norm_A*norm_x/norm_b)
		rtol	= btol + atol*norm_A*norm_x/norm_b
		if (verbose>1) logging_lsqr("iter", iter, test1, rtol, test2, atol, test3, ctol, norm_A)

		if (iter < miniter) {
			continue
		}

		if (test1 == . | test2 == .) {
			msg = sprintf("{err}@ LSQR stopped; missing values in rel_res or rel_neq\n")
			sol.stop_code = 12 // missing values
			break
		}
		else if (test1 <= rtol) {
			msg = sprintf("{txt}   - LSQR converged in %g iterations (criterion: relative residual)\n", iter)
			sol.stop_code = 2 // consistent systems (with exact solutions)
			break
		}
		else if (test2 <= atol) {
			msg = sprintf("{txt}   - LSQR converged in %g iterations (criterion: normal equation)\n", iter)
			sol.stop_code = 3 // inconsistent systems (least squares)
			break
		}
		else if (test3 <= ctol) {
			msg = sprintf("{err}@ LSQR stopped; A is ill-conditioned\n")
			sol.stop_code = 11 // ill-conditioned
			break
		}
	}
	if (sol.stop_code == 0) {
		assert_msg(iter == maxiter+1)
		iter = maxiter
		msg = sprintf("{err}@ LSQR stopped; maximum number of iterations reached\n")
		sol.stop_code = 13 // max-iter reached
	}
	if (verbose>1) logging_lsqr("end")
	if (verbose>0 | (verbose>-1 & sol.stop_code>=10)) printf(msg)
	
	sol.data = b - A.mult(x) // BUGBUG: this will use A LOT of space for a sec: 1) "x" 2) A.mult(x) 3) b 4) the substraction; sadly we don't have in-place substract
	swap(sol.alphas, x) // BUGBUG: we need to adjust the alphas to undo any preconditioning ; see https://web.stanford.edu/group/SOL/software/lsmr/
	sol.iteration_count = iter
	return(sol)
}


`Void' logging_lsqr(`String' step, | `Integer' iter, `Real' test1, `Real' rtol, `Real' test2, `Real' atol, `Real' test3, `Real' ctol, `Real' norm_A){
	`Integer' 	col
	`String' 	table_row, color1, color2, color3
	// See "help smcl##ascii"
	if (step == "start") {
		printf("\n{txt}## Solving linear system via LSQR:\n")
		printf(" {txt}{c TLC}{hline 5}{c TT}{hline 22}{c TT}{hline 22}{c TT}{hline 22}{c TT}{hline 10}{c TRC}\n")
		printf(" {txt}{c |}{space 5}{c |}        Rule 1        {c |}        Rule 2        {c |}        Rule 3        {c |}{space 10}{c |}\n")
		printf(" {txt}{c |}  i  {c LT}{hline 12}{c TT}{hline 9}{c +}{hline 12}{c TT}{hline 9}{c +}{hline 12}{c TT}{hline 9}{c RT}   ||A||  {c |}\n")
		printf(" {txt}{c |}{space 5}{c |}  rel.res.  {c |}   rtol  {c |}  rel.neq.  {c |}   atol  {c |}  1/cond(A) {c |}   ctol  {c |}{space 10}{c |}\n")
		printf(" {txt}{c LT}{hline 5}{c +}{hline 12}{c +}{hline 9}{c +}{hline 12}{c +}{hline 9}{c +}{hline 12}{c +}{hline 9}{c +}{hline 10}{c RT}\n")
	}
	else if (step == "end") {
		printf(" {txt}{c BLC}{hline 5}{c BT}{hline 12}{c BT}{hline 9}{c BT}{hline 12}{c BT}{hline 9}{c BT}{hline 12}{c BT}{hline 9}{c BT}{hline 10}{c BRC}\n")
	}
	else {
		col = 0
		color1 = (test1 <= rtol) ? "{res}" : "{txt}"
		color2 = (test2 <= atol) ? "{res}" : "{txt}"
		color3 = (test3 <= ctol) ? "{res}" : "{txt}"
		table_row = sprintf("{txt} {c |} %3.0f {col %3.0f}{c |}", iter, col = col + 8)
		
		table_row = table_row + sprintf("%s%11.0e{txt}{col %3.0f}{c |}", color1, test1, col = col + 13)
		table_row = table_row + sprintf("%6.0e{col %3.0f}{c |}", rtol, col = col + 10)
				
		table_row = table_row + sprintf("%s%11.0e{txt}{col %3.0f}{c |}", color2, test2, col = col + 13)
		table_row = table_row + sprintf("%6.0e{col %3.0f}{c |}", atol, col = col + 10)
		
		table_row = table_row + sprintf("%s%11.0e{txt}{col %3.0f}{c |}", color3, test3, col = col + 13)
		table_row = table_row + sprintf("%6.0e{col %3.0f}{c |}", ctol, col = col + 10)
		
		table_row = table_row + sprintf("%9.1f {col %3.0f}{c |}", norm_A, col = col + 11)
		printf(table_row + "\n")
	}
	displayflush()
}

end


// --------------------------------------------------------------------------
// Main code for MAP solver
// --------------------------------------------------------------------------

mata:

`Void' map_solver(`FixedEffects' S,
					  `Matrix' data,
					  `Integer' poolsize,
					  `FunctionP' fun_accel,
					  `FunctionP' fun_transform)
{
	`Integer'				i
	
		S.solution.converged = 0 // converged will get updated by check_convergence()

	// Speedup for constant-only case (no fixed effects)
	if (S.G==1 & S.factors[1].method=="none" & !S.factors[1].num_slopes & !(S.storing_alphas & S.factors[1].save_fe)) {
		assert(S.factors[1].num_levels == 1)
		// Not using poolsize here as this case is unlikely to happen with very large datasets
		data = data :- mean(data, S.factors[1].weights)
		S.solution.converged = 1
		S.solution.iteration_count = 1
	}
	else {
		if (poolsize >= cols(data)) {
			data = (*fun_accel)(S, data, fun_transform)
		}
		else {
			assert(S.solution.converged == 0) // check_convergence() sets converged=1 so we need to undo it
			for (i=1; i<=S.solution.K; i++) {
				S.solution.converged = 0
				data[., i] = (*fun_accel)(S, data[., i], fun_transform)
			}
		}
	}

	swap(S.solution.data, data)
}

end
// --------------------------------------------------------------------------
// MAP Acceleration Schemes
// --------------------------------------------------------------------------

mata:

`Variables' function accelerate_test(`FixedEffects' S, `Variables' y, `FunctionP' T) {
	`Integer'				iter, g
	`Variables'				resid
	`FE_Factor'				f
	pragma unset resid

	assert_msg(S.solution.converged == 0, "solution had already converged")

	for (iter=1; iter<=S.maxiter; iter++) {
		for (g=1; g<=S.G; g++) {
			f = S.factors[g]
			if (g==1) resid = y - f.panelmean(y, 1)[f.levels, .]
			else resid = resid - f.panelmean(resid, 1)[f.levels, .]
		}
		if (check_convergence(S, iter, resid, y)) break
		y = resid
	}
	return(resid)
}

// --------------------------------------------------------------------------

`Variables' function accelerate_none(`FixedEffects' S, `Variables' y, `FunctionP' T) {
	`Integer'				iter
	`Variables'				resid
	pragma unset resid
	assert_msg(S.solution.converged == 0, "solution had already converged")

	for (iter=1; iter<=S.maxiter; iter++) {
		(*T)(S, y, resid) // Faster version of "resid = S.T(y)"
		if (check_convergence(S, iter, resid, y)) break
		y = resid
	}
	return(resid)
}
// --------------------------------------------------------------------------

// Start w/out acceleration, then switch to CG
`Variables' function accelerate_hybrid(`FixedEffects' S, `Variables' y, `FunctionP' T) {
	`Integer' iter, accel_start
	`Variables' resid
	pragma unset resid

	accel_start = 6
	assert_msg(S.solution.converged == 0, "solution had already converged")

	for (iter=1; iter<=accel_start; iter++) {
		(*T)(S, y, resid) // Faster version of "resid = S.T(y)"
		if (check_convergence(S, iter, resid, y)) break
		y = resid
	}

	T = &transform_sym_kaczmarz() // Override

	return(accelerate_cg(S, y, T))
}

// --------------------------------------------------------------------------
// Memory cost is approx = 4*size(y) (actually 3 since y is already there)
// But we need to add maybe 1 more due to u:*v
// And I also need to check how much does project and T use..
// Double check with a call to memory

// For discussion on the stopping criteria, see the following presentation:
// Arioli & Gratton, "Least-squares problems, normal equations, and stopping criteria for the conjugate gradient method". URL: https://www.stfc.ac.uk/SCD/resources/talks/Arioli-NAday2008.pdf

// Basically, we will use the Hestenes and Stiefel rule

`Variables' function accelerate_cg(`FixedEffects' S, `Variables' y, `FunctionP' T) {
	// BUGBUG iterate the first 6? without acceleration??
	`Integer'	iter, d, Q
	`Variables'		r, u, v
	`RowVector' alpha, beta, ssr, ssr_old, improvement_potential
	`Matrix' recent_ssr
	pragma unset r
	pragma unset v

	assert_msg(S.solution.converged == 0, "solution had already converged")
	Q = cols(y)
	
	d = 1 // BUGBUG Set it to 2/3 // Number of recent SSR values to use for convergence criteria (lower=faster & riskier)
	// A discussion on the stopping criteria used is described in
	// http://scicomp.stackexchange.com/questions/582/stopping-criteria-for-iterative-linear-solvers-applied-to-nearly-singular-system/585#585

	improvement_potential = weighted_quadcolsum(S, y, y)
	recent_ssr = J(d, Q, .)
	
	(*T)(S, y, r, 1) // this counts as "iteration 0"
	ssr = weighted_quadcolsum(S, r, r) // cross(r,r) when cols(y)==1 // BUGBUG maybe diag(quadcross()) is faster?
	u = r

	for (iter=1; iter<=S.maxiter; iter++) {
		(*T)(S, u, v, 1) // This is the hottest loop in the entire program
		alpha = safe_divide( ssr , weighted_quadcolsum(S, u, v) )
		recent_ssr[1 + mod(iter-1, d), .] = alpha :* ssr
		improvement_potential = improvement_potential - alpha :* ssr
		y = y - alpha :* u
		//if (S.compute_rre & !S.prune) reghdfe_rre_benchmark(y[., 1], S.rre_true_residual, S.rre_depvar_norm)
		r = r - alpha :* v
		ssr_old = ssr
		if (S.verbose>=5) r
		ssr = weighted_quadcolsum(S, r, r)
		beta = safe_divide( ssr , ssr_old) // Fletcher-Reeves formula, but it shouldn't matter in our problem
		u = r + beta :* u
		// Convergence if sum(recent_ssr) > tol^2 * improvement_potential
		if ( check_convergence(S, iter, colsum(recent_ssr), improvement_potential, "hestenes") ) {
			break
		}
	}
	return(y)
}

// --------------------------------------------------------------------------

`Variables' function accelerate_sd(`FixedEffects' S, `Variables' y, `FunctionP' T) {
	`Integer'	iter, g
	`Variables' proj
	`RowVector' t
	pragma unset proj

	assert_msg(S.solution.converged == 0, "solution had already converged")

	for (iter=1; iter<=S.maxiter; iter++) {
		(*T)(S, y, proj, 1)
		if (check_convergence(S, iter, y-proj, y)) break
		t = safe_divide( weighted_quadcolsum(S, y, proj) , weighted_quadcolsum(S, proj, proj) )
		if (uniform(1,1)<0.1) t = 1 // BUGBUG: Does this REALLY help to randomly unstuck an iteration?

		y = y - t :* proj
		//if (S.compute_rre & !S.prune) reghdfe_rre_benchmark(y[., 1], S.rre_true_residual, S.rre_depvar_norm)

		if (S.storing_alphas) {
			for (g=1; g<=S.G; g++) {
				if (S.factors[g].save_fe) {
					S.factors[g].alphas = S.factors[g].alphas + t :* S.factors[g].tmp_alphas
				}
			}
		}
	}
	
	// Clean up temp object
	if (S.storing_alphas) {
		for (g=1; g<=S.G; g++) {
			if (S.factors[g].save_fe) {
				S.factors[g].tmp_alphas = J(0, 0, .)
			}
		}
	}

	return(y-proj)
}

// --------------------------------------------------------------------------
// This is method 3 of Macleod (1986), a vector generalization of the Aitken-Steffensen method
// Also: "when numerically computing the sequence.. stop..  when rounding errors become too 
// important in the denominator, where the ^2 operation may cancel too many significant digits"
// Note: Sometimes the iteration gets "stuck"; can we unstuck it with adding randomness
// in the accelerate decision? There should be a better way.. (maybe symmetric kacz instead of standard one?)

`Variables' function accelerate_aitken(`FixedEffects' S, `Variables' y, `FunctionP' T) {
	`Integer'	iter
	`Variables'		resid, y_old, delta_sq
	`Boolean'	accelerate
	`RowVector' t
	pragma unset resid

	assert_msg(S.solution.converged == 0, "solution had already converged")
	y_old = J(rows(y), cols(y), .)

	for (iter=1; iter<=S.maxiter; iter++) {
		
		(*T)(S, y, resid)
		accelerate = iter>=S.accel_start & !mod(iter,S.accel_freq)

		// Accelerate
		if (accelerate) {
			delta_sq = resid - 2 * y + y_old // = (resid - y) - (y - y_old) // Equivalent to D2.resid
			// t is just (d'd2) / (d2'd2)
			t = safe_divide( weighted_quadcolsum(S,  (resid - y) , delta_sq) ,  weighted_quadcolsum(S, delta_sq , delta_sq) )
			resid = resid - t :*  (resid - y)
		}


		// Only check converge on non-accelerated iterations
		// BUGBUG: Do we need to disable the check when accelerating?
		//if (S.compute_rre & !S.prune) reghdfe_rre_benchmark(resid[., 1], S.rre_true_residual, S.rre_depvar_norm)
		if (check_convergence(S, iter, resid, y)) break
		y_old = y // y_old is resid[iter-2]
		y = resid // y is resid[iter-1]
	}
	return(resid)
}

// --------------------------------------------------------------------------

`Boolean' check_convergence(`FixedEffects' S, `Integer' iter, `Variables' y_new, `Variables' y_old,| `String' method) {
	`Boolean'	is_last_iter
	`Real'		update_error
	`Real'		eps_threshold

	// max() ensures that the result when bunching vars is at least as good as when not bunching
	if (args()<5) method = "vectors"

	if (S.G==1 & !S.storing_alphas) {
		// Shortcut for trivial case (1 FE)
		update_error = 0
	}
	else if (method=="vectors") {
		update_error = max(mean(reldif(y_new, y_old), S.weights))
	}
	else if (method=="hestenes") {
		// If the regressor is perfectly explained by the absvars, we can have SSR very close to zero but negative
		// (so sqrt is missing)

		eps_threshold = 1e-15 // 10 * epsilon(1) ; perhaps too aggressive and should be 1e-14 ?
		if (S.verbose > 0 & all(y_new :< eps_threshold)) {
			printf("{txt} note: eps. is very close to zero (%g), so hestenes assumed convergence to avoid numerical precision errors\n", min(y_new))
		}
		update_error = safe_divide(edittozerotol(y_new, eps_threshold ),
		                           editmissing(y_old, epsilon(1)),
		                           epsilon(1) )
		update_error = sqrt(max(update_error))
	}
	else {
		exit(error(100))
	}

	assert_msg(!missing(update_error), "update error is missing")

	S.solution.converged = S.solution.converged + (update_error <= S.tolerance)
	is_last_iter = iter==S.maxiter
	
	if (S.solution.converged >= S.min_ok) {
		S.solution.iteration_count = max((iter, S.solution.iteration_count))
		S.solution.accuracy = max((S.solution.accuracy, update_error))
		if (S.verbose==1) printf("{txt}   converged in %g iterations (last error =%3.1e)\n", iter, update_error)
		if (S.verbose>1) printf("\n{txt}   - Converged in %g iterations (last error =%3.1e)\n", iter, update_error)
	}
	else if (is_last_iter & S.abort) {
		printf("\n{err}convergence not achieved in %g iterations (last error=%e); try increasing maxiter() or decreasing tol().\n", S.maxiter, update_error)
		exit(430)
	}
	else {
		if ((S.verbose>=2 & S.verbose<=3 & mod(iter,1)==0) | (S.verbose==1 & mod(iter,1)==0)) {
			printf("{res}.{txt}")
			displayflush()
		}
		if ( (S.verbose>=2 & S.verbose<=3 & mod(iter,100)==0) | (S.verbose==1 & mod(iter,100)==0) ) {
			printf("{txt}%9.1f\n      ", update_error/S.tolerance)
		}

		if (S.verbose==4 & method!="hestenes") printf("{txt} iter={res}%4.0f{txt}\tupdate_error={res}%-9.6e\n", iter, update_error)
		if (S.verbose==4 & method=="hestenes") printf("{txt} iter={res}%4.0f{txt}\tupdate_error={res}%-9.6e  {txt}norm(ssr)={res}%g\n", iter, update_error, norm(y_new))
		
		if (S.verbose>=5) {
			printf("\n{txt} iter={res}%4.0f{txt}\tupdate_error={res}%-9.6e{txt}\tmethod={res}%s\n", iter, update_error, method)
			"old:"
			y_old
			"new:"
			y_new
		}
	}
	return(S.solution.converged >= S.min_ok)
}

// --------------------------------------------------------------------------

`Matrix' weighted_quadcolsum(`FixedEffects' S, `Matrix' x, `Matrix' y) {
	// BUGBUG: override S.has_weights with pruning
	// One approach is faster for thin matrices
	// We are using cross instead of quadcross but it should not matter for this use
	if (S.has_weights) {
		if (cols(x) < 14) {
			return(quadcross(x :* y, S.weights)')
		}
		else {
			return(diagonal(quadcross(x, S.weights, y))')
		}
	}
	else {
		if (cols(x) < 25) {
			return(diagonal(quadcross(x, y))')
		}
		else {
			return(colsum(x :* y))
		}
	}
}


// RRE benchmarking
//	|| yk - y || / || y || === || ek - e || / || y ||
`Real' reghdfe_rre_benchmark(`Vector' resid, `Vector' true_resid, `Real' norm_y) {
	`Real' ans
	ans = norm(resid - true_resid) / norm_y
	return(ans)
}

end
// --------------------------------------------------------------------------
// Transformations: Compute RESIDUALS, not projections
// --------------------------------------------------------------------------

mata:

`Void' function transform_cimmino(`FixedEffects' S, `Variables' y, `Variables' ans,| `Boolean' get_proj) {
	`Integer' 	g
	if (args()<4 | get_proj==.) get_proj = 0
	assert_boolean(get_proj)

	if (get_proj & S.G == 2) {
		ans = (S.project_one_fe(y, 1) + S.project_one_fe(y, 2)) / S.G
	}
	else if (get_proj & S.G == 3) {
		ans = (S.project_one_fe(y, 1) + S.project_one_fe(y, 2) + S.project_one_fe(y, 3)) / S.G
	}
	else {
		// General case
		ans = S.project_one_fe(y, 1)
		for (g=2; g<=S.G; g++) {
			ans = ans + S.project_one_fe(y, g)
		}
		ans = get_proj ? ans / S.G : y - ans / S.G
	}
}

// --------------------------------------------------------------------------

`Void' function transform_kaczmarz(`FixedEffects' S, `Variables' y, `Variables' ans,| `Boolean' get_proj) {
	`Integer' 	g

	if (args()<4 | get_proj==.) get_proj = 0

	ans = y - S.project_one_fe(y, 1)
	for (g=2; g<=S.G; g++) {
		ans = ans - S.project_one_fe(ans, g)
	}
	if (get_proj) ans = y - ans
}

// --------------------------------------------------------------------------
// This seems slower than kaczmarz (sym kaczmarz!); not used currently
`Void' function transform_rand_kaczmarz(`FixedEffects' S, `Variables' y, `Variables' ans,| `Boolean' get_proj) {
	`Integer' 	g
	`Vector' rand
	if (args()<4 | get_proj==.) get_proj = 0
	rand = sort( ( (1::S.G) , uniform(S.G,1) ) , 2 )[.,1]
	ans = y - S.project_one_fe(y, rand[1])
	for (g=2; g<=S.G; g++) {
		ans = ans - S.project_one_fe(ans, rand[g])
	}
	for (g=S.G-1; g>=1; g--) {
		ans = ans - S.project_one_fe(ans, rand[g])
	}
	if (get_proj) ans = y - ans
}

// --------------------------------------------------------------------------

`Void' function transform_sym_kaczmarz(`FixedEffects' S, `Variables' y, `Variables' ans,| `Boolean' get_proj) {
	`Integer' 	g
	if (args()<4 | get_proj==.) get_proj = 0
	ans = y - S.project_one_fe(y, 1)
	for (g=2; g<=S.G; g++) {
		ans = ans - S.project_one_fe(ans, g)
	}
	for (g=S.G-1; g>=1; g--) {
		ans = ans - S.project_one_fe(ans, g)
	}
	if (get_proj) ans = y - ans
}

end
	
// --------------------------------------------------------------------------
// Code for Multiprocessing (Parallel Processing)
// --------------------------------------------------------------------------

/* WARNINGS:
- There are several bugs when saving/loading mata classes into disk
  so we need to be very careful to avoid them. These include:

  1) No associative arrays i.e. asarray()
  	For instance, we must set
  	HDFE.factors[g].vl = HDFE.factors[g].extra = .

  2) If the class has transmorphic attributes, Stata might crash depending on the value
  	 (if it's a number it seems it doesn't crash, but if its a string it might)
*/

mata:

// --------------------------------------------------------------------------
// Clean up data for parallel
// --------------------------------------------------------------------------
`Void' cleanup_for_parallel(`FixedEffects' HDFE)
{
	`Integer'				g

	if (HDFE.verbose>0) printf("\n{txt}# Cleaning up the HDFE object so it can be saved/loaded from disk\n\n")
	
	// extra is an asarray() used by reghdfe v5. Unless we remove the asarray, Stata will hard crash when saving or loading the HDFE object
	// vl is an asarray() used by fcollapse. Unless we remove the asarray, Stata will hard crash when saving or loading the HDFE object

	for (g=1; g<=HDFE.G; g++) {
		HDFE.factors[g].cleanup_before_saving()
		//HDFE.factors[g].bg.PF1 = NULL
		//HDFE.factors[g].bg.PF2 = NULL
	}

	// if (HDFE.parallel_opts == "") HDFE.parallel_opts = J(0, 0, "")
	// if (HDFE.parallel_dir == "") HDFE.parallel_dir = J(0, 0, "")

	//HDFE.weight_type = J(0, 0, "")
	//HDFE.weight_var = J(0, 0, "")
}


// --------------------------------------------------------------------------
// Save data for parallel processing
// --------------------------------------------------------------------------
`Void' save_before_parallel(`String' parallel_dir,
							`FixedEffects' HDFE,
							`Matrix' data)
{
	`Integer'               verbose
	`Integer'               parallel_poolsize, parallel_numproc
	`Integer'               num_cols, left_index, right_index, remaining_cols, remaining_workers
	`String'				fn
	`Integer'               fh

	verbose = HDFE.verbose
	num_cols = cols(data)
	assert(HDFE.parallel_maxproc > 1 | HDFE.parallel_force == 1) // We should never run only ONE worker process

	// parallel_numproc can be below parallel_maxproc if we don't have enough variables

	parallel_poolsize = ceil(num_cols / HDFE.parallel_maxproc)
	if (verbose > 0) printf("\n{txt}## [Parallel] Loading and partialling %g variables using up to %g worker processes\n", num_cols, HDFE.parallel_maxproc)
	if (verbose > 0) printf("{txt}              Each process will work in blocks of %g-%g variables\n", parallel_poolsize-1, parallel_poolsize)
	if (verbose > 0) printf("{txt}              Temporary files will be saved in %s\n", parallel_dir)

	// Save HDFE object
	mkdir(parallel_dir, 1) // Need to create it before calling -parallel_map- (1=Public)
	fn = pathjoin(parallel_dir, "data0.tmp")
	fh = fopen(fn, "w")
	fputmatrix(fh, HDFE)
	fclose(fh)
	if (verbose > 0) printf("{txt}              - HDFE object saved in %s\n", fn)

	HDFE.parallel_poolsizes = J(1, 0, .)

	// Save data objects
	for (left_index=parallel_numproc=1; left_index<=num_cols; parallel_numproc++) {
		remaining_cols = num_cols - left_index + 1
		remaining_workers = HDFE.parallel_maxproc - parallel_numproc + 1
		parallel_poolsize = ceil(remaining_cols / remaining_workers)
		right_index = min((left_index+parallel_poolsize-1, num_cols))

		fn = pathjoin(parallel_dir, sprintf("data%f.tmp", parallel_numproc))
		fh = fopen(fn, "w")
		fputmatrix(fh, data[., left_index..right_index])
		fclose(fh)

		if (verbose > 0) printf("{txt}              - Data block #%f with %f cols saved in %s\n", parallel_numproc, right_index-left_index+1, fn)
		left_index = right_index + 1
		HDFE.parallel_numproc = parallel_numproc
		HDFE.parallel_poolsizes = HDFE.parallel_poolsizes , parallel_poolsize
	}
	if (verbose > 0) printf("\n")
	assert(HDFE.parallel_numproc <= HDFE.parallel_maxproc)
}


// --------------------------------------------------------------------------
// Load, partial out, and save data
// --------------------------------------------------------------------------
`Void' worker_partial_out(`String' hdfe_fn, `String' data_fn)
{
	`FixedEffects'			HDFE
	`Matrix'				data
	`Integer'               fh

	fh = fopen(hdfe_fn, "r")
	HDFE = fgetmatrix(fh, 1)
	fclose(fh)
	
	fh = fopen(data_fn, "r")
	data = fgetmatrix(fh, 1)
	fclose(fh)

	inner_worker_partial_out(HDFE, data)

	unlink(data_fn)
	fh = fopen(data_fn, "w")
	fputmatrix(fh, HDFE.solution.data)
	fclose(fh)

	// Uncomment to debug
	//stata(sprintf("sleep %4.0f", runiform(1,1)*10000))
}


`Void' inner_worker_partial_out(`FixedEffects' HDFE,
								`Matrix' data)
{
	// Based on ::partial_out()

	`FunctionP'				fun_transform
	`FunctionP'				fun_accel
	`String'				technique
	`String'				transform
	`String'				acceleration
	`Integer'				verbose
	`Integer'				poolsize
	`Solution'				temp_sol
	`Integer'				i
	`Boolean'				solver_failed

	technique = HDFE.technique
	transform = HDFE.transform
	acceleration = HDFE.acceleration
	verbose = HDFE.verbose
	poolsize = HDFE.poolsize

	if (technique == "map") {
		// Load transform pointer
		if (transform=="cimmino") fun_transform = &transform_cimmino()
		if (transform=="kaczmarz") fun_transform = &transform_kaczmarz()
		if (transform=="symmetric_kaczmarz") fun_transform = &transform_sym_kaczmarz()
		if (transform=="random_kaczmarz") fun_transform = &transform_rand_kaczmarz() // experimental

		// Pointer to acceleration routine
		if (acceleration=="test") fun_accel = &accelerate_test()
		if (acceleration=="none") fun_accel = &accelerate_none()
		if (acceleration=="conjugate_gradient") fun_accel = &accelerate_cg()
		if (acceleration=="steepest_descent") fun_accel = &accelerate_sd()
		if (acceleration=="aitken") fun_accel = &accelerate_aitken()
		if (acceleration=="hybrid") fun_accel = &accelerate_hybrid()

		if (verbose>0) printf("{txt}   - Running solver (acceleration={res}%s{txt}, transform={res}%s{txt} tol={res}%-1.0e{txt})\n", acceleration, transform, HDFE.tolerance)
		if (verbose==1) printf("{txt}   - Iterating:")
		if (verbose>1) printf("{txt}      ")

		map_solver(HDFE, data, poolsize, fun_accel, fun_transform) // , maxiter, tolerance, verbose)

	}
	else if (technique == "lsmr" | technique == "lsqr") {

		if (verbose > 0) printf("{txt}   - Partialling out (%s) in a pools up to size %g\n", strupper(technique), poolsize)
		for (i=1; i<=cols(data); i++) {
			if (technique == "lsmr") {
				temp_sol = lsmr(HDFE, data[., i], HDFE.miniter, HDFE.maxiter, HDFE.tolerance, HDFE.tolerance, ., verbose)
			}
			else {
				temp_sol = lsqr(HDFE, data[., i], HDFE.miniter, HDFE.maxiter, HDFE.tolerance, HDFE.tolerance, ., verbose)
			}
			assert(temp_sol.stop_code > 0)
			solver_failed = (temp_sol.stop_code >= 10)
			if (solver_failed) printf("{err}convergence not achieved in %s iterations (stop code=%g); try increasing maxiter() or decreasing tol().\n", strtrim(sprintf("%8.0fc", temp_sol.iteration_count)), temp_sol.stop_code)
			if (solver_failed & temp_sol.stop_code == 13 & HDFE.abort==0) solver_failed = 0 // Don't exit if we set abort=0 and reached maxiter
			if (solver_failed) {
				data = J(rows(data), cols(data), .)
				exit(430)
			}

			//solution.stop_code = max(( solution.stop_code , temp_sol.stop_code )) // higher number is generally worse
			//solution.converged = solution.stop_code < 10
			//solution.iteration_count = max(( solution.iteration_count , temp_sol.iteration_count ))
			data[., i] = temp_sol.data
		}
		swap(HDFE.solution.data, data)
	}
	else {
		_assert_abort(90, "ALGORITHM NOT CURRENTLY IMPLEMENTED", 1)
	}

}


`Void' parallel_combine(`FixedEffects' HDFE)
{
	`Integer'				num_cols, proc, left_index, right_index
	`String'				fn
	`Integer'				fh
	`Matrix'				data_slice

	// Trick to get # of columns
	num_cols = cols(HDFE.solution.means)
	assert(num_cols == cols(HDFE.solution.norm2))

	assert(cols(HDFE.parallel_poolsizes) == HDFE.parallel_numproc)

	HDFE.solution.data = J(HDFE.N, num_cols, .)

	// Join back results
	left_index = 1
	for (proc=1; proc<=HDFE.parallel_numproc; proc++)
	{
		right_index = left_index + HDFE.parallel_poolsizes[proc] - 1
		if (HDFE.verbose > 0) printf("{txt}              - Loading data block #%f with %f cols from %s\n", proc, right_index-left_index+1, fn)

		fn = pathjoin(HDFE.parallel_dir, sprintf("data%f.tmp", proc))
		fh = fopen(fn, "r")
		data_slice = fgetmatrix(fh, 1)
		fclose(fh)

		assert_msg(cols(data_slice) == right_index - left_index + 1)

		HDFE.solution.data[., left_index..right_index] = data_slice
		left_index = right_index + 1
	}

	if (HDFE.verbose == 0) printf(`"{txt}({browse "http://scorreia.com/research/hdfe.pdf":MWFE estimator} converged in %s iteration%s)\n"', "??", "s")
	if (HDFE.verbose > 0) printf("\n") // add space

	assert_msg(!hasmissing(HDFE.solution.data), "error partialling out; missing values found")
}

end

exit
