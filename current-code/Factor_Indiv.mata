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
