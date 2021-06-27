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
