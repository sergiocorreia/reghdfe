// --------------------------------------------------------------------------
// FixedEffects main class
// --------------------------------------------------------------------------

mata:

class FixedEffects
{
	// Factors
	`Integer'               G                   // Number of sets of FEs
	`Integer'               N                   // number of obs
	`Integer'               M                   // Sum of all possible FE coefs
	`Factors'               factors
	`Vector'                sample
	`Varlist'               absvars
	`Varlist'               ivars
	`Varlist'               cvars
	`Boolean'               has_intercept
	`RowVector'             intercepts
	`RowVector'             num_slopes
	`Integer'               num_singletons
	`Boolean'               save_any_fe
	`Boolean'               save_all_fe
	`Varlist'               targets
	`RowVector'             save_fe

	// Constant-related (also see -has_intercept-)
	`Boolean'               report_constant
	`Boolean'               compute_constant

	// Optimization options
	`Real'                  tolerance
	`Real'                  extra_tolerance		// Try to achieve this tol if it only takes a few more iters: ceil(10%)
	`Integer'               maxiter
	`String'                transform           // Kaczmarz Cimmino Symmetric_kaczmarz (k c s)
	`String'                acceleration        // Acceleration method. None/No/Empty is none\
	`Integer'               accel_start         // Iteration where we start to accelerate // set it at 6? 2?3?
	`string'                slope_method
	`Boolean'               prune               // Whether to recursively prune degree-1 edges
	`Boolean'               abort               // Raise error if convergence failed?
	`Integer'               accel_freq          // Specific to Aitken's acceleration
	`Boolean'               storing_alphas      // 1 if we should compute the alphas/fes
	`Real'                  conlim              // specific to LSMR
	`Real'                  btol                // specific to LSMR
	`Boolean'				always_run_lsmr_preconditioner
	`Integer'				min_ok

	// Optimization objects
	`BipartiteGraph'        bg                  // Used when pruning 1-core vertices
	`Vector'                pruned_weight       // temp. weight for the factors that were pruned
	`Integer'               prune_g1            // Factor 1/2 in the bipartite subgraph that gets pruned
	`Integer'               prune_g2            // Factor 2/2 in the bipartite subgraph that gets pruned
	`Integer'               num_pruned          // Number of vertices (levels) that were pruned

	// Misc
	`Integer'               verbose
	`Boolean'               timeit
	`Boolean'               compact
	`Integer'               poolsize
	`Boolean'               store_sample
	`Real'                  finite_condition
	`Real'                  compute_rre         // Relative residual error: || e_k - e || / || e ||
	`Real'                  rre_depvar_norm
	`Vector'                rre_varname
	`Vector'                rre_true_residual
	`String'                panelvar
	`String'                timevar

	`RowVector'             not_basevar         // Boolean vector indicating whether each regressor is or not a basevar
	`String'                fullindepvars       // indepvars including basevars

	// Weight-specific
	`Boolean'               has_weights
	`Variable'              weight              // unsorted weight
	`String'                weight_var          // Weighting variable
	`String'                weight_type         // Weight type (pw, fw, etc)

	// Absorbed degrees-of-freedom computations
	`Integer'               G_extended          // Number of intercepts plus slopes
	`Integer'               df_a_redundant      // e(mobility)
	`Integer'               df_a_initial
	`Integer'               df_a                // df_a_inital - df_a_redundant
	`Vector'                doflist_M
	`Vector'                doflist_K
	`Vector'                doflist_M_is_exact
	`Vector'                doflist_M_is_nested
	`Vector'                is_slope
	`Integer'               df_a_nested // Redundant due to bein nested; used for: r2_a r2_a_within rmse

	// VCE and cluster variables
	`String'                vcetype
	`Integer'               num_clusters
	`Varlist'               clustervars
	`Varlist'               base_clustervars
	`String'                vceextra

	// Regression-specific
	`String'                varlist             // y x1 x2 x3
	`String'                depvar              // y
	`String'                indepvars           // x1 x2 x3
	`String'                tousevar
	
	`Boolean'               drop_singletons
	`String'                absorb              // contents of absorb()
	`String'                select_if           // If condition
	`String'                select_in           // In condition
	`String'                model               // ols, iv
	`String'                summarize_stats
	`Boolean'               summarize_quietly
	`StringRowVector'       dofadjustments // firstpair pairwise cluster continuous
	`Varname'               groupvar
	`String'                residuals
	`Variable'              residuals_vector
	`RowVector'             kept // 1 if the regressors are not deemed as omitted (by partial_out+cholsolve+invsym)
	`String'                diopts

	// Output
	`String'                cmdline
	`String'                subcmd
	`String'                title
	`Boolean'               converged
	`Integer'               iteration_count // e(ic)
	`Varlist'               extended_absvars
	`String'                notes
	`Integer'               df_r
	`Integer'               df_m
	`Integer'               N_clust
	`Integer'               N_clust_list
	`Real'                  rss
	`Real'                  rmse
	`Real'                  F
	`Real'                  tss
	`Real'                  tss_within
	`Real'                  sumweights
	`Real'                  r2
	`Real'                  r2_within
	`Real'                  r2_a
	`Real'                  r2_a_within
	`Real'                  ll
	`Real'                  ll_0
	`Real'                  accuracy
	`RowVector'             means
	`RowVector'				all_stdevs

	// Methods
	`Void'                  new()
	`Void'                  destroy()
	`Void'                  load_weights() // calls update_sorted_weights, etc.
	`Void'                  update_sorted_weights()
	`Void'                  update_cvar_objects()
	`Matrix'                partial_out()
	`Matrix'                partial_out_pool()
	`Void'                  _partial_out()
	`Variables'             project_one_fe()
	`Void'                  prune_1core()
	`Void'                  _expand_1core()
	`Void'                  estimate_dof()
	`Void'                  estimate_cond()
	`Void'                  save_touse()
	`Void'                  store_alphas()
	`Void'                  save_variable()
	`Void'                  post_footnote()
	`Void'                  post()
	`FixedEffects'          reload() // create new instance of object

	// LSMR-Specific Methods
	`Real'                  lsmr_norm()
	`Vector'                lsmr_A_mult()
	`Vector'                lsmr_At_mult()
}    


// Set default value of properties
`Void' FixedEffects::new()
{
	num_singletons = .
	sample = J(0, 1, .)
	weight = 1 // set to 1 so cross(x, S.weight, y)==cross(x, y)

	verbose = 0
	timeit = 0
	compact = 0
	poolsize = .
	finite_condition = .
	residuals = ""
	residuals_vector = .
	panelvar = timevar = ""
	iteration_count = 0
	accuracy = -1 // Epsilon at the time of convergence

	// Optimization defaults
	slope_method = "invsym"
	maxiter = 1e4
	tolerance = 1e-8
	transform = "symmetric_kaczmarz"
	acceleration = "conjugate_gradient"
	accel_start = 6
	conlim = 1e+8 // lsmr only
	btol = 1e-8 // lsmr only (note: atol is just tolerance)
	always_run_lsmr_preconditioner = 0
	min_ok = 1

	prune = 0
	converged = 0
	abort = 1
	storing_alphas = 0
	report_constant = compute_constant = 1

	// Specific to Aitken:
	accel_freq = 3

	not_basevar = J(1, 0, .)

	means = all_stdevs = J(1, 0, .) // necessary with pool() because we append to it
	kept = J(1, 0, .) // necessary with pool() because we append to it
}


`Void' FixedEffects::destroy()
{
	// stata(sprintf("cap drop %s", tousevar))
}


// This adds/removes weights or changes their type
`Void' FixedEffects::load_weights(`String' weighttype, `String' weightvar, `Variable' weight, `Boolean' verbose)
{
	`Integer'				g
	`FactorPointer'         pf
	`Matrix'                precond // used for lsmr
	`Varname'               cvars_g
   
	this.has_weights = (weighttype != "" & weightvar != "")
	if (this.verbose > 0 & verbose > 0 & this.has_weights) printf("{txt}## Loading weights [%s=%s]\n", weighttype, weightvar)

	// Update main properties
	this.weight_var = weightvar
	this.weight_type = weighttype

	// Update booleans
	for (g=1; g<=this.G; g++) {
		asarray(this.factors[g].extra, "has_weights", this.has_weights)
	}

	// Optionally load weight from dataset
	if (this.has_weights & weight==J(0,1,.)) {
		weight = st_data(this.sample, this.weight_var)
	}

	// Update weight vectors
	if (this.has_weights) {
		if (this.verbose > 0 & verbose > 0) printf("{txt}## Sorting weights for each absvar\n")
		this.update_sorted_weights(weight)
	}
	else {
		// If no weights, clear this up
		this.weight = 1 // same as defined by new()
		for (g=1; g<=this.G; g++) {
			asarray(this.factors[g].extra, "weights", .)
			asarray(this.factors[g].extra, "weighted_counts", .)
		}
	}
	
	// Update cvar objects (do AFTER updating weights!)
	// (this is meaningless with iweights)
	if (weighttype != "iweight") this.update_cvar_objects()

	// Preconditioners for LSMR
	if (acceleration=="lsmr" | always_run_lsmr_preconditioner) {

		// Compute M
		M = 0
		for (g=1; g<=G; g++) {
			M = M + factors[g].num_levels * (intercepts[g] + num_slopes[g])
		}

		// Preconditioner
		for (g=1; g<=G; g++) {
			pf = &(factors[g])
			if (intercepts[g]) {
				precond = has_weights ? asarray((*pf).extra, "weighted_counts") : (*pf).counts
				asarray((*pf).extra, "precond_intercept", sqrt(1 :/ precond))
			}

			if (num_slopes[g]) {
				cvars_g = tokens(this.cvars[g])
				precond = st_data(this.sample, cvars_g)
				precond = reghdfe_panel_precondition(precond, (*pf))
				asarray((*pf).extra, "precond_slopes", precond)
			}

			precond = .
		}
	}

}


// This just updates the weight but doesn't change the type or variable of the weight
`Void' FixedEffects::update_sorted_weights(`Variable' weight)
{
	`Integer'               g
	`Real'                  min_w
	`Variable'              w
	`FactorPointer'         pf

	assert_msg(!hasmissing(weight), "weights can't be missing")
	this.weight = weight
	assert(rows(weight)==rows(sample))
	if (verbose > 0) printf("{txt}   - loading %s weight from variable %s\n", weight_type, weight_var)
	for (g=1; g<=G; g++) {
		if (verbose > 0) printf("{txt}   - sorting weight for factor {res}%s{txt}\n", absvars[g])
		pf = &(factors[g])
		w = (*pf).sort(weight)

		// Rescale weights so there are no weights below 1
		if (weight_type != "fweight") {
			min_w = colmin(w)
			if (min_w < 1e-6) min_w = 1e-6 // Prevent bugs if a weight is very close to zero
			//assert_msg(min_w > 0, "weights must be positive")
			//if (min_w <= 0) printf("{err} not all weights are positive\n")
			if (0 < min_w & min_w < 1) {
				w = w :/ min_w
			}
		}

		asarray((*pf).extra, "weights", w)
		asarray((*pf).extra, "weighted_counts", `panelsum'(w, (*pf).info))
	}
}


`Void' FixedEffects::update_cvar_objects()
{
	`Integer'               g
	`FactorPointer'         pf

	for (g=1; g<=G; g++) {
		pf = &(factors[g])
		// Update mean(z; w) and inv(z'z; w) where z is a slope variable and w is the weight
		if (num_slopes[g]) {
			if (verbose > 0) printf("{txt}   - precomputing cvar objects for factor {res}%s{txt}\n", absvars[g])
			if (intercepts[g]) {
			    asarray((*pf).extra, "xmeans",
			            panelmean(asarray((*pf).extra, "x"), *pf))
			}
			asarray((*pf).extra, "inv_xx", precompute_inv_xx(*pf, intercepts[g]))
		}
	}
}


`Variables' FixedEffects::partial_out(`Anything' data,
									| `Boolean' save_tss,
									  `Boolean' standardize_data,
									  `Boolean' first_is_depvar)
{
	// -data- is either a varlist or a matrix
	`Variables'             y
	`Varlist'               vars
	`Integer'               i
	`Integer'               k

	if (args()<2 | save_tss==.) save_tss = 0
	if (args()<3 | standardize_data==.) standardize_data = 1
	if (args()<4 | first_is_depvar==.) first_is_depvar = 1

	if (eltype(data) == "string") {
		vars = tokens(invtokens(data)) // tweak to allow string scalars and string vectors
		k = cols(vars)

		if (poolsize < k) {
			if (verbose > 0) printf("\n{txt}## Loading and partialling out %g variables in blocks of %g\n\n", k, poolsize)
			if (timeit) timer_on(50)
			partial_out_pool(vars, save_tss, standardize_data, first_is_depvar, poolsize, y=.)
			if (timeit) timer_off(50)
		}
		else {
			if (verbose > 0) printf("\n{txt}## Partialling out %g variables: {res}%s{txt}\n\n", cols(vars), invtokens(vars))
			if (verbose > 0) printf("{txt}   - Loading variables into Mata\n")
			if (timeit) timer_on(50)
			_st_data_wrapper(sample, invtokens(vars), y=., verbose)
			if (timeit) timer_off(50)

			// Workaround to odd Stata quirk
			if (timeit) timer_on(51)
			if (cols(y) > cols(vars)) {
				printf("{err}(some empty columns were added due to a bug/quirk in {bf:st_data()}; %g cols created instead of %g for {it:%s}; running slower workaround)\n", cols(y), cols(vars), invtokens(vars))
				partial_out_pool(vars, save_tss, standardize_data, first_is_depvar, 1, y=.)
			}
			else {
				_partial_out(y, save_tss, standardize_data, first_is_depvar)
			}
			if (timeit) timer_off(51)
			
		}
	}
	else {
		if (verbose > 0) printf("\n{txt}## Partialling out %g variables\n\n", cols(data))
		if (timeit) timer_on(54)
		_partial_out(y=data, save_tss, standardize_data, first_is_depvar)
		if (timeit) timer_off(54)
	}

	if (verbose==0) printf(`"{txt}({browse "http://scorreia.com/research/hdfe.pdf":MWFE estimator} converged in %s iteration%s)\n"', strofreal(iteration_count), iteration_count > 1 ? "s" : "s")
	return(y)
}



`Variables' FixedEffects::partial_out_pool(`Anything' vars,
										   `Boolean' save_tss,
										   `Boolean' standardize_data,
										   `Boolean' first_is_depvar,
										   `Integer' step,
										   `Variables' y)
{
	`Variables'             part_y
	`Integer'               i, j, ii
	`Integer'               k
	`StringRowVector'       keepvars

	k = cols(vars)
	assert(step > 0)
	assert(step < k)
	y = J(rows(sample), 0, .)

	for (i=1; i<=k; i=i+step) {
		
		j = i + step - 1
		if (j>k) j = k

		// Load data
		_st_data_wrapper(sample, vars[i..j], part_y=., verbose)

		if (cols(part_y) > j - i + 1) {
			printf("{err}(some empty columns were added due to a bug/quirk in {bf:st_data()}; running slower workaround)\n")
			if (timeit) timer_on(51)
			part_y = J(rows(sample), 0, .)
			for (ii=i; ii<=j; ii++) {
				part_y = part_y, st_data(sample, vars[ii])
			}
			if (timeit) timer_off(51)
		}

		// Drop loaded vars as quickly as possible
		if (compact) {
			// st_dropvar(vars[i..j]) // bugbug what if repeated??
			keepvars = base_clustervars , timevar, panelvar, (j == k ? "" : vars[j+1..k])
			keepvars = tokens(invtokens(keepvars))
			if (cols(keepvars)) {
				stata(sprintf("fvrevar %s, list", invtokens(keepvars)))
				stata(sprintf("keep %s", st_global("r(varlist)")))
			}
			else {
				stata("clear")
			}
		}

		_partial_out(part_y, save_tss, standardize_data, first_is_depvar)
		y = y, part_y
		part_y = .
	}
}


`Void' FixedEffects::store_alphas(`Anything' d_varname)
{
	`Integer'               g, i, j
	`StringRowVector'       varlabel
	`Variable'              d
	`RowVector'             tmp_stdev

	if (verbose > 0) printf("\n{txt}## Storing estimated fixed effects\n\n")

	// -d- can be either the data or the variable name

	// Load -d- variable
	if (eltype(d_varname) == "string") {
		if (verbose > 0) printf("{txt}   - Loading d = e(depvar) - xb - e(resid)\n")
		d = st_data(sample, d_varname)
	}
	else {
		d = d_varname
	}
	assert(!missing(d))

	// Create empty alphas
	if (verbose > 0) printf("{txt}   - Initializing alphas\n")
	for (g=j=1; g<=G; g++) {
		if (!save_fe[g]) continue
		asarray(factors[g].extra, "alphas", J(factors[g].num_levels, intercepts[g] + num_slopes[g], 0))
		asarray(factors[g].extra, "tmp_alphas", J(factors[g].num_levels, intercepts[g] + num_slopes[g], 0))
	}

	// Fill out alphas
	if (verbose > 0) printf("{txt}   - Computing alphas\n")
	storing_alphas = 1
	converged = 0
	d = accelerate_sd(this, d, &transform_kaczmarz())
	storing_alphas = 0

	if (verbose > 0) printf("{txt}   - SSR of d wrt FEs: %g\n", quadcross(d,d))

	// Store alphas in dataset
	if (verbose > 0) printf("{txt}   - Creating varlabels\n")
	for (g=j=1; g<=G; g++) {
		if (!save_fe[g]) {
			j = j + intercepts[g] + num_slopes[g]
			continue
		}
		varlabel = J(1, intercepts[g] + num_slopes[g], "")
		for (i=1; i<=cols(varlabel); i++) {
			varlabel[i] = sprintf("[FE] %s", extended_absvars[j])
			j++
		}

		if (num_slopes[g]) {
			if (verbose > 0) printf("{txt}   - Recovering unstandardized variables\n")
			tmp_stdev = asarray(factors[g].extra, "x_stdevs")
			if (intercepts[g]) tmp_stdev = 1, tmp_stdev

			// We need to *divide* the coefs by the stdev, not multiply!
			asarray(factors[g].extra, "alphas",
				asarray(factors[g].extra, "alphas") :/ tmp_stdev
			)
		}

		if (verbose > 0) printf("{txt}   - Storing alphas in dataset\n")
		save_variable(targets[g], asarray(factors[g].extra, "alphas")[factors[g].levels, .], varlabel)
		asarray(factors[g].extra, "alphas", .)
		asarray(factors[g].extra, "tmp_alphas", .)
	}
}


`Void' FixedEffects::_partial_out(`Variables' y,
								| `Boolean' save_tss,
								  `Boolean' standardize_data,
								  `Boolean' first_is_depvar,
								  `Boolean' flush)
{
	`RowVector'             stdevs, needs_zeroing, kept2
	`FunctionP'             funct_transform, func_accel // transform
	`Real'                  y_mean, collinear_tol
	`Vector'                lhs
	`Vector'                alphas
	`StringRowVector'       vars
	`Integer'               i

	if (args()<2 | save_tss==.) save_tss = 0
	if (args()<3 | standardize_data==.) standardize_data = 1
	if (args()<4 | first_is_depvar==.) first_is_depvar = 1
	if (args()<5 | flush==.) flush = 0 // whether or not to flush the values of means & kept

	assert(anyof((0, 1, 2), standardize_data)) // 0=Don't standardize; 1=Std. and REVERT after partial; 2=Std., partial, and KEEP STANDARDIZED

	if (flush) {
		iteration_count = 0
		accuracy = -1
		means = stdevs = J(1, 0, .)
		kept = J(1, 0, .)
	}

	// Shortcut for trivial case (1 FE)
	if (G==1) acceleration = "none"

	// Solver Warnings
	if (transform=="kaczmarz" & acceleration=="conjugate_gradient") {
		printf("{err}(WARNING: convergence is {bf:unlikely} with transform=kaczmarz and accel=CG)\n")
	}

	// Load transform pointer
	if (transform=="cimmino") funct_transform = &transform_cimmino()
	if (transform=="kaczmarz") funct_transform = &transform_kaczmarz()
	if (transform=="symmetric_kaczmarz") funct_transform = &transform_sym_kaczmarz()
	if (transform=="random_kaczmarz") funct_transform = &transform_rand_kaczmarz()

	// Pointer to acceleration routine
	if (acceleration=="test") func_accel = &accelerate_test()
	if (acceleration=="none") func_accel = &accelerate_none()
	if (acceleration=="conjugate_gradient") func_accel = &accelerate_cg()
	if (acceleration=="steepest_descent") func_accel = &accelerate_sd()
	if (acceleration=="aitken") func_accel = &accelerate_aitken()
	if (acceleration=="hybrid") func_accel = &accelerate_hybrid()

	// Compute TSS of depvar
	if (timeit) timer_on(60)
	if (save_tss & tss==.) {
		lhs = y[., 1]
		if (has_intercept) {
			y_mean = mean(lhs, weight)
			tss = crossdev(lhs, y_mean, weight, lhs, y_mean) // Sum of w[i] * (y[i]-y_mean) ^ 2
		}
		else {
			tss = cross(lhs, weight, lhs) // Sum of w[i] * y[i] ^ 2
		}
		lhs = .
		if (weight_type=="aweight" | weight_type=="pweight") tss = tss * rows(y) / sum(weight)
	}
	if (timeit) timer_off(60)


	// Compute 2-norm of each var, to see if we need to drop as regressors
	kept2 = diagonal(cross(y, y))'

	// Compute and save means of each var
	means = means , ( compute_constant ? mean(y, weight) : J(1, cols(y), 1) )

	// Intercept LSMR case
	if (acceleration=="lsmr") {
		// RRE benchmarking
		if (compute_rre) rre_depvar_norm = norm(y[., 1])
		if (cols(y)==1) {
			y = lsmr(this, y, alphas=.)
			alphas = . // or return them!
		}
		else {
			for (i=1; i<=cols(y); i++) {
				y[., i] = lsmr(this, y[., i], alphas=.)
			}
			alphas = .
		}
	}
	else {

		// Standardize variables
		if (timeit) timer_on(61)
		if (standardize_data) {
			if (verbose > 0) printf("{txt}   - Standardizing variables\n")
			stdevs = reghdfe_standardize(y)
			all_stdevs = all_stdevs, stdevs
			kept2 = kept2 :/ stdevs :^ 2
		}
		if (timeit) timer_off(61)

		// RRE benchmarking
		if (compute_rre) {
			rre_true_residual = rre_true_residual / (standardize_data ? stdevs[1] : 1)
			rre_depvar_norm = norm(y[., 1])
		}

		// Solve
		if (verbose>0) printf("{txt}   - Running solver (acceleration={res}%s{txt}, transform={res}%s{txt} tol={res}%-1.0e{txt})\n", acceleration, transform, tolerance)
		if (verbose==1) printf("{txt}   - Iterating:")
		if (verbose>1) printf("{txt}      ")
		converged = 0 // converged will get updated by check_convergence()

		if (timeit) timer_on(62)
		if (G==1 & factors[1].method=="none" & num_slopes[1]==0 & !(storing_alphas & save_fe[1])) {
			// Speedup for constant-only case (no fixed effects)
			assert(factors[1].num_levels == 1)
			iteration_count = 1
			accuracy = 0
			if (standardize_data == 1) {
				y = stdevs :* y :- stdevs :* mean(y, has_weights ? asarray(factors[1].extra, "weights") : 1) // Undoing standardization
			}
			else {
				y = y :- mean(y, has_weights ? asarray(factors[1].extra, "weights") : 1)
			}
		}
		else {
			if (standardize_data == 1) {
				y = (*func_accel)(this, y, funct_transform) :* stdevs // Undoing standardization
			}
			else {
				y = (*func_accel)(this, y, funct_transform) // 'this' is like python's self
			}
		}
		if (timeit) timer_off(62)
		
		if (prune) {
			assert_msg(G==2, "prune option requires only two FEs")
			if (timeit) timer_on(63)
			_expand_1core(y)
			if (timeit) timer_off(63)
		}
	}

	assert_msg(!hasmissing(y), "error partialling out; missing values found")

	// Standardizing makes it hard to detect values that are perfectly collinear with the absvars
	// in which case they should be 0.00 but they end up as e.g. 1e-16
	// EG: reghdfe price ibn.foreign , absorb(foreign)

	// This will edit to zero entire columns where *ALL* values are very close to zero
	if (timeit) timer_on(64)
	vars = cols(varlist) > 1 ? varlist : tokens(varlist)
	if (cols(vars)!=cols(y)) vars ="variable #" :+ strofreal(1..cols(y))
	collinear_tol = min(( 1e-6 , tolerance / 10))

	kept2 = (diagonal(cross(y, y))' :/ kept2) :> (collinear_tol)
	if (first_is_depvar & kept2[1]==0) {
		kept2[1] = 1
		if (verbose > -1) printf("{txt}warning: %s might be perfectly explained by fixed effects (tol =%3.1e)\n", vars[1], collinear_tol)
	}
	needs_zeroing = `selectindex'(!kept2)
	if (cols(needs_zeroing)) {
		y[., needs_zeroing] = J(rows(y), cols(needs_zeroing), 0)
		for (i=1; i<=cols(vars); i++) {
			if (!kept2[i] & verbose>-1 & (i > 1 | !first_is_depvar)) {
				printf("{txt}note: {res}%s{txt} is probably collinear with the fixed effects (all partialled-out values are close to zero; tol =%3.1e)\n", vars[i], collinear_tol)
			}
		}
	}

	kept = kept, kept2
	if (timeit) timer_off(64)
}


`Variables' FixedEffects::project_one_fe(`Variables' y, `Integer' g)
{
	`Factor'                f
	`Boolean'               store_these_alphas
	`Matrix'                alphas, proj_y

	// Cons+K+W, Cons+K, K+W, K, Cons+W, Cons = 6 variants

	f = factors[g]
	store_these_alphas = storing_alphas & save_fe[g]
	if (store_these_alphas) assert(cols(y)==1)

	if (num_slopes[g]==0) {
		if (store_these_alphas) {
			alphas = panelmean(f.sort(y), f)
			asarray(factors[g].extra, "tmp_alphas", alphas)
			return(alphas[f.levels, .])
		}
		else {
			if (cols(y)==1 & f.num_levels > 1) {
				return(panelmean(f.sort(y), f)[f.levels])
			}
			else {
				return(panelmean(f.sort(y), f)[f.levels, .])
			}
		}
	}
	else {
		// This includes both cases, with and w/out intercept (## and #)
		if (store_these_alphas) {
			alphas = J(f.num_levels, intercepts[g] + num_slopes[g], .)
			proj_y = panelsolve_invsym(f.sort(y), f, intercepts[g], alphas)
			asarray(factors[g].extra, "tmp_alphas", alphas)
			return(proj_y)
		}
		else {
			return(panelsolve_invsym(f.sort(y), f, intercepts[g]))
		}
	}
}


`Void' FixedEffects::estimate_dof()
{
	`Boolean'               has_int
	`Integer'               g, h                        // index FEs (1..G)
	`Integer'               num_intercepts              // Number of absvars with an intercept term
	`Integer'               i_cluster, i_intercept, j_intercept
	`Integer'               i                           // index 1..G_extended
	`Integer'               j
	`Integer'               bg_verbose                  // verbose level when calling BipartiteGraph()
	`Integer'               m                           // Mobility groups between a specific pair of FEs
	`RowVector'             SubGs
	`RowVector'             offsets, idx, zeros, results
	`Matrix'                tmp
	`Variables'             data
	`DataCol'               cluster_data
	`String'                absvar, clustervar
	`Factor'                F
	`BipartiteGraph'        BG
	`Integer'               pair_count
	
	if (verbose > 0) printf("\n{txt}## Estimating degrees-of-freedom absorbed by the fixed effects\n\n")

	// Count all FE intercepts and slopes
	SubGs = intercepts + num_slopes
	G_extended = sum(SubGs)
	num_intercepts = sum(intercepts)
	offsets = runningsum(SubGs) - SubGs :+ 1 // start of each FE within the extended list
	idx = `selectindex'(intercepts) // Select all FEs with intercepts
	if (verbose > 0) printf("{txt}   - there are %f fixed intercepts and slopes in the %f absvars\n", G_extended, G)

	// Initialize result vectors and scalars
	doflist_M_is_exact = J(1, G_extended, 0)
	doflist_M_is_nested = J(1, G_extended, 0)
	df_a_nested = 0

	// (1) M will hold the redundant coefs for each extended absvar (G_extended, not G)
	doflist_M = J(1, G_extended, 0)
	assert(0 <= num_clusters & num_clusters <= 10)
	if (num_clusters > 0 & anyof(dofadjustments, "clusters")) {

		// (2) (Intercept-Only) Look for absvars that are clustervars
		for (i_intercept=1; i_intercept<=length(idx); i_intercept++) {
			g = idx[i_intercept]
			i = offsets[g]
			absvar = invtokens(tokens(ivars[g]), "#")
			if (anyof(clustervars, absvar)) {
				doflist_M[i] = factors[g].num_levels
				df_a_nested = df_a_nested + doflist_M[i]
				doflist_M_is_exact[i] = doflist_M_is_nested[i] = 1
				idx[i_intercept] = 0
				if (verbose > 0) printf("{txt} - categorical variable {res}%s{txt} is also a cluster variable, so it doesn't reduce DoF\n", absvar)
			}
		}
		idx = select(idx, idx)

		// (3) (Intercept-Only) Look for absvars that are nested within a clustervar
		for (i_cluster=1; i_cluster<= num_clusters; i_cluster++) {
			cluster_data = .
			if (!length(idx)) break // no more absvars to process
			for (i_intercept=1; i_intercept<=length(idx); i_intercept++) {

				g = idx[i_intercept]
				i = offsets[g]
				absvar = invtokens(tokens(ivars[g]), "#")
				clustervar = clustervars[i_cluster]
				if (doflist_M_is_exact[i]) continue // nothing to do

				if (cluster_data == .) {
					if (strpos(clustervar, "#")) {
						clustervar = subinstr(clustervars[i_cluster], "#", " ", .)
						F = factor(clustervar, sample, ., "", 0, 0, ., 0)
						cluster_data = F.levels
						F = Factor() // clear
					}
					else {
						cluster_data = __fload_data(clustervar, sample, 0)
					}
				}

				if (factors[g].nested_within(cluster_data)) {
					doflist_M[i] = factors[g].num_levels
					doflist_M_is_exact[i] = doflist_M_is_nested[i] = 1
					df_a_nested = df_a_nested + doflist_M[i]
					idx[i_intercept] = 0
					if (verbose > 0) printf("{txt} - categorical variable {res}%s{txt} is nested within a cluster variable, so it doesn't reduce DoF\n", absvar)
				}
			}
			idx = select(idx, idx)
		}
		cluster_data = . // save memory
	} // end of the two cluster checks (absvar is clustervar; absvar is nested within clustervar)


	// (4) (Intercept-Only) Every intercept but the first has at least one redundant coef.
	if (length(idx) > 1) {
		if (verbose > 0) printf("{txt}   - there is at least one redundant coef. for every set of FE intercepts after the first one\n")
		doflist_M[offsets[idx[2..length(idx)]]] = J(1, length(idx)-1, 1) // Set DoF loss of all intercepts but the first one to 1
	}

	// (5) (Intercept-only) Mobility group algorithm
	// Excluding those already solved, the first absvar is exact

	if (length(idx)) {
		i = idx[1]
		doflist_M_is_exact[i] = 1
	}

	// Compute number of dijsoint subgraphs / mobility groups for each pair of remaining FEs
	if (anyof(dofadjustments, "firstpair") | anyof(dofadjustments, "pairwise")) {
		BG = BipartiteGraph()
		bg_verbose = max(( verbose - 1 , 0 ))
		pair_count = 0

		for (i_intercept=1; i_intercept<=length(idx)-1; i_intercept++) {
			for (j_intercept=i_intercept+1; j_intercept<=length(idx); j_intercept++) {
				g = idx[i_intercept]
				h = idx[j_intercept]
				i = offsets[h]
				BG.init(&factors[g], &factors[h], bg_verbose)
				m = BG.init_zigzag()
				++pair_count
				if (verbose > 0) printf("{txt}   - mobility groups between FE intercepts #%f and #%f: {res}%f\n", g, h, m)
				doflist_M[i] = max(( doflist_M[i] , m ))
				if (j_intercept==2) doflist_M_is_exact[i] = 1
				if (pair_count & anyof(dofadjustments, "firstpair")) break
			}
			if (pair_count & anyof(dofadjustments, "firstpair")) break
		}
		BG = BipartiteGraph() // clear
	}
	// TODO: add group3hdfe

	// (6) See if cvars are zero (w/out intercept) or just constant (w/intercept)
	if (anyof(dofadjustments, "continuous")) {
		for (i=g=1; g<=G; g++) {
			// If model has intercept, redundant cvars are those that are CONSTANT
			// Without intercept, a cvar has to be zero within a FE for it to be redundant
			// Since S.fes[g].x are already demeaned IF they have intercept, we don't have to worry about the two cases
			has_int = intercepts[g]
			if (has_int) i++
			if (!num_slopes[g]) continue

			data = asarray(factors[g].extra, "x")
			assert(num_slopes[g]==cols(data))
			results = J(1, cols(data), 0)
			// float(1.1) - 1 == 2.384e-08 , so let's pick something bigger, 1e-6
			zeros = J(1, cols(data), 1e-6)
			// This can be speed up by moving the -if- outside the -for-
			for (j = 1; j <= factors[g].num_levels; j++) {
				tmp = colminmax(panelsubmatrix(data, j, factors[g].info))
				if (has_int) {
					results = results + ((tmp[2, .] - tmp[1, .]) :<= zeros)
				}
				else {
					results = results + (colsum(abs(tmp)) :<= zeros)
				}
			}
			data = .
			if (sum(results)) {
				if (has_int  & verbose) printf("{txt}   - the slopes in the FE #%f are constant for {res}%f{txt} levels, which don't reduce DoF\n", g, sum(results))
				if (!has_int & verbose) printf("{txt}   - the slopes in the FE #%f are zero for {res}%f{txt} levels, which don't reduce DoF\n", g, sum(results))
				doflist_M[i..i+num_slopes[g]-1] = results
			}
			i = i + num_slopes[g]
		}
	}

	// Store results (besides doflist_..., etc.)
	doflist_K = J(1, G_extended, .)
	for (g=1; g<=G; g++) {
		i = offsets[g]
		j = g==G ? G_extended : offsets[g+1]
		doflist_K[i..j] = J(1, j-i+1, factors[g].num_levels)
	}
	df_a_initial = sum(doflist_K)
	df_a_redundant = sum(doflist_M)
	df_a = df_a_initial - df_a_redundant
}



`Void' FixedEffects::prune_1core()
{
	// Note that we can't prune degree-2 nodes, or the graph stops being bipartite
	`Integer'               i, j, g
	`Vector'                subgraph_id
	
	`Vector'                idx
	`RowVector'             i_prune

	// For now; too costly to use prune for G=3 and higher
	// (unless there are *a lot* of degree-1 vertices)
	if (G!=2) return //assert_msg(G==2, "G==2") // bugbug remove?

	// Abort if the user set HDFE.prune = 0
	if (!prune) return

	// Pick two factors, and check if we really benefit from pruning
	prune = 0
	i_prune = J(1, 2, 0)
	for (g=i=1; g<=2; g++) {
		//if (intercepts[g] & !num_slopes[g] & factors[g].num_levels>100) {
		if (intercepts[g] & !num_slopes[g]) {
			i_prune[i++] = g // increments at the end
			if (i > 2) { // success!
				prune = 1
				break
			}
		}
	}

	if (!prune) return

	// for speed, the factor with more levels goes first
	i = i_prune[1]
	j = i_prune[2]
	//if (factors[i].num_levels < factors[j].num_levels) swap(i, j) // bugbug uncomment it!
	prune_g1 = i
	prune_g2 = j

	bg = BipartiteGraph()
	bg.init(&factors[prune_g1], &factors[prune_g2], verbose)
	(void) bg.init_zigzag(1) // 1 => save subgraphs into bg.subgraph_id
	bg.compute_cores()
	bg.prune_1core(weight)
	num_pruned = bg.N_drop
}

// bugbug fix or remove this fn altogether
`Void' FixedEffects::_expand_1core(`Variables' y)
{
	y = bg.expand_1core(y)
}


`Void' FixedEffects::save_touse(| `Varname' touse, `Boolean' replace)
{
	`Integer'               idx
	`Vector'                mask

	// Set default arguments
	if (args()<1 | touse=="") {
		assert(tousevar != "")
		touse = tousevar
	}
	// Note that args()==0 implies replace==1 (else how would you find the name)
	if (args()==0) replace = 1
	if (args()==1 | replace==.) replace = 0

	if (verbose > 0) printf("\n{txt}## Saving e(sample)\n")

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
}


`Void' FixedEffects::save_variable(`Varname' varname,
								   `Variables' data,
								 | `Varlist' varlabel)
{
	`RowVector'               idx
	`Integer'               i
	idx = st_addvar("double", tokens(varname))
	st_store(sample, idx, data)
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

	st_global("e(absvars)", invtokens(absvars))
	text = invtokens(extended_absvars)
	text = subinstr(text, "1.", "")
	st_global("e(extended_absvars)", text)

	// Absorbed degrees-of-freedom table
	table = (doflist_K \ doflist_M \ (doflist_K-doflist_M) \ !doflist_M_is_exact \ doflist_M_is_nested)'
	rowstripe = extended_absvars'
	rowstripe = J(rows(table), 1, "") , extended_absvars' // add equation col
	colstripe = "Categories" \ "Redundant" \ "Num Coefs" \ "Exact?" \ "Nested?" // colstripe cannot have dots on Stata 12 or earlier
	colstripe = J(cols(table), 1, "") , colstripe // add equation col
	st_matrix("e(dof_table)", table)
	st_matrixrowstripe("e(dof_table)", rowstripe)
	st_matrixcolstripe("e(dof_table)", colstripe)

	st_numscalar("e(ic)", iteration_count)
	st_numscalar("e(drop_singletons)", drop_singletons)
	st_numscalar("e(num_singletons)", num_singletons)
	st_numscalar("e(N_full)", st_numscalar("e(N)") + num_singletons)

	// Prune specific
	if (prune==1) {
		st_numscalar("e(pruned)", 1)
		st_numscalar("e(num_pruned)", num_pruned)
	}

	if (!missing(finite_condition)) st_numscalar("e(finite_condition)", finite_condition)
}


`Void' FixedEffects::post()
{
	`String'        text
	`Integer'       i

	post_footnote()

	// ---- constants -------------------------------------------------------

	st_global("e(predict)", "reghdfe_p")
	st_global("e(estat_cmd)", "reghdfe_estat")
	st_global("e(footnote)", "reghdfe_footnote")
	//st_global("e(marginsok)", "")
	st_global("e(marginsnotok)", "Residuals SCore")
	st_numscalar("e(df_m)", df_m)


	assert(title != "")
	text = sprintf("HDFE %s", title)
	st_global("e(title)", text)
	
	text = sprintf("Absorbing %g HDFE %s", G, plural(G, "group"))
	st_global("e(title2)", text)
	
	st_global("e(model)", model)
	st_global("e(cmdline)", cmdline)

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
	
	if (!missing(N_clust)) {
		st_numscalar("e(N_clust)", N_clust)
		for (i=1; i<=num_clusters; i++) {
			text = sprintf("e(N_clust%g)", i)
			st_numscalar(text, N_clust_list[i])
		}
		text = "Statistics robust to heteroskedasticity"
		st_global("e(title3)", text)
	}

	if (!missing(sumweights)) st_numscalar("e(sumweights)", sumweights)

	st_numscalar("e(report_constant)", compute_constant & report_constant)


	// ---- .options properties ---------------------------------------------

	st_global("e(depvar)", depvar)
	st_global("e(indepvars)", invtokens(indepvars))

	if (!missing(N_clust)) {
		st_numscalar("e(N_clustervars)", num_clusters)
		st_global("e(clustvar)", invtokens(clustervars))
		for (i=1; i<=num_clusters; i++) {
			text = sprintf("e(clustvar%g)", i)
			st_global(text, clustervars[i])
		}
	}

	if (residuals != "") {
		st_global("e(resid)", residuals)
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


// --------------------------------------------------------------------------
// Recreate HDFE object
// --------------------------------------------------------------------------
`FixedEffects' FixedEffects::reload(`Boolean' copy)
{
	`FixedEffects' ans
	assert(copy==0 | copy==1)
	
	// Trim down current object as much as possible
	// this. is optional but useful for clarity
	if (copy==0) {
		this.factors = Factor()
		this.sample = .
		this.bg = BipartiteGraph()
		this.pruned_weight = .
		this.rre_varname = .
		this.rre_true_residual = .
	}

	// Initialize new object
	ans = fixed_effects(this.absorb, this.tousevar, this.weight_type, this.weight_var, this.drop_singletons, this.verbose)

	// Fill out new object with values of current one
	ans.depvar = this.depvar
	ans.indepvars = this.indepvars
	ans.varlist = this.varlist
	ans.model = this.model
	ans.vcetype = this.vcetype
	ans.num_clusters = this.num_clusters
	ans.clustervars = this.clustervars
	ans.base_clustervars = this.base_clustervars
	ans.vceextra = this.vceextra
	ans.summarize_stats = this.summarize_stats
	ans.summarize_quietly = this.summarize_quietly
	ans.notes = this.notes
	ans.store_sample = this.store_sample
	ans.timeit = this.timeit
	ans.compact = this.compact
	ans.poolsize = this.poolsize
	ans.diopts = this.diopts

	ans.fullindepvars = this.fullindepvars
	ans.not_basevar = this.not_basevar

	ans.compute_constant = this.compute_constant
	ans.report_constant = this.report_constant
	ans.tolerance = this.tolerance
	ans.save_any_fe = this.save_any_fe

	ans.slope_method = this.slope_method
	ans.maxiter = this.maxiter
	ans.transform = this.transform
	ans.acceleration = this.acceleration
	ans.accel_start = this.accel_start
	ans.conlim = this.conlim
	ans.btol = this.btol
	ans.min_ok = this.min_ok
	ans.prune = this.prune
	ans.always_run_lsmr_preconditioner = this.always_run_lsmr_preconditioner

	return(ans)
}


// --------------------------------------------------------------------------
// Estimate finite condition number of the graph Laplacian
// --------------------------------------------------------------------------
`Void' FixedEffects::estimate_cond()
{
	`Matrix'                D1, D2, L
	`Vector'                lambda
	`RowVector'             tmp
	`Factor'                F12

	if (finite_condition!=-1) return

	if (verbose > 0) printf("\n{txt}## Computing finite condition number of the Laplacian\n\n")

	if (verbose > 0) printf("{txt}   - Constructing vectors of levels\n")
	F12 = join_factors(factors[1], factors[2], ., ., 1)
	
	// Non-sparse (lots of memory usage!)
	if (verbose > 0) printf("{txt}   - Constructing design matrices\n")
	D1 = designmatrix(F12.keys[., 1])
	D2 = designmatrix(F12.keys[., 2])
	assert_msg(rows(D1)<=2000, "System is too big!")
	assert_msg(rows(D2)<=2000, "System is too big!")

	if (verbose > 0) printf("{txt}   - Constructing graph Laplacian\n")
	L =   D1' * D1 , - D1' * D2 \
		- D2' * D1 ,   D2' * D2
	if (verbose > 0) printf("{txt}   - L is %g x %g \n", rows(L), rows(L))
		
	if (verbose > 0) printf("{txt}   - Computing eigenvalues\n")
	assert_msg(rows(L)<=2000, "System is too big!")
	eigensystem(L, ., lambda=.)
	lambda = Re(lambda')

	if (verbose > 0) printf("{txt}   - Selecting positive eigenvalues\n")
	lambda = edittozerotol(lambda, 1e-8)
	tmp = select(lambda,  edittozero(lambda, 1))
	tmp = minmax(tmp)
	finite_condition = tmp[2] / tmp[1]

	if (verbose > 0) printf("{txt}   - Finite condition number: {res}%s{txt}\n", strofreal(finite_condition))
}


`Real' FixedEffects::lsmr_norm(`Matrix' x)
{
	assert(cols(x)==1 | rows(x)==1)
	// BUGBUG: what if we have a corner case where there are as many obs as params?
	if (has_weights & cols(x)==1 & rows(x)==rows(weight)) {
		return(sqrt(quadcross(x, weight, x)))
	}
	else if (cols(x)==1) {
		return(sqrt(quadcross(x, x)))
	}
	else {
		return(sqrt(quadcross(x', x')))
	}
}


// Ax: given the coefs 'x', return the projection 'Ax'
`Vector' FixedEffects::lsmr_A_mult(`Vector' x)
{
	`Integer' g, k, idx_start, idx_end, i
	`Vector' ans
	`FactorPointer'         pf

	ans = J(N, 1, 0)
	idx_start = 1

	for (g=1; g<=G; g++) {
		pf = &(factors[g])
		k = (*pf).num_levels

		if (intercepts[g]) {
			idx_end = idx_start + k - 1
			ans = ans + (x[|idx_start, 1 \ idx_end , 1 |] :* asarray((*pf).extra, "precond_intercept") )[(*pf).levels, .]
			idx_start = idx_end + 1
		}

		if (num_slopes[g]) {
			for (i=1; i<=num_slopes[g]; i++) {
				idx_end = idx_start + k - 1
				ans = ans + x[|idx_start, 1 \ idx_end , 1 |][(*pf).levels] :* asarray((*pf).extra, "precond_slopes")
				idx_start = idx_end + 1
			}
		}

	}
	//assert(!missing(ans))
	return(ans)
}


// A'x: Compute the FEs and store them in a big stacked vector
`Vector' FixedEffects::lsmr_At_mult(`Vector' x)
{
	`Integer' m, g, i, idx_start, idx_end, k
	`Vector' ans
	`FactorPointer'         pf
	`Vector' alphas
	`Matrix' tmp_alphas

	alphas = J(M, 1, .)
	idx_start = 1

	for (g=1; g<=G; g++) {
		pf = &(factors[g])
		k = (*pf).num_levels

		if (intercepts[g]) {
			idx_end = idx_start + k - 1
			if (has_weights) {
				alphas[| idx_start , 1 \ idx_end , 1 |] = `panelsum'((*pf).sort(x :* weight), (*pf).info) :* asarray((*pf).extra, "precond_intercept")
			}
			else {
				alphas[| idx_start , 1 \ idx_end , 1 |] = `panelsum'((*pf).sort(x), (*pf).info) :* asarray((*pf).extra, "precond_intercept")
			}
			idx_start = idx_end + 1
		}

		if (num_slopes[g]) {
			idx_end = idx_start + k * num_slopes[g] - 1
			if (has_weights) {
				tmp_alphas = `panelsum'((*pf).sort(x :* weight :* asarray((*pf).extra, "precond_slopes")), (*pf).info)
			}
			else {
				tmp_alphas = `panelsum'((*pf).sort(x :* asarray((*pf).extra, "precond_slopes")), (*pf).info)
			}
			alphas[| idx_start , 1 \ idx_end , 1 |] = vec(tmp_alphas)
			idx_start = idx_end + 1
		}
	}
	//assert(!missing(alphas))
	return(alphas)
}

end
