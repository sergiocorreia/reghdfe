mata:
mata set matastrict on

//  Parse absvars and initialize the almost empty MapProblem struct
`Problem' function map_init()
{
	`Integer'		g, G, num_slopes, has_intercept, i, H, j
	`Problem' 		S
	`Boolean'		auto_target // Automatically assign target names to all FEs
	`Varname'		basetarget
	`Varlist'		target, original_absvars, extended_absvars
	`String'		equation_d
	`Boolean'		equation_d_valid
	pointer(`FE') 	fe

	S.weightvar = S.weighttype = S.weights = ""
	S.verbose = 0
	S.transform = "symmetric_kaczmarz" // cimmino ?
	S.acceleration = "conjugate_gradient"
	S.tolerance = 1e-8
	S.maxiterations = 1e4
	S.accel_start = 6
	S.poolsize = 10

	// If clustering by timevar or panelvar and VCE is HAC, we CANNOT touch the clustervars to create compact ids!
	S.timevar = ""
	S.panelvar = ""
	S.vce_is_hac = 0

	// Specific to Aitken:
	S.accel_freq = 3
	S.pause_length = 20
	S.bad_loop_threshold = 1
	S.stuck_threshold = 5e-3
	S.N = .

	S.groupvar = "" // Initialize as empty
	S.sortedby = "" // Initialize as empty (prevents bugs if we change the dataset before map_precompute)


	S.keepsingletons = 0
	S.G = G = st_numscalar("r(G)")
	S.C = 0
	S.clustervars = S.clustervars_original = J(0,0,"")
	S.fes = FixedEffect(G)
	S.will_save_fe = 0
	auto_target = st_numscalar("r(savefe)")
	assert(auto_target==1 | auto_target==0)
	if (auto_target) stata(sprintf("cap drop __hdfe*__*"))

	original_absvars = extended_absvars = J(1, G, "")

	for (g=1; g<=G; g++) {
		fe = &(S.fes[g])
		// recall a->b is the same as (*a).b
		num_slopes = st_numscalar(sprintf("r(num_slopes%f)",g))
		has_intercept = st_numscalar(sprintf("r(has_intercept%1.0f)",g))

		fe->order = g
		fe->num_slopes = num_slopes
		fe->has_intercept = has_intercept
		fe->varlabel = st_global(sprintf("r(varlabel%f)",g))
		fe->ivars = tokens(st_global( sprintf("r(ivars%f)",g) ))
		fe->cvars = tokens(st_global( sprintf("r(cvars%f)",g) ))
		fe->idvarname = sprintf("__ID%f__", g)
		fe->levels = .
		fe->target = J(0,0,"")
		fe->is_clustervar = 0
		fe->in_clustervar = 0
		fe->nesting_clustervar = .

		extended_absvars[g] = original_absvars[g] = invtokens(fe->ivars, "#")
		if (num_slopes>0) {
			original_absvars[g] = original_absvars[g] + (has_intercept ? "##c." : "#c.")
			original_absvars[g] = original_absvars[g] + (num_slopes==1 ? (fe->cvars) : "("+invtokens(fe->cvars)+")")
			
			extended_absvars[g] = (has_intercept ? extended_absvars[g] + " " : "") + invtokens(extended_absvars[g] + "#c." :+ (fe->cvars))
		}
			
		basetarget = st_global(sprintf("r(target%f)",g))
		if (basetarget=="" & auto_target) basetarget = sprintf("__hdfe%f__", g)
		if (basetarget!="") {
			S.will_save_fe = 1
			target = J(1, has_intercept + num_slopes, basetarget)
			if (has_intercept) stata(sprintf("confirm new variable %s", target[1]))
			for (i=1+has_intercept; i<=length(target); i++) {
				target[i] = target[i] + sprintf("_Slope%f", i-has_intercept)
				stata(sprintf("confirm new variable %s", target[i]))
			}
			fe->target = target
		}
	}

	equation_d_valid = 1
	equation_d = ""
	for (g=1; g<=S.G; g++) {
		H = S.fes[g].has_intercept + S.fes[g].num_slopes
		if (length(S.fes[g].target)==0) {
			equation_d_valid = 0
			break
		}
		assert(length(S.fes[g].target==H))

		j = 0
		for (i=1; i<=H;i++) {
			equation_d = equation_d + sprintf("%s%s", equation_d=="" ? "" : " + ", S.fes[g].target[i])
			if (i>1 | !S.fes[g].has_intercept) j++ // j is the cvar counter
			if (j>0) equation_d = equation_d + sprintf(" * %s", S.fes[g].cvars[j])
		}
	}
	if (equation_d_valid) st_global("r(equation_d)", invtokens(equation_d) )
	
	st_numscalar("r(will_save_fe)", S.will_save_fe)
	st_global("r(original_absvars)", invtokens(original_absvars) )
	st_global("r(extended_absvars)", invtokens(extended_absvars) )
	return(S)
}

void function map_init_clustervars(`Problem' S, `String' clustervars) {
	S.clustervars = S.clustervars_original = tokens(clustervars)
	S.C = length(S.clustervars)
}

void function map_init_weights(`Problem' S, `Varname' weightvar, `String' weighttype) {
	assert_msg(weightvar!="" & weighttype!="", "map_init_weights() requires weight var and type")
	stata(sprintf("confirm numeric variable %s, exact", weightvar))
	assert_msg(anyof(("fweight", "pweight", "aweight"), weighttype), "wrong weight type")
	S.weightvar = weightvar
	S.weighttype = weighttype
	S.weights = sprintf("[%s=%s]", weighttype, weightvar)
}

void function map_init_keepvars(`Problem' S, `Varname' keepvars) {
	if (keepvars!="") stata(sprintf("confirm numeric variable %s, exact", keepvars))
	S.keepvars = tokens(keepvars)
}

void function map_init_transform(`Problem' S, `String' transform) {
	transform = strlower(transform)
	// Convert abbreviations
	if (strpos("cimmino", transform)==1) transform = "cimmino"
	if (strpos("kaczmarz", transform)==1) transform = "kaczmarz"
	if (strpos("symmetric_kaczmarz", transform)==1) transform = "symmetric_kaczmarz"
	if (strpos("rand_kaczmarz", transform)==1) transform = "random_kaczmarz" // experimental
	if (strpos("random_kaczmarz", transform)==1) transform = "random_kaczmarz" // experimental
	assert_msg(anyof(("cimmino","kaczmarz","symmetric_kaczmarz", "random_kaczmarz"),transform), "invalid transform")
	S.transform = transform
}

void function map_init_verbose(`Problem' S, `Integer' verbose) {
	assert_msg(round(verbose)==verbose, "verbose must be an integer")
	assert_msg(0<=verbose & verbose<=5, "verbose must be between 0 and 5")
	S.verbose = verbose
}

void function map_init_timeit(`Problem' S, `Integer' timeit) {
	assert_msg(timeit==0 | timeit==1, "timeit must be 0 or 1")
	S.timeit = timeit
}

void function map_init_poolsize(`Problem' S, `Integer' poolsize) {
	assert_msg(round(poolsize)==poolsize & poolsize>0, "poolsize must be a positive integer")
	S.poolsize = poolsize
}

void function map_init_panelvar(`Problem' S, `Varname' panelvar) {
	S.panelvar = panelvar
}

void function map_init_timevar(`Problem' S, `Varname' timevar) {
	S.timevar = timevar
}

void function map_init_groupvar(`Problem' S, `Varname' groupvar) {
	S.groupvar = groupvar
}

void function map_init_acceleration(`Problem' S, `String' acceleration) {
	acceleration = strlower(acceleration)
	// Convert abbreviations
	if (strpos("conjugate_gradient", acceleration)==1 | acceleration=="cg") acceleration = "conjugate_gradient"
	if (strpos("steepest_descent", acceleration)==1 | acceleration=="sd") acceleration = "steepest_descent"
	if (strpos("aitken", acceleration)==1) acceleration = "aitken"
	if (strpos("hybrid", acceleration)==1) acceleration = "hybrid" // experimental
	if (acceleration=="no" | acceleration=="none" | acceleration=="off") acceleration = "none"
	assert_msg(anyof(("conjugate_gradient","steepest_descent", "aitken", "none", "hybrid"),acceleration), "invalid acceleration")
	S.acceleration = acceleration
}

void function map_init_tolerance(`Problem' S, `Real' tolerance) {
	assert_msg(1e-16<=tolerance & tolerance<1, "tolerance must be in the range [1e-16, 1).")
	S.tolerance = tolerance
}

void function map_init_maxiterations(`Problem' S, `Integer' maxiterations) {
	assert_msg(round(maxiterations)==maxiterations, "maxiterations must be an integer")
	assert_msg(maxiterations>0, "maxiterations must be positive")
	S.maxiterations = maxiterations
}

void function map_init_keepsingletons(`Problem' S, `Boolean' keepsingletons) {
	assert_msg(keepsingletons==0 | keepsingletons==1, "keepsingletons must be 0 or 1")
	S.keepsingletons = keepsingletons
}

void function map_init_vce_is_hac(`Problem' S, `Boolean' vce_is_hac) {
	assert_msg(vce_is_hac==0 | vce_is_hac==1, "vce_is_hac must be 0 or 1")
	S.vce_is_hac = vce_is_hac
}

end
