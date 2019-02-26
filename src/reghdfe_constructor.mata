// --------------------------------------------------------------------------
// FixedEffects constructor (also precomputes factors)
// --------------------------------------------------------------------------

mata:

`FixedEffects' fixed_effects(`Varlist' absvars,
						   | `Varname' touse,
							 `String' weighttype,
							 `Varname' weightvar,
							 `Boolean' drop_singletons,
							 `Boolean' verbose)
{
	`FixedEffects'          S
	`Varname'               absvar, cvars
	`Integer'               i, j, g, gg, remaining
	`Vector'                idx
	`Integer'               spaces
	`Integer'               num_singletons_i
	`Variables'             cvar_data
	`FactorPointer'         pf

	// Set default value of arguments
	if (args()<2) touse = ""
	if (args()<3) weighttype = ""
	if (args()<4) weightvar = ""
	if (args()<5 | drop_singletons==.) drop_singletons = 1
	if (args()<6 | verbose==.) verbose = 0
	
	S = FixedEffects()
	S.verbose = verbose
	S.drop_singletons = drop_singletons

	// Parse absvars
	if (S.verbose > 0) printf("\n{txt}## Parsing absvars and HDFE options\n")
	
	if (touse == "") touse = st_tempname()
	st_global("reghdfe_touse", touse)
	stata(`"reghdfe_parse "' + absvars)
	S.sample = `selectindex'(st_data(., touse))
	S.tousevar = touse // useful if later on we want to clone the HDFE object
	st_global("reghdfe_touse", "")

	if (st_global("s(residuals)") != "") S.residuals = st_global("s(residuals)")
	if (st_global("s(verbose)")!="") S.verbose = verbose = strtoreal(st_global("s(verbose)"))
	if (st_global("s(drop_singletons)")!="") S.drop_singletons = drop_singletons = strtoreal(st_global("s(drop_singletons)"))
	assert(S.verbose < .)
	assert(S.drop_singletons==0 | S.drop_singletons==1)

	if (S.verbose > 0) stata("sreturn list")
	S.G = strtoreal(st_global("s(G)"))
	S.absorb = absvars // useful if later on we want to clone the HDFE object
	S.absvars = tokens(st_global("s(absvars)"))
	S.has_intercept = strtoreal(st_global("s(has_intercept)"))
	S.save_any_fe = strtoreal(st_global("s(save_any_fe)"))
	S.save_all_fe = strtoreal(st_global("s(save_all_fe)"))
	S.ivars = tokens(st_global("s(ivars)"))
	S.cvars = tokens(st_global("s(cvars)"))
	S.targets = strtrim(tokens(st_global("s(targets)")))
	S.intercepts = strtoreal(tokens(st_global("s(intercepts)")))
	S.num_slopes = strtoreal(tokens(st_global("s(num_slopes)")))
	S.save_fe = S.targets :!= ""
	S.report_constant = strtoreal(st_global("s(report_constant)"))
	S.always_run_lsmr_preconditioner = strtoreal(st_global("s(precondition)"))

	// Ensure that S.report_constant and S.has_intercept are 0/1
	assert(anyof((0,1), S.has_intercept))
	assert(anyof((0,1), S.report_constant))
	S.compute_constant = S.has_intercept & S.report_constant

	if (st_global("s(tolerance)") != "") S.tolerance = strtoreal(st_global("s(tolerance)"))
	if (st_global("s(maxiter)") != "") S.maxiter = strtoreal(st_global("s(maxiter)"))
	if (st_global("s(prune)") != "") S.prune = strtoreal(st_global("s(prune)"))
	if (st_global("s(transform)") != "") S.transform = st_global("s(transform)")
	if (st_global("s(acceleration)") != "") S.acceleration = st_global("s(acceleration)")

	// Override LSMR if G=1
	if (S.G==1 & S.acceleration=="lsmr") S.acceleration = "conjugate_gradient"

	S.dofadjustments = tokens(st_global("s(dofadjustments)"))
	S.groupvar = st_global("s(groupvar)")
	if (st_global("s(finite_condition)")=="1") S.finite_condition = -1 // signal to compute it
	S.compute_rre = (st_global("s(compute_rre)")=="1")
	if (S.compute_rre) S.rre_varname = st_global("s(rre)")
	
	S.poolsize = strtoreal(st_global("s(poolsize)"))

	if (S.verbose > -1 & !S.has_intercept) printf("{txt}(warning: no intercepts terms in absorb(); regression lacks constant term)\n")

	S.extended_absvars = tokens(st_global("s(extended_absvars)"))
	S.tss = .

	assert(1<=S.G)
	if (S.G>10) printf("{txt}(warning: absorbing %2.0f dimensions of fixed effects; check that you really want that)\n", S.G)
	assert(S.G == cols(S.ivars))
	assert(S.G == cols(S.cvars))
	assert(S.G == cols(S.targets))
	assert(S.G == cols(S.intercepts))
	assert(S.G == cols(S.num_slopes))

	// Fill out object
	S.G = cols(S.absvars)
	S.factors = Factor(S.G)

	assert_msg(anyof(("", "fweight", "pweight", "aweight", "iweight"), weighttype), "wrong weight type")
	S.weight_type = weighttype
	S.weight_var = weightvar

	S.num_singletons = 0
	if (drop_singletons) {
		num_singletons_i = 0
		if (weighttype=="fweight" | weighttype=="iweight") {
			S.weight = st_data(S.sample, weightvar) // just to use it in F.drop_singletons()
		}
	}


	// (1) create the factors and remove singletons
	remaining = S.G
	i = 0
	if (S.verbose > 0) {
		printf("\n{txt}## Initializing Mata object for %g fixed effects\n\n", S.G)
		spaces = max((0, max(strlen(S.absvars))-4))
		printf("{txt}   {c TLC}{hline 4}{c TT}{hline 3}{c TT}{hline 1}%s{hline 6}{c TT}{hline 6}{c TT}{hline 9}{c TT}{hline 11}{c TT}{hline 12}{c TT}{hline 9}{c TT}{hline 14}{c TRC}\n", "{hline 1}" * spaces)
		printf("{txt}   {c |}  i {c |} g {c |} %s Name {c |} Int? {c |} #Slopes {c |}    Obs.   {c |}   Levels   {c |} Sorted? {c |} #Drop Singl. {c |}\n", " " * spaces)
		printf("{txt}   {c LT}{hline 4}{c +}{hline 3}{c +}{hline 1}%s{hline 6}{c +}{hline 6}{c +}{hline 9}{c +}{hline 11}{c +}{hline 12}{c +}{hline 9}{c +}{hline 14}{c RT}\n", "{hline 1}" * spaces)
		displayflush()
	}

	while (remaining) {
		++i
		g = 1 + mod(i-1, S.G)
		absvar = S.absvars[g]
		
		if (S.verbose > 0) {
			printf("{txt}   {c |} %2.0f {c |} %1.0f {c |} {res}%s{txt} {c |} ", i, g, (spaces+5-strlen(absvar)) * " " + absvar)
			printf("{txt}{%s}%3s{txt}  {c |}    %1.0f    {c |}", S.intercepts[g] ? "txt" : "err", S.intercepts[g] ? "Yes" : "No", S.num_slopes[g])
			displayflush()
		}

		if (S.verbose > 0) {
			printf("{res}%10.0g{txt} {c |}", rows(S.sample))
			displayflush()
		}

		if (rows(S.sample) < 2) {
			if (S.verbose > 0) printf("\n")
			exit(error(2001))
		}

		if (i<=S.G) {
			if (S.ivars[g] == "_cons" & S.G == 1) {
				// Special case without any fixed effects

				S.factors[g] = Factor()
				pf = &(S.factors[g])
				(*pf).num_obs = (*pf).counts = rows(S.sample)
				(*pf).num_levels = 1
				//(*pf).levels = . // Not filled to save space
				(*pf).levels = J(rows(S.sample), 1, 1)
				(*pf).is_sorted = 1
				(*pf).method = "none"

				// The code below is equivalent but 3x slower
				// S.factors[g] = _factor(J(rows(S.sample),1,1), 1, ., "hash0", ., 1, ., 0)
			}
			else {
				// We don't need to save keys (or sort levels but that might change estimates of FEs)
				S.factors[g] = factor(S.ivars[g], S.sample, ., "", ., 1, ., 0)
			}
		}

		if (S.verbose > 0) {
			printf(" {res}%10.0g{txt} {c |} %7s {c |}", S.factors[g].num_levels, S.factors[g].is_sorted ? "Yes" : "No")
			displayflush()
		}
 
		if (drop_singletons) {
			
			if (weighttype=="fweight") {
				idx = S.factors[g].drop_singletons(S.weight)
			}
			else if (weighttype=="iweight") {
				idx = S.factors[g].drop_singletons(S.weight, 1) // zero_threshold==1
			}
			else {
				idx = S.factors[g].drop_singletons()
			}

			num_singletons_i = rows(idx)
			S.num_singletons = S.num_singletons + num_singletons_i
			if (S.verbose > 0) {
				printf(" %10.0g   {c |}", num_singletons_i)
				displayflush()
			}

			if (num_singletons_i==0) {
				--remaining
			}
			else {
				remaining = S.G - 1
				
				// sample[idx] = . // not allowed in Mata; instead, make 0 and then select()
				S.sample[idx] = J(rows(idx), 1, 0)
				S.sample = select(S.sample, S.sample)

				for (j=i-1; j>=max((1, i-remaining)); j--) {
					gg = 1 + mod(j-1, S.G)
					S.factors[gg].drop_obs(idx)
					if (S.verbose > 0) printf("{res} .")
				}
			}
		}
		else {
			if (S.verbose > 0) printf("      n/a     {c |}")
			--remaining
		}
		if (S.verbose > 0) printf("\n")
	}
	if (S.verbose > 0) {
		printf("{txt}   {c BLC}{hline 4}{c BT}{hline 3}{c BT}{hline 1}%s{hline 6}{c BT}{hline 6}{c BT}{hline 9}{c BT}{hline 11}{c BT}{hline 12}{c BT}{hline 9}{c BT}{hline 14}{c BRC}\n", "{hline 1}" * spaces)
	}

	if ( drop_singletons & S.num_singletons>0 & S.verbose>-1 | S.factors[1].num_obs<2) {
		if (weighttype=="iweight") {
			// PPML-specific
			printf(`"{txt}(dropped %s observations that are either {browse "http://scorreia.com/research/singletons.pdf":singletons} or {browse "http://scorreia.com/research/separation.pdf":separated} by a fixed effect)\n"', strofreal(S.num_singletons))
		}
		else {
			printf(`"{txt}(dropped %s {browse "http://scorreia.com/research/singletons.pdf":singleton observations})\n"', strofreal(S.num_singletons))
		}
	}

	if (S.factors[1].num_obs < 2) {
		exit(error(2001))
	}

	S.N = S.factors[1].num_obs // store number of obs.
	assert(S.N = S.factors[S.G].num_obs)
	assert(S.N > 1)


	// (2) run *.panelsetup() after the sample is defined
	if (S.verbose > 0) printf("\n{txt}## Initializing panelsetup() for each fixed effect\n\n")
	for (g=1; g<=S.G; g++) {
		absvar = S.absvars[g]
		if (S.verbose > 0) printf("{txt}   - panelsetup({res}%s{txt})\n", absvar)
		S.factors[g].panelsetup()
	}

	// (3) load cvars
	if (sum(S.num_slopes)) {
		if (S.verbose > 0) printf("\n{txt}## Loading slope variables\n\n")
		for (g=1; g<=S.G; g++) {
			cvars = tokens(S.cvars[g])
			if (S.num_slopes[g]) {
				// Load, standardize, sort by factor and store
				// Don't precompute aux objects (xmeans, inv_xx) as they depend on the weights
				// and will be computed on step (5)
				if (S.verbose > 0) printf("{txt}   - cvars({res}%s{txt})\n", invtokens(cvars))
				pf = &(S.factors[g])
				cvar_data = (*pf).sort(st_data(S.sample, cvars))
				asarray((*pf).extra, "x_stdevs", reghdfe_standardize(cvar_data))
				asarray((*pf).extra, "x", cvar_data)
			}
		}
		cvar_data = .
	}

	// (4) prune edges of degree-1
	// S.prune = 0 // bugbug
	if (S.prune) S.prune_1core()

	// (5) load weight
	S.load_weights(weighttype, weightvar, J(0,1,.), 1) // update S.has_weights, S.factors, etc.

	// Save "true" residuals for RRE
	if (S.compute_rre) {
		assert_msg(S.rre_varname != "")
		S.rre_true_residual = st_data(S.sample, S.rre_varname)
	}

	return(S)
}

end
