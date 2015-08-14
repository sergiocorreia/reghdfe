mata:
mata set matastrict on

void function map_solve(`Problem' S, `Varlist' vars, | `Varlist' newvars, `Boolean' restore_dta) {

	`Integer' i, j, h, Q
	`Group' y, residuals
	`Varlist'	chars
	real rowvector stdevs
	`FunctionPointer' transform, accelerate

	if (S.verbose>0) printf("{txt}{bf:mata: map_solve()}\n")
	assert_msg(S.N!=., "map_solve() needs to be run after map_precompute()")
	assert_msg(S.N==st_nobs(), "dataset cannot change after map_precompute()")
	S.storing_betas = 0

	if (S.verbose>0) printf("{txt} - Loading variables into Mata\n")
	vars = tokens(vars)
	Q = cols(vars)

	// Store chars var[name] that contain the original varname (e.g. L.var)
	chars = J(1, Q, "")
	for (i=1;i<=Q;i++) {
		chars[i] = st_global(sprintf("%s[name]", vars[i]))
	}

	// Optional arguments
	newvars = (args()<3 | newvars=="") ? vars : tokens(newvars)
	assert_msg(length(vars)==length(newvars), "map_solve error: newvars must have the same size as vars")
	if (args()<4) restore_dta = 0 // restore dataset (only for hdfe.ado)
	assert_msg(restore_dta==0 | restore_dta==1, "restore_dta must be 0 or 1")
	
	// Solver Warnings
	if (S.transform=="kaczmarz" & S.acceleration=="conjugate_gradient") {
		printf("{err}(WARNING: convergence is {bf:unlikely} with transform=kaczmarz and accel=CG)\n")
	}

	// Load transform pointer
	if (S.transform=="cimmino") transform = &transform_cimmino()
	if (S.transform=="kaczmarz") transform = &transform_kaczmarz()
	if (S.transform=="symmetric_kaczmarz") transform = &transform_sym_kaczmarz()
	if (S.transform=="random_kaczmarz") transform = &transform_rand_kaczmarz()

	// Pointer to acceleration routine
	if (S.acceleration=="none") accelerate = &accelerate_none()
	if (S.acceleration=="conjugate_gradient") accelerate = &accelerate_cg()
	if (S.acceleration=="steepest_descent") accelerate = &accelerate_sd()
	if (S.acceleration=="aitken") accelerate = &accelerate_aitken()
	if (S.acceleration=="hybrid") accelerate = &accelerate_hybrid()

	// Shortcut for trivial case (1 FE)
	if (S.G==1) accelerate = &accelerate_none()

	// Iterate on variables grouped by poolsize: subset = data[i..j]
	S.num_iters_max = 0
	if (restore_dta) residuals = J(S.N, 0, .)

	for (i=1;i<=Q;i=i+S.poolsize) {
		
		j = i + S.poolsize - 1
		if (j>Q) j = Q

		// Load data
		y = st_data(., vars[i..j])
		
		// Drop loaded vars as quickly as possible (also they might not even be -double-)
		st_dropvar(vars[i..j])
		
		// Standardize variables
		if (S.verbose>0) printf("{txt} - Standardizing variables\n")
		stdevs = J(1, cols(y), .)
		for (h=1; h<=cols(y); h++) {
			stdevs[h] = max((  sqrt(quadvariance(y[., h])) , sqrt(epsilon(1)) ))
		}
		y = y :/ stdevs

		// Solve!		
		if (S.verbose>0) printf("{txt} - Solving problem (acceleration={res}%s{txt}, transform={res}%s{txt} tol={res}%-1.0e{txt} poolsize={res}%f{txt} varsize={res}%f{txt})\n", S.acceleration, S.transform, S.tolerance, S.poolsize, cols(y))
		if (S.verbose>1) printf("{txt} - Variables: {res}" + invtokens(vars[i..j])+"{txt}\n")
		y = (*accelerate)(S, y, transform) :* stdevs
		if (S.num_iters_last_run>S.num_iters_max) S.num_iters_max = S.num_iters_last_run
		
		// Return residuals
		if (restore_dta) {
			residuals = residuals , y
		}
		else {
			if (S.verbose>1) printf("{txt} - Saving transformed variables\n")
			st_store(., st_addvar("double", newvars[i..j]), y)
		}
		y = J(0,0,.) // clear space (not really needed)

	}

	// this is max(iter) across all vars
	if (S.verbose==0) printf("{txt}(converged in %g iterations)\n", S.num_iters_max)

	// Restore dataset; used by hdfe.ado with generate() instead of clear
	if (restore_dta) {
		stata("restore")
		st_store(S.uid, st_addvar("double", newvars), residuals)
		residuals = J(0,0,.) // clear space (not really needed)
	}

	// Add chars with original names to new variables
	for (i=1;i<=Q;i++) {
		st_global(sprintf("%s[name]", newvars[i]), chars[i])
	}

}
end
