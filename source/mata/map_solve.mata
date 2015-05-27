mata:
mata set matastrict on

void function map_solve(`Problem' S, `Varlist' vars,
		| `Varlist' newvars, `Varlist' partial, `Boolean' save_fe) {
	`Integer' i, Q, Q_partial, offset, g
	`Group' y
	`FunctionPointer' transform, accelerate
	real rowvector stdevs
	`Varlist'	target
	`Varlist'	chars

	if (S.verbose>0) printf("{txt}{bf:mata: map_solve()}\n")
	assert_msg(S.N!=., "map_solve() needs to be run after map_precompute()")
	assert_msg(S.N==st_nobs(), "dataset cannot change after map_precompute()")

	// Load data
	// BUGBUG: This will use 2x memory for a while; partition the copy+drop based on S.groupsize?
	if (S.verbose>0) printf("{txt} - Loading variables into Mata\n")
	vars = tokens(vars)
	y = st_data(., vars)
	Q = cols(y)

	// Store chars var[name] that contain the original varname (e.g. L.var)
	chars = J(1, Q, "")
	for (i=1;i<=Q;i++) {
		chars[i] = st_global(sprintf("%s[name]", vars[i]))
	}

	st_dropvar(vars) // We need the new ones on double precision

	if (args()<3 | newvars=="") newvars = vars
	assert_msg(length(vars)==length(newvars), "map_solve error: newvars must have the same size as vars")

	// Load additional partialled-out regressors
	Q_partial = 0
	if (args()==4 & partial!="") {
		vars = tokens(partial)
		Q_partial = cols(vars)
		y = y , st_data(., vars)
		st_dropvar(vars) // We need the new ones on double precision		
	}

	// Storing FEs and returning them requires 6 changes
	// 1) Extend the S and FE structures (add S.storing_betas, FE.alphas FE.tmp_alphas)
	// 2) Allocate them here
	// 3) Return results at the end of this function
	// 4) Within the accelerators, modify the SD to update the alphas
	// 5) Within map_projection, add a conditional to update tmp_alphas if needed
	S.storing_betas = 0
	if (args()<5) save_fe = 0
	assert_msg(save_fe==0 | save_fe==1, "map_solve error: save_fe must be either 0 or 1")
	if (save_fe) {
		assert_msg(partial=="", "map_solve error: partial must be empty if save_fe==1")
		assert_msg(length(vars)==1, "map_solve error: only one variable allowed if save_fe==1")
		if (S.verbose>0) printf("{txt} - Allocating objects to save the fixed effect estimates\n")
		S.storing_betas = 1
		for (g=1; g<=S.G; g++) {
			if (length(S.fes[g].target)>0) {
				S.fes[g].alphas = S.fes[g].tmp_alphas = 
					J(S.fes[g].levels, S.fes[g].has_intercept + S.fes[g].num_slopes, 0)
			}
		}
	}

	// Standardize all variables
	if (S.verbose>0) printf("{txt} - Standardizing variables\n")
	stdevs = J(1,cols(y),.)
	for (i=1; i<=cols(y); i++) {
		stdevs[i] = max((  sqrt(quadvariance(y[., i])) , sqrt(epsilon(1)) ))
	}
	y = y :/ stdevs

	// TODO: Report PEAK MEMORY that will be used
	// EG: de por si uso 2*size(y)
	// Luego bajo y asi que me quedo con y + ..
	// Contar cuantos vectores creo en los aceleradores y en los proyectores

	if (S.verbose>0) printf("{txt} - Solving problem (acceleration={res}%s{txt}, transform={res}%s{txt} tol={res}%-1.0e{txt} poolsize={res}%f{txt} varsize={res}%f{txt})\n", save_fe ? "steepest_descent" : S.acceleration, save_fe ? "kaczmarz" : S.transform, S.tolerance, S.groupsize, cols(y))

	// Warnings
	if (S.transform=="kaczmarz" & S.acceleration=="conjugate_gradient") {
		printf("{err}(warning: convergence is {bf:unlikely} with transform=kaczmarz and accel=CG)\n")
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

	// Call acceleration routine
	if (save_fe) {
		y = accelerate_sd(S, y, &transform_kaczmarz()) :* stdevs // Only these were modified to save FEs
		S.num_iters_max = S.num_iters_last_run
	}
	else if (S.groupsize>=cols(y)) {
		y = (*accelerate)(S, y, transform) :* stdevs
		S.num_iters_max = S.num_iters_last_run
	}
	else {
		S.num_iters_last_run = 0
		for (i=1;i<=cols(y);i=i+S.groupsize) {
			offset = min((i + S.groupsize - 1, cols(y)))
			if (S.verbose>1) printf("{txt} - Variables: {res}" + invtokens(vars[i..offset])+"{txt}\n")
			y[., i..offset] = (*accelerate)(S, y[., i..offset], transform) :* stdevs[i..offset]
			if (S.num_iters_last_run>S.num_iters_max) S.num_iters_max = S.num_iters_last_run
		}
	}

	// this is max(iter)0 for all vars
	if (S.verbose==0) printf("{txt}(converged in %g iterations)\n", S.num_iters_last_run)

	// Partial-out variables
	assert(Q_partial==0) // DISABLED FOR NOW DUE TO MEMORY ISSUES
	// if (Q_partial>0) {
	// 	if (S.verbose>1) printf("{txt} - Partialling out variables\n")
	// 	assert(cols(y)==Q+Q_partial)
	// 	y = y[., 1..Q] - y[., (Q+1)..cols(y)] * qrsolve(y[., (Q+1)..cols(y)] , y[., 1..Q])
	// 	stdevs =  stdevs[1..Q]
	// }

	// Store variables in dataset; do it by blocks to avoid 2x memory consumption
	if (S.verbose>1) printf("{txt} - Saving transformed variables\n")
	i = 1
	while (cols(y)>0) {
		if (S.groupsize>=cols(y)) {
			st_store(., st_addvar("double", newvars[i..length(newvars)]), y)
			y = J(0,0,.) // clear space
		}
		else {
			st_store(., st_addvar("double", newvars[i..(i+S.groupsize-1)]), y[., 1..(S.groupsize)])
			y = y[., (S.groupsize+1)..cols(y)] // clear space
		}
		i = i + S.groupsize
	}

	for (i=1;i<=Q;i++) {
		st_global(sprintf("%s[name]", newvars[i]), chars[i])
	}

	// Store FEs
	if (save_fe) {
		if (S.verbose>1) printf("{txt} - Saving fixed effects\n")
		for (g=1; g<=S.G; g++) {
			target = S.fes[g].target
			if (length(target)>0) {
				S.fes[g].tmp_alphas = J(0,0,.)
				S.fes[g].alphas = S.fes[g].alphas[ st_data(., S.fes[g].idvarname) , . ] :* stdevs
			}
		}
	}

}
end
