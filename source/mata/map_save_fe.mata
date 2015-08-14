mata:
mata set matastrict on

void function map_save_fe(`Problem' S, `Varname' sum_fe) {

	// Storing FEs and returning them requires 6 changes
	// 1) Extend the S and FE structures (add S.storing_betas, FE.alphas FE.tmp_alphas)
	// 2) Allocate them here
	// 3) Return results at the end of this function
	// 4) Within the accelerators, modify the SD to update the alphas
	// 5) Within map_projection, add a conditional to update tmp_alphas if needed

	`Integer' 	g
	`Vector' 	y
	`Varlist'	target

	if (S.verbose>0) printf("{txt}{bf:mata: map_save_fe()}\n")
	assert_msg(S.N!=., "map_save_fe() needs to be run after map_precompute()")
	assert_msg(S.N==st_nobs(), "dataset cannot change after map_precompute()")

	if (S.verbose>0) printf("{txt} - Loading variable into Mata\n")
	assert_msg(length(tokens(sum_fe))==1, "sum_fe should be a varname, not a varlist")
	y = st_data(., sum_fe)
	
	st_dropvar(sum_fe)
	
	if (S.verbose>0) printf("{txt} - Allocating objects to save the fixed effect estimates\n")
	S.storing_betas = 1
	for (g=1; g<=S.G; g++) {
		if (length(S.fes[g].target)>0) {
			S.fes[g].alphas = S.fes[g].tmp_alphas = 
				J(S.fes[g].levels, S.fes[g].has_intercept + S.fes[g].num_slopes, 0)
		}
	}

	// No need to standardize b/c it's only one variable
	if (S.verbose>0) printf("{txt} - Solving problem (acceleration={res}%s{txt}, transform={res}%s{txt} tol={res}%-1.0e{txt} poolsize={res}%f{txt} varsize={res}%f{txt})\n", "steepest_descent", "kaczmarz", S.tolerance, S.poolsize, 1)

	y = accelerate_sd(S, y, &transform_kaczmarz())
	if (S.verbose==0) printf("{txt}(converged in %g iterations)\n", S.num_iters_last_run)
	
	if (S.verbose>1) printf("{txt} - Saving transformed variable\n")
	st_store(., st_addvar("double", sum_fe), y)
	y = J(0,0,.) // clear space

	// Store FEs
	if (S.verbose>1) printf("{txt} - Saving fixed effects\n")
	for (g=1; g<=S.G; g++) {
		target = S.fes[g].target
		if (length(target)>0) {
			S.fes[g].tmp_alphas = J(0,0,.)
			S.fes[g].alphas = S.fes[g].alphas[ st_data(., S.fes[g].idvarname) , . ]
		}
	}
}
end
