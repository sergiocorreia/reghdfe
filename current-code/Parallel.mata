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
