mata:

// --------------------------------------------------------------------------
// Acceleration Schemes
// --------------------------------------------------------------------------

`Variables' function accelerate_test(`FixedEffects' S, `Variables' y, `FunctionP' T) {
	`Integer'	iter, g
	`Variables'		resid
	`Factor' f
	pragma unset resid

	assert(S.converged == 0)

	for (iter=1; iter<=S.maxiter; iter++) {
		for (g=1; g<=S.G; g++) {
			f = S.factors[g]
			if (g==1) resid = y - panelmean(f.sort(y), f)[f.levels, .]
			else resid = resid - panelmean(f.sort(resid), f)[f.levels, .]
		}
		if (check_convergence(S, iter, resid, y)) break
		y = resid
	}
	return(resid)
}

// --------------------------------------------------------------------------

`Variables' function accelerate_none(`FixedEffects' S, `Variables' y, `FunctionP' T) {
	`Integer'	iter
	`Variables'		resid
	pragma unset resid

	assert(S.converged == 0)

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
	assert(S.converged == 0)

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

	assert(S.converged == 0)
	if (S.timeit) timer_on(70)
	Q = cols(y)
	
	d = 1 // BUGBUG Set it to 2/3 // Number of recent SSR values to use for convergence criteria (lower=faster & riskier)
	// A discussion on the stopping criteria used is described in
	// http://scicomp.stackexchange.com/questions/582/stopping-criteria-for-iterative-linear-solvers-applied-to-nearly-singular-system/585#585

	if (S.timeit) timer_on(73)
	improvement_potential = weighted_quadcolsum(S, y, y)
	recent_ssr = J(d, Q, .)
	if (S.timeit) timer_off(73)
	
	if (S.timeit) timer_on(71)
	(*T)(S, y, r, 1)
	if (S.timeit) timer_off(71)
	if (S.timeit) timer_on(73)
	ssr = weighted_quadcolsum(S, r, r) // cross(r,r) when cols(y)==1 // BUGBUG maybe diag(quadcross()) is faster?
	u = r
	if (S.timeit) timer_off(73)

	for (iter=1; iter<=S.maxiter; iter++) {
		if (S.timeit) timer_on(71)
		(*T)(S, u, v, 1) // This is the hottest loop in the entire program
		if (S.timeit) timer_off(71)
		if (S.timeit) timer_on(73)
		alpha = safe_divide( ssr , weighted_quadcolsum(S, u, v) )
		if (S.timeit) timer_off(73)
		if (S.timeit) timer_on(74)
		recent_ssr[1 + mod(iter-1, d), .] = alpha :* ssr
		improvement_potential = improvement_potential - alpha :* ssr
		y = y - alpha :* u
		if (S.timeit) timer_off(74)
		if (S.timeit) timer_on(75)
		if (S.compute_rre & !S.prune) reghdfe_rre_benchmark(y[., 1], S.rre_true_residual, S.rre_depvar_norm)
		r = r - alpha :* v
		ssr_old = ssr
		if (S.timeit) timer_off(75)
		if (S.timeit) timer_on(73)
		if (S.verbose>=5) r
		ssr = weighted_quadcolsum(S, r, r)
		beta = safe_divide( ssr , ssr_old) // Fletcher-Reeves formula, but it shouldn't matter in our problem
		if (S.timeit) timer_off(73)
		u = r + beta :* u
		// Convergence if sum(recent_ssr) > tol^2 * improvement_potential
		if (S.timeit) timer_on(76)
		if ( check_convergence(S, iter, colsum(recent_ssr), improvement_potential, "hestenes") ) {
			break
			if (S.timeit) timer_off(76)
		}
		if (S.timeit) timer_off(76)
	}
	if (S.timeit) timer_off(70)
	return(y)
}

// --------------------------------------------------------------------------

`Variables' function accelerate_sd(`FixedEffects' S, `Variables' y, `FunctionP' T) {
	`Integer'	iter, g
	`Variables' proj
	`RowVector' t
	pragma unset proj

	assert(S.converged == 0)

	for (iter=1; iter<=S.maxiter; iter++) {
		(*T)(S, y, proj, 1)
		if (check_convergence(S, iter, y-proj, y)) break
		t = safe_divide( weighted_quadcolsum(S, y, proj) , weighted_quadcolsum(S, proj, proj) )
		if (uniform(1,1)<0.1) t = 1 // BUGBUG: Does this REALLY help to randomly unstuck an iteration?

		y = y - t :* proj
		if (S.compute_rre & !S.prune) reghdfe_rre_benchmark(y[., 1], S.rre_true_residual, S.rre_depvar_norm)

		if (S.storing_alphas) {
			for (g=1; g<=S.G; g++) {
				//g, ., ., t
				//asarray(S.factors[g].extra, "alphas"), asarray(S.factors[g].extra, "tmp_alphas")
				if (S.save_fe[g]) {
					asarray(S.factors[g].extra, "alphas",
					    asarray(S.factors[g].extra, "alphas") +
					    t :* asarray(S.factors[g].extra, "tmp_alphas")
					)
				}
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

	assert(S.converged == 0)
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
		// if (check_convergence(S, iter, accelerate? resid :* .  : resid, y)) break
		if (S.compute_rre & !S.prune) reghdfe_rre_benchmark(resid[., 1], S.rre_true_residual, S.rre_depvar_norm)
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
		update_error = max(mean(reldif(y_new, y_old), S.weight))
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

	S.converged = S.converged + (update_error <= S.tolerance)
	is_last_iter = iter==S.maxiter
	
	if (S.converged >= S.min_ok) {
		S.iteration_count = max((iter, S.iteration_count))
		S.accuracy = max((S.accuracy, update_error))
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
	return(S.converged >= S.min_ok)
}

// --------------------------------------------------------------------------

`Matrix' weighted_quadcolsum(`FixedEffects' S, `Matrix' x, `Matrix' y) {
	// BUGBUG: override S.has_weights with pruning
	// One approach is faster for thin matrices
	// We are using cross instead of quadcross but it should not matter for this use
	if (S.has_weights) {
		if (cols(x) < 14) {
			return(quadcross(x :* y, S.weight)')
		}
		else {
			return(diagonal(quadcross(x, S.weight, y))')
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
