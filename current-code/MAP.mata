// --------------------------------------------------------------------------
// Main code for MAP solver
// --------------------------------------------------------------------------

mata:

`Void' map_solver(`FixedEffects' S,
					  `Matrix' data,
					  `Integer' poolsize,
					  `FunctionP' fun_accel,
					  `FunctionP' fun_transform)
{
	`Integer'				i
	
		S.solution.converged = 0 // converged will get updated by check_convergence()

	// Speedup for constant-only case (no fixed effects)
	if (S.G==1 & S.factors[1].method=="none" & !S.factors[1].num_slopes & !(S.storing_alphas & S.factors[1].save_fe)) {
		assert(S.factors[1].num_levels == 1)
		// Not using poolsize here as this case is unlikely to happen with very large datasets
		data = data :- mean(data, S.factors[1].weights)
		S.solution.converged = 1
		S.solution.iteration_count = 1
	}
	else {
		if (poolsize >= cols(data)) {
			data = (*fun_accel)(S, data, fun_transform)
		}
		else {
			assert(S.solution.converged == 0) // check_convergence() sets converged=1 so we need to undo it
			for (i=1; i<=S.solution.K; i++) {
				S.solution.converged = 0
				data[., i] = (*fun_accel)(S, data[., i], fun_transform)
			}
		}
	}

	swap(S.solution.data, data)
}

end
