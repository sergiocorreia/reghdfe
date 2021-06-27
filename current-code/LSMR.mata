mata:

// --------------------------------------------------------------------------
// LSMR: least squares solver of Ax=b
// --------------------------------------------------------------------------

// Reference: 	http://web.stanford.edu/group/SOL/software/lsmr/
// - Julia:		https://github.com/JuliaMath/IterativeSolvers.jl/blob/master/src/lsmr.jl
// - Python:		https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.lsmr.html
// - Matlab:	 	https://github.com/timtylin/lsmr-SLIM/blob/master/lsmr.m
// - Fotran:		https://web.stanford.edu/group/SOL/software/lsmr/f90/lsmrf90.zip
// (note that most implementations are extremely similar to each other, up to comments)

// Usage:
// x = lsmr(A, B | , converged=?, maxiter=?, atol, btol, conlim, verbose)

// TODO:
// - Add preconditioning [PARTLY]
// - Allow A.weights [DONE]
// - Allow matrix B [DONE] -> REMOVED
// - Allow initial values for -x-
//	 This would be VERY (!!) useful when computing regressions with multiple LHSs ..
//	 like "reghdfe y x1, ... savecache(xyz)" and then reghdfe y x1 ... "usecache()".
//	 We just feed the last Xs as initial values. They don't even have to be the same varibles!
//	 So if e.g. we now include x2, we can always reuse what we had for y and x1.
// - Allow partial reorthogonalization
// - Standardize inputs???
// - Output: we want several objects:
//		- r=b-Ax
//		- in general we don't care about x, but maybe sometimes?? not sure
//		- standard errors of some xs?
//		- estimate of the ||A|| which the lsqr paper says might be useful for debugging
//		- estimate of the condition number of A (for which of the partial outs?)
//		- ...
//	  Problem: If we compute "r" maybe then compute ||r|| directly instead of relying on the product of sines formula
// - General: add -fast- option that replaces QR by the faster alternative (at the end, not part of lsmr)

// Concepts:

// residual: r = b - Ax
// relative residual: norm(r) / norm(b)  --> useful measure of how are we doing (reported by Matlab)


`Solution' lsmr(`FixedEffects' A,
			    `Matrix' b,
				`Integer' miniter,
			  | `Integer' maxiter,
			    `Real' atol,
			    `Real' btol,
			    `Real' conlim, // stop iterations if estimate of cond(A) exceeds conlim; intended to regularize ill-conditioned systems
			    `Integer' verbose)
{

	// Type definitions
	`Integer' 	m, n, iter
	`Vector' 	x
	`Solution'	sol

	// Type definitions (LSMR)
	`Real'		α, β, ρ, θ, ζ, α_bar, ρ_bar, θ_bar, ζ_bar
	`Real'		last_ρ, last_ρ_bar
	`Vector' 	u, v, h, h_bar
	`Matrix'	P, P_bar, P_tilde // Givens rotations (2x2 matrices defined by c and s)

	// Type definitions (residual norm and stopping criteria)
	`Real' 		β_hat, β_umlaut, ρ_dot, θ_tilde, β_dot, β_tilde, τ_tilde, τ_dot, ρ_tilde
	`Real' 		last_ζ, last_θ_tilde
	`Real' 		norm_A 		// Estimate of Frobenius norm ||A|| = sqrt(trace(A'A)); shouldn't exceed sqrt(n) if columns of A have been scaled to length 1
	`Real' 		norm_r 		// Estimate of ||r||
	`Real' 		norm_At_r 	// Estimate of ||A'r||
	`Real' 		norm_x 		// Estimate of ||x||
	`Real' 		cond_A 		// Estimate of cond(A)
	`Real' 		norm_b 		// Actual ||b||
	`Real' 		max_ρ_bar, min_ρ_bar, norm_A2
	`Real' 		ρ_temp
	`Real'		rel_res, rtol, rel_neq

	`String'	msg

	assert(A.weights != .)

	// Initialize constants and other scalars
	m = A.num_rows()
	n = A.num_cols()
	sol = Solution()

	// Default parameters
	if (miniter == .) miniter = 0
	if (maxiter == .) maxiter = 1e5  // Julia uses "m", Matlab uses "min(m, 20)" (recall that in exact arithmetic LSQR takes at most "m" iterations)
	if (atol == .) atol = 1e-6
	if (btol == .) btol = 1e-6
	if (conlim == .) conlim = 1e8
	if (verbose == .) verbose = 0

	// Tolerances cannot exceed roundoff error (machine precision)
	if (atol < epsilon(1)) atol = epsilon(1)
	if (btol < epsilon(1)) btol = epsilon(1)
	if (conlim == 0 | conlim > 1 / epsilon(1)) conlim = 1 / epsilon(1)

	// Sanity checks
	assert_msg(m == rows(b), "Matrix A and vector b not conformable")
	assert_msg(inrange(miniter, 0, 1e10), "miniter outside valid range: [0, 1e+10]")
	assert_msg(inrange(maxiter, 1, 1e10), "maxiter outside valid range: [1, 1e+10]")
	assert_msg(miniter <= maxiter, "miniter should be below maxiter")
	assert_msg(inrange(atol, epsilon(1), 1e-1), "atol outside valid range: [1e-16, 0.1]")
	assert_msg(inrange(btol, epsilon(1), 1e-1), "btol outside valid range: [1e-16, 0.1]")
	assert_msg(inrange(conlim, 1, 1e+16), "conlim outside valid range: [1, 1e+16]")
	assert_in(verbose, (-1, 0, 1, 2, 3, 4), "-verbose- must be an integer between -1 and 4")

	if (verbose>1) logging_lsmr("start")

	// Initial LSMR iteration
	u = normalize(b, A.weights, β=., "Failed βu = b")
	// BUGBUG: if β=0 then b=0 and we stop as any -x- solves this system
	v = normalize(A.mult_transpose(u), 1, α=., "Failed αv = A'u")
	// BUGBUG: if α=0 then normal equations already hold and we have found the solution
	
	// Initialize variables for LSMR
	h = v
	h_bar = x = J(n, 1, 0)
	α_bar = α
	ζ_bar = α * β
	ρ = ρ_bar = 1
	P_bar = I(2)

	// Initialize variables for estimation of residual norm ||r||
	β_umlaut = β
	β_dot = τ_tilde = θ_tilde = ζ = 0
	ρ_dot = 1
	ρ_tilde = 1 // This is not in the paper, but on iter=1 ρ_tilde=hypot(ρ_dot, 0)=1
	norm_b = vector_norm(b, A.weights, "Failed ||b||")

	if (norm_b == 0) {
		if (verbose>1) logging_lsmr("end")
		if (verbose>0) printf("{txt}   - LSMR converged in 0 iterations (trivial solution: 'b' is zero)\n")
		sol.stop_code = 1 // missing values
		sol.alphas = J(n, 1, 0) // could also just do "sol.alphas = x"
		sol.data = b
		sol.iteration_count = 0
		return(sol)
	}

	// Initialize variables for estimation of ||A|| and cond(A)
	norm_A = cond_A = norm_x = -1 // BUGBUG do we need to pre-initialize all these???? NOT SURE!
	norm_A2 = α ^ 2
	max_ρ_bar = 0
	min_ρ_bar= 1e100

	// Iterate Golub-Kahan
	for (iter=1; iter<=maxiter; iter++) {

		// Keep track of past values for Golub-Kahan
		last_ρ = ρ
		last_ρ_bar = ρ_bar

		// Keep track of past values for residual norm
		last_ζ = ζ
		last_θ_tilde = θ_tilde

		// 1) LSMR Algorithm:

		// Bidiagonalization
		assert_msg(!hasmissing(v), sprintf("iter !%g", iter)) // bugbug todo remove?
		u = normalize(A.mult(v) - α * u, A.weights, β=., "Failed βu = Av - αu")
		// BUGBUG We can stop if beta=0, as normal eqns hold!!!
		// if reorthogonalize: Store v in v_stack
		v = normalize(A.mult_transpose(u) - β * v, 1, α=., "Failed αv = A'u - βv")
		// if (reorthogonalize) reorthogonalize(v, V, iter)
		// BUGBUG We can stop if alpha=0, as normal eqns hold!!!

		// Construct rotation P from (α_bar, β)
		P = givens_rotation(α_bar, β, ρ=.)

		// Apply rotation P: (α, 0)  -->  (α_bar, θ)
		assign(P * (α \ 0), α_bar=., θ=.)

		// Apply rotation P_bar: (c_bar * ρ, θ)  -->  (ρ_temp, θ_bar)
		assign(P_bar * (ρ \ 0), ρ_temp=., θ_bar=.) // need to do this before overwriting P_bar below

		// Construct rotation P_bar from (ρ_temp, θ)
		P_bar = givens_rotation(ρ_temp, θ, ρ_bar=.)
		
		// Apply rotation P_bar: (0, ζ_bar) -->  (ζ_bar, ζ)
		assign(P_bar * (0 \ ζ_bar), ζ_bar=., ζ=.)

		// Update h, h_bar, x
		h_bar = h - (θ_bar * ρ) / (last_ρ * last_ρ_bar) * h_bar
		x = x + ζ / (ρ * ρ_bar) * h_bar
		h = v - (θ / ρ) * h

		// 2) Compute residual norm ||r||:

		// Apply rotation P: (0, β_umlaut)  -->  (β_umlaut, β_hat)
		assign(P * (0 \ β_umlaut), β_umlaut=., β_hat=.)

		if (iter > 1) {
			// Construct rotation P_tilde from (ρ_dot, θ_bar)
			P_tilde = givens_rotation(ρ_dot, θ_bar, ρ_tilde=.)
			
			// Apply rotation P_tilde: (ρ_bar, 0)  -->  (ρ_dot, θ_tilde)
			assign(P_tilde * (ρ_bar \ 0), ρ_dot=., θ_tilde=.)

			// Apply rotation P_tilde: (β_hat, β_dot)  -->  (β_dot, β_tilde)
			assign(P_tilde * (β_hat \ β_dot), β_dot=., β_tilde=.) // Note that "β_tilde" is never used
		}

		// Update t_tilde by forward substitution
		τ_tilde = (last_ζ - last_θ_tilde * τ_tilde) / ρ_tilde
		τ_dot = (ζ - θ_tilde * τ_tilde) / ρ_dot

		// 3) Estimate norms:

		// Estimate ||r||
		norm_r = hypot(β_dot - τ_dot, β_umlaut)
		
		// Estimate ||A||
		norm_A2 = norm_A2 + β ^ 2
		norm_A = sqrt(norm_A2)
		norm_A2 = norm_A2 + α ^ 2

		// Estimate cond(A)
		max_ρ_bar = max((max_ρ_bar, last_ρ_bar))
		if (iter > 1) min_ρ_bar = min((min_ρ_bar, last_ρ_bar))
		cond_A = max((max_ρ_bar, ρ_temp)) / min((min_ρ_bar, ρ_temp))
		
		// Estimate ||error|| = ||A'r||
		norm_At_r = abs(ζ_bar)

		// Estimate ||x||
		norm_x = vector_norm(x, 1, "Failed ||x||")

		// 4) Test for convergence
		rel_res = norm_r / norm_b
		rtol = btol + atol * norm_A * norm_x / norm_b
		rel_neq = norm_At_r / norm_A / norm_x
		if (verbose>1) logging_lsmr("iter", iter, rel_res, rtol, rel_neq, atol, 1/cond_A, 1/conlim, norm_A)

		if (iter < miniter) {
			continue
		}

		if (rel_res == . | rel_neq == .) {
			msg = sprintf("{err}@ LSMR stopped; missing values in rel_res or rel_neq\n")
			sol.stop_code = 12 // missing values
			break
		}
		else if (rel_res <= rtol) {
			msg = sprintf("{txt}   - LSMR converged in %g iterations (criterion: relative residual)\n", iter)
			sol.stop_code = 2 // consistent systems (with exact solutions)
			break
		}
		else if (rel_neq <= atol) {
			msg = sprintf("{txt}   - LSMR converged in %g iterations (criterion: normal equation)\n", iter)
			sol.stop_code = 3 // inconsistent systems (least squares)
			break
		}
		else if (cond_A >= conlim) {
			msg = sprintf("{err}@ LSMR stopped; A is ill-conditioned\n")
			sol.stop_code = 11 // ill-conditioned
			break
		}
	}

	if (sol.stop_code == 0) {
		assert_msg(iter == maxiter+1)
		iter = maxiter
		msg = sprintf("{err}@ LSMR stopped; maximum number of iterations reached\n")
		sol.stop_code = 13 // max-iter reached
	}

	if (verbose>1) logging_lsmr("end")
	if (verbose>0 | (verbose>-1 & sol.stop_code>=10)) printf(msg)
	
	sol.data = b - A.mult(x) // BUGBUG: this will use A LOT of space for a sec: 1) "x" 2) A.mult(x) 3) b 4) the substraction; sadly we don't have in-place substract
	swap(sol.alphas, x) // BUGBUG: we need to adjust the alphas to undo any preconditioning ; see https://web.stanford.edu/group/SOL/software/lsmr/
	sol.iteration_count = iter
	return(sol)
}


`Void' logging_lsmr(`String' step, | `Integer' iter, `Real' rel_res, `Real' rtol, `Real' rel_neq, `Real' atol, `Real' invcond, `Real' ctol, `Real' norm_A)
{
	`Integer' 	col
	`String' 	table_row, color1, color2, color3

	// See "help smcl##ascii"
	if (step == "start") {
		printf("\n{txt}## Solving linear system via LSMR:\n")
		printf(" {txt}{c TLC}{hline 5}{c TT}{hline 22}{c TT}{hline 22}{c TT}{hline 22}{c TT}{hline 10}{c TRC}\n")
		printf(" {txt}{c |}{space 5}{c |}        Rule 1        {c |}        Rule 2        {c |}        Rule 3        {c |}{space 10}{c |}\n")
		printf(" {txt}{c |}  i  {c LT}{hline 12}{c TT}{hline 9}{c +}{hline 12}{c TT}{hline 9}{c +}{hline 12}{c TT}{hline 9}{c RT}   ||A||  {c |}\n")
		printf(" {txt}{c |}{space 5}{c |}  rel.res.  {c |}   rtol  {c |}  rel.neq.  {c |}   atol  {c |}  1/cond(A) {c |}   ctol  {c |}{space 10}{c |}\n")
		printf(" {txt}{c LT}{hline 5}{c +}{hline 12}{c +}{hline 9}{c +}{hline 12}{c +}{hline 9}{c +}{hline 12}{c +}{hline 9}{c +}{hline 10}{c RT}\n")
	}
	else if (step == "end") {
		printf(" {txt}{c BLC}{hline 5}{c BT}{hline 12}{c BT}{hline 9}{c BT}{hline 12}{c BT}{hline 9}{c BT}{hline 12}{c BT}{hline 9}{c BT}{hline 10}{c BRC}\n")
	}
	else {
		col = 0
		color1 = (rel_res <= rtol) ? "{res}" : "{txt}"
		color2 = (rel_neq <= atol) ? "{res}" : "{txt}"
		color3 = (invcond <= ctol) ? "{res}" : "{txt}"
		table_row = sprintf("{txt} {c |} %3.0f {col %3.0f}{c |}", iter, col = col + 8)
		
		table_row = table_row + sprintf("%s%11.0e{txt}{col %3.0f}{c |}", color1, rel_res, col = col + 13)
		table_row = table_row + sprintf("%6.0e{col %3.0f}{c |}", rtol, col = col + 10)
		//table_row = table_row + sprintf("%8.1f{col %3.0f}{c |}", 100 * rtol / rel_res, col = col + 10)
		
		table_row = table_row + sprintf("%s%11.0e{txt}{col %3.0f}{c |}", color2, rel_neq, col = col + 13)
		table_row = table_row + sprintf("%6.0e{col %3.0f}{c |}", atol, col = col + 10)
		//table_row = table_row + sprintf("%8.1f{col %3.0f}{c |}", 100 * atol / rel_neq, col = col + 10)
		
		table_row = table_row + sprintf("%s%11.0e{txt}{col %3.0f}{c |}", color3, invcond, col = col + 13)
		table_row = table_row + sprintf("%6.0e{col %3.0f}{c |}", ctol, col = col + 10)
		//table_row = table_row + sprintf("%8.1f{col %3.0f}{c |}", 100 * ctol / invcond, col = col + 10)
		
		table_row = table_row + sprintf("%9.1f {col %3.0f}{c |}", norm_A, col = col + 11)
		printf(table_row + "\n")
	}

	displayflush()
}

end
