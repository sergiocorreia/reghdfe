mata:

// --------------------------------------------------------------------------
// LSQR: least squares solver of Ax=b
// --------------------------------------------------------------------------

// - Method derived from: https://web.stanford.edu/group/SOL/software/lsqr/lsqr-toms82a.pdf
// - Code uses scaffolding created by SAC in LSMR.mata
// - Most code adapted from Julia LSQR code 
// - Julia code: https://github.com/JuliaMath/IterativeSolvers.jl/blob/master/src/lsqr.jl
// - Note: If atol & btol are 1e-9, residual norm should be accurate to about 9 digits


`Solution' lsqr(`FixedEffects' A, 
				`Matrix' b,
				`Integer' miniter,
			  | `Integer' maxiter,
				`Real' atol,
				`Real' btol,
				`Real' ctol,
				`Integer' verbose)
{


	// Essential type definitions
	`Integer' 	m, n, iter
	`Vector' 	x
	`Solution'	sol

	// Type definitions for LSQR
	`Real'		α, β, ρ, ρ1, θ, ρ_bar, φ, φ_bar, φ_bar1, ζ, ζ_bar
	`Real'		rhs, δ, γ_bar, γ, cs2, sn2, extra_var
	`Vector' 	u, v, w
	`Matrix'	P, P_bar

	// Type definitions for residual norm and stopping criteria
	`Real' 		norm_r, norm_A, norm_b, norm_x, norm_r1, norm_r2, norm_dd, norm_xx, cond_A
	`Real' 		norm_Ar, res1, res2, r1sq, test1, test2 
	`Real' 		ρ_temp, w_ρ, cs1, sn1, ψ, τ, t1, t2, test3, rtol 

	`String'	msg
	
	assert(A.weights != .)

	// Initialize constants and other scalars
	m = A.num_rows()
	n = A.num_cols()
	sol = Solution()

	// Default parameters
	if (miniter == .)	miniter = 0
	if (maxiter == .) 	maxiter = 1e5
	if (atol == .) 		atol = 1e-6
	if (btol == .) 		btol = 1e-6
	if (ctol == .) 		ctol = 1e-6
	if (verbose == .) 	verbose = 0

	// Tolerances cannot exceed roundoff error (machine precision)
	if (atol < epsilon(1)) atol = epsilon(1)
	if (btol < epsilon(1)) btol = epsilon(1)
	if (ctol < epsilon(1)) ctol = epsilon(1)

	// Sanity checks
	assert_msg(m == rows(b), "Matrix A and vector b not conformable")
	assert_msg(inrange(miniter, 0, 1e10), "miniter outside valid range: [0, 1e+10]")
	assert_msg(inrange(maxiter, 1, 1e10), "maxiter outside valid range: [1, 1e+10]")
	assert_msg(miniter <= maxiter, "miniter should be below maxiter")
	assert_msg(inrange(atol, 1e-16, 1e-1), "atol outside valid range: [1e-16, 0.1]")
	assert_msg(inrange(btol, 1e-16, 1e-1), "btol outside valid range: [1e-16, 0.1]")
	assert_msg(inrange(ctol, 1e-16, 1e-1), "ctol outside valid range: [1e-16, 0.1]")
	assert_in(verbose, (-1, 0, 1, 2, 3, 4), "-verbose- must be an integer between -1 and 4")
	if (verbose>1) logging_lsqr("start")
	// if (verbose > 0) printf("\n{txt}## Solving linear system via LSQR:\n")

	// Initial LSQR iteration
	u = normalize(b, A.weights, β=., "Failed βu = b") // βu = b 
	v = normalize(A.mult_transpose(u), A.weights, α=., "Failed αv = A'u") // αv = A'u 

	// New section -- if norm(b) is zero, just exit 
	norm_b = vector_norm(b, A.weights, "Failed ||b||")
	if (norm_b == 0) {
		if (verbose>1) logging_lsqr("end")
		if (verbose>0) printf("{txt}   - LSQR converged in 0 iterations (trivial solution: 'b' is zero)\n")
		sol.stop_code = 1 // missing values
		sol.alphas = J(n, 1, 0) // could also just do "sol.alphas = x"
		sol.data = b
		sol.iteration_count = 0
		return(sol)
	}
	
	// Initialize vars for LSQR
	x 		= J(n, 1, 0)
	w 		= v 
	norm_A 	= cond_A = norm_dd = res2 = norm_x = norm_xx = sn2 = ζ = 0
	cs2 	= -1
	norm_Ar = α * β
	ρ_bar 	= α
	φ_bar	= norm_b = norm_r = norm_r1 = norm_r2 = β
	
	// Start the loop
	for (iter=1; iter<=maxiter; iter++) {

		// Lanczos process: generate vector v and scalars α, β 

		// Step 3: Bidiagonalization --> copied code from Sergio 
		u = normalize(A.mult(v) - α * u, A.weights, β=., "Failed βu = Av - αu") // βu = Av - αu
		v = normalize(A.mult_transpose(u) - β * v, A.weights, α=., "Failed αv = A'u - βv") // αv = A'u - βv 
		
		// Step 4: Construct & apply next orthogonal transformation
		// The following section is the calculation of x we need 
 
		φ_bar1  = φ_bar 
		ρ 		= hypot(ρ_bar, β)

		// givens rotations function returns (c, -s \ s, c)
		P = givens_rotation(ρ_bar, β, ρ=.)

		// apply rotations
		// multiple rotations here because I need some extra values (compared to LSMR)
		assign(P * (0 \ -α),     θ = ., ρ_bar = .)

		assign(P * (φ_bar1 \ 0), φ = ., φ_bar = .) // we need φ for updating t1, τ 
		assign(P * (φ \ 0)     , τ = ., extra_var = .) // need τ in the calculation of norm_Ar 


		// Calculate some vars following Givens Rotation 
		t1 		= φ / ρ
		t2 		=-θ / ρ

		// update w, x 
		x 		= t1 * w + x 
		w 		= t2 * w + v 

		// Now we move to constructing norms that are required for stopping criteria 
		// Use a plane rotation on the right to elimatine the 
		// super-diagonal element (θ) of the upper-bidiaganol matrix to estimate norm(x)

		// Want to use Givens Rotations instead of keeping cs2, sn2 
		γ_bar = -ρ	
		P_bar 	= givens_rotation(γ_bar, θ, γ = .)


		assign(P_bar * (0 \ -ρ), δ = ., γ_bar = .)


		rhs 	= φ - (δ * ζ)
		ζ_bar	= rhs / γ_bar

		norm_x 	= sqrt(norm_xx + (ζ_bar^2))
		γ 		= hypot(γ_bar, θ)

		ζ 		= rhs / γ
		norm_xx = norm_xx + ζ^2 

		// Estimate cond(A_bar), ||r_bar||, and ||A_bar'r_bar||
		w_ρ     = w * (1/ρ)
		norm_dd = norm_dd + norm(w_ρ)

		norm_A 	= sqrt(norm_A^2 + α^2 + β^2) // can't use hypot() as written bc of 3 args 
		
		cond_A 	= norm_A * sqrt(norm_dd)
		res1 	= φ_bar ^ 2	
		norm_Ar = α * abs(τ)
		norm_r 	= sqrt(res1)

		// Distinguish between norms 
		r1sq 	= norm_r^2
		norm_r1 = sqrt(r1sq)
		norm_r2 = norm_r

		// Use these norms to estimate other quantities
		// These quantities --> 0 near a solution 
		test1	= norm_r / norm_b
		test2 	= norm_Ar / (norm_A * norm_r)
		test3 	= 1 / cond_A
		t1 		= test1 / (1 + norm_A*norm_x/norm_b)
		rtol	= btol + atol*norm_A*norm_x/norm_b
		if (verbose>1) logging_lsqr("iter", iter, test1, rtol, test2, atol, test3, ctol, norm_A)

		if (iter < miniter) {
			continue
		}

		if (test1 == . | test2 == .) {
			msg = sprintf("{err}@ LSQR stopped; missing values in rel_res or rel_neq\n")
			sol.stop_code = 12 // missing values
			break
		}
		else if (test1 <= rtol) {
			msg = sprintf("{txt}   - LSQR converged in %g iterations (criterion: relative residual)\n", iter)
			sol.stop_code = 2 // consistent systems (with exact solutions)
			break
		}
		else if (test2 <= atol) {
			msg = sprintf("{txt}   - LSQR converged in %g iterations (criterion: normal equation)\n", iter)
			sol.stop_code = 3 // inconsistent systems (least squares)
			break
		}
		else if (test3 <= ctol) {
			msg = sprintf("{err}@ LSQR stopped; A is ill-conditioned\n")
			sol.stop_code = 11 // ill-conditioned
			break
		}
	}
	if (sol.stop_code == 0) {
		assert_msg(iter == maxiter+1)
		iter = maxiter
		msg = sprintf("{err}@ LSQR stopped; maximum number of iterations reached\n")
		sol.stop_code = 13 // max-iter reached
	}
	if (verbose>1) logging_lsqr("end")
	if (verbose>0 | (verbose>-1 & sol.stop_code>=10)) printf(msg)
	
	sol.data = b - A.mult(x) // BUGBUG: this will use A LOT of space for a sec: 1) "x" 2) A.mult(x) 3) b 4) the substraction; sadly we don't have in-place substract
	swap(sol.alphas, x) // BUGBUG: we need to adjust the alphas to undo any preconditioning ; see https://web.stanford.edu/group/SOL/software/lsmr/
	sol.iteration_count = iter
	return(sol)
}


`Void' logging_lsqr(`String' step, | `Integer' iter, `Real' test1, `Real' rtol, `Real' test2, `Real' atol, `Real' test3, `Real' ctol, `Real' norm_A){
	`Integer' 	col
	`String' 	table_row, color1, color2, color3
	// See "help smcl##ascii"
	if (step == "start") {
		printf("\n{txt}## Solving linear system via LSQR:\n")
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
		color1 = (test1 <= rtol) ? "{res}" : "{txt}"
		color2 = (test2 <= atol) ? "{res}" : "{txt}"
		color3 = (test3 <= ctol) ? "{res}" : "{txt}"
		table_row = sprintf("{txt} {c |} %3.0f {col %3.0f}{c |}", iter, col = col + 8)
		
		table_row = table_row + sprintf("%s%11.0e{txt}{col %3.0f}{c |}", color1, test1, col = col + 13)
		table_row = table_row + sprintf("%6.0e{col %3.0f}{c |}", rtol, col = col + 10)
				
		table_row = table_row + sprintf("%s%11.0e{txt}{col %3.0f}{c |}", color2, test2, col = col + 13)
		table_row = table_row + sprintf("%6.0e{col %3.0f}{c |}", atol, col = col + 10)
		
		table_row = table_row + sprintf("%s%11.0e{txt}{col %3.0f}{c |}", color3, test3, col = col + 13)
		table_row = table_row + sprintf("%6.0e{col %3.0f}{c |}", ctol, col = col + 10)
		
		table_row = table_row + sprintf("%9.1f {col %3.0f}{c |}", norm_A, col = col + 11)
		printf(table_row + "\n")
	}
	displayflush()
}

end

