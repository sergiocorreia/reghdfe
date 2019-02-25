mata:

// --------------------------------------------------------------------------
// LSMR estimation: Solve Ax=b with LS (ignore consistent case) (A?=y) (Z=D?)
// --------------------------------------------------------------------------
// Source: http://web.stanford.edu/group/SOL/software/lsmr/
// Code based on https://github.com/timtylin/lsmr-SLIM/blob/master/lsmr.m
// Copyright (BSD2): https://github.com/timtylin/lsmr-SLIM/blob/master/license.txt

// Requirements
// A(x, 1) = Ax     Projections "xβ"
// A(x, 2) = A'x    Sum of y by group; panelmean() if dummies and w/precond

`Vector' lsmr(`FixedEffects' S, `Vector' b, `Vector' x) {
    `Real' eps
    `Integer' iter // m, n
    `Real' beta, zetabar, alphabar, rho, rhobar, cbar, sbar
    `Real' betadd, betad, rhodold, tautildeold, thetatilde, zeta, d
    `Real' normA2, maxrbar, minrbar
    `Real' normb, normr
    `Real' test1, test2, test3
    `Vector' u, v, h, hbar

    `Real' alpha, alphahat, lambda, chat, shat, rhoold, c, s, thetanew, rhobarold, zetaold, stildeold
    `Real' thetabar, rhotemp, betaacute, betacheck, betahat, thetatildeold, rhotildeold, ctildeold, taud
    `Real' normA, normAr, condA, normx, rtol

    assert(cols(b)==1)
    if (S.verbose > 0) printf("\n{txt}## Computing LSMR\n\n")

    // Constants
    eps = epsilon(1)
    
    lambda = 0 // not used
    S.converged = 0

    beta = S.lsmr_norm(b)
    assert_msg(beta < . , "beta is missing")
    u = (beta > eps) ? (b / beta) : b
    v = S.lsmr_At_mult(u) // v = (*A)(u, 2)
    assert_msg(!missing(v), "-v- missing")
    // m = rows(v) // A is m*n
    // n = rows(u)

    alpha = S.lsmr_norm(v)
    assert_msg(alpha < . , "alpha is missing")
    if (alpha > eps) v = v / alpha

    // Initialize variables for 1st iteration.
    zetabar = alpha * beta
    alphabar = alpha
    rho = rhobar = cbar = 1
    sbar = 0

    h = v
    hbar = J(rows(h), 1, 0) // remove this
    //x = J(rows(h), 1, 0)

    // Initialize variables for estimation of ||r||
    betadd = beta
    betad = 0
    rhodold = 1
    tautildeold = 0
    thetatilde = 0
    zeta = 0
    d = 0

    // Initialize variables for estimation of ||A|| and cond(A)
    normA2  = alpha ^ 2
    maxrbar = 0
    minrbar = 1e+100

    // Items for use in stopping rules.
    normb = beta
    normr = beta

    // Exit if b=0 or A'b = 0.
    normAr = alpha * beta
    if (normAr == 0) {
        "DONE -> UPDATE THIS STOPPING CONDITION"
        return
    }

    if (S.verbose > 1) {
        "< < < <"
        test1 = 1
        test2 = alpha / beta
        printf(" %10.3e %10.3e\n", normr, normAr )
        printf("  %8.1e %8.1e\n" , test1, test2  )
        "> > > > "
    }

    // Main loop

    for (iter=1; iter<=S.maxiter; iter++) {

        // Update (1) βu = Av - αu (2) αv = A'u - βv
        u = S.lsmr_A_mult(v) - alpha * u // u = (*A)(v, 1) - alpha * u

        //"hash of u"
        //hash1(round(u*1e5))
        //u[1..5]

        beta = S.lsmr_norm(u)
        if (beta > eps) u = u / beta

        v = S.lsmr_At_mult(u) - beta * v // v = (*A)(u, 2) - beta * v
        alpha  = S.lsmr_norm(v)
        if (alpha > eps) v = v / alpha
        
        // α and β are now on iteration {k+1}

        // Construct rotation Qhat_{k, 2k+1}
        alphahat = S.lsmr_norm((alphabar, lambda))
        assert_msg(alphahat < . , "alphahat is missing")
        chat = alphabar / alphahat
        shat = lambda / alphahat

        // Use a plane rotation (Q_i) to turn B_i to R_i.
        rhoold = rho
        rho = norm((alphahat, beta))
        c = alphahat / rho
        s = beta / rho
        thetanew = s * alpha
        alphabar = c * alpha

        // Use a plane rotation (Qbar_i) to turn R_i^T to R_i^bar.
        rhobarold = rhobar
        zetaold = zeta
        thetabar = sbar * rho
        rhotemp = cbar * rho
        rhobar = norm((cbar * rho, thetanew))
        cbar = cbar * rho / rhobar
        sbar = thetanew / rhobar
        zeta =   cbar * zetabar
        zetabar = -sbar * zetabar

        // Update h, h_hat, x
        hbar = iter > 1 ? h - (thetabar * rho / (rhoold * rhobarold)) * hbar : h
        assert_msg(!missing(hbar), "hbar missing")
        x = iter > 1 ? x + (zeta / (rho * rhobar)) * hbar  : (zeta / (rho * rhobar)) * hbar
        assert_msg(!missing(x), "x missing")
        h = v - (thetanew / rho) * h

        // Estimate of ||r||
        
        // Apply rotation Qhat_{k,2k+1}
        betaacute =  chat * betadd
        betacheck = -shat * betadd

        // Apply rotation Q_{k,k+1}
        betahat   =  c * betaacute;
        betadd    = -s * betaacute;
          
        // Apply rotation Qtilde_{k-1}
        // betad = betad_{k-1} here
        thetatildeold = thetatilde
        rhotildeold = norm((rhodold, thetabar))
        ctildeold = rhodold / rhotildeold
        stildeold = thetabar / rhotildeold
        thetatilde = stildeold * rhobar
        rhodold = ctildeold * rhobar
        betad = -stildeold * betad + ctildeold * betahat
    
        // betad   = betad_k here
        // rhodold = rhod_k  here
        tautildeold   = (zetaold - thetatildeold * tautildeold) / rhotildeold
        taud          = (zeta - thetatilde * tautildeold) / rhodold
        d             = d + betacheck^2
        normr         = sqrt(d + (betad - taud)^2 + betadd^2)
        
        // Estimate ||A||.
        normA2        = normA2 + beta^2
        normA         = sqrt(normA2)
        normA2        = normA2 + alpha^2
    
        // Estimate cond(A)
        maxrbar = max((maxrbar,rhobarold))
        if (iter > 1) minrbar = min((minrbar,rhobarold))
        condA = max((maxrbar,rhotemp)) / min((minrbar,rhotemp))

        // Test for convergence.

        // Compute norms for convergence testing.
        normAr  = abs(zetabar)
        normx   = S.lsmr_norm(x)

        // Now use these norms to estimate certain other quantities,
        // some of which will be small near a solution.
        test1   = normr  / normb
        test2   = normAr / (normA*normr)
        test3   =      1 / condA
        rtol    = S.btol + S.tolerance *normA*normx / normb
    
        // The following tests guard against extremely small values of
        // atol, btol or ctol.  (The user may have set any or all of
        // the parameters atol, btol, conlim  to 0.)
        // The effect is equivalent to the normAl tests using
        // atol = eps,  btol = eps,  conlim = 1/eps.

        // Allow for tolerances set by the user.

        if  (test3 <= 1 / S.conlim) S.converged = 3
        if  (test2 <= S.tolerance) S.converged = 2
        if  (test1 <= rtol) S.converged = 1

        if (S.verbose > 1) {
            printf(" - Convergence: %g\n", S.converged)
            "iter normr normAr"
            iter, normr, normAr
            "test1 test2 test3"
            test1, test2, test3
            "criteria1 criteria2 criteria3"
            1/S.conlim , S.tolerance, rtol
            ">>>"
        }
        
        if (S.compute_rre & !S.prune) {
            reghdfe_rre_benchmark(b - S.lsmr_A_mult(x), S.rre_true_residual, S.rre_depvar_norm)
        }
        
        if (S.converged) break
    }

    if (!S.converged) {
        printf("\n{err}convergence not achieved in %g iterations (last error=%e); try increasing maxiter() or decreasing tol().\n", S.maxiter, test2)
        exit(430)
    }

    S.iteration_count = max((iter, S.iteration_count))

    u = b - S.lsmr_A_mult(x)
    return(u)
}
end
