// --------------------------------------------------------------------------
// Driscoll-Kraay Standard Errors
// --------------------------------------------------------------------------
// Driscoll-Kraay (1998) standard errors are robust to both cross-sectional
// and serial correlation. They apply a Newey-West style HAC correction to
// cross-sectionally averaged moment conditions.
//
// This is asymptotically equivalent to clustering on the time variable,
// but with a HAC correction for autocorrelation in the time-aggregated moments.
//
// References:
// - Driscoll & Kraay (1998), "Consistent Covariance Matrix Estimation with
//   Spatially Dependent Panel Data", REStat
// - Hoechle (2007), "Robust Standard Errors for Panel Regressions with
//   Cross-Sectional Dependence", Stata Journal

mata:


`Void' reghdfe_vce_dkraay(`FixedEffects' S,
                          `Solution' sol,
                          `Matrix' D,
                          `Matrix' X,
                          `Variable' w,
                          `String' vce_mode,
                          `Integer' bw)
{
    `Matrix'                h_t, M, Omega_j
    `Integer'               T, t, j, K, lags
    `Real'                  weight_j, dof_adj
    `Vector'                resid
    `Factor'                F_time
    `Vector'                sorted_resid
    `Matrix'                sorted_X
    `Integer'               N_clust

    // 0. Validate that we have a time variable
    assert_msg(S.timevar != "", "Driscoll-Kraay standard errors require a tsset time variable")

    // 1. Create factor for time variable (same approach as cluster VCE)
    F_time = factor(S.timevar, S.sample, 0, "", 1, 1, ., 0)
    T = F_time.num_levels

    // Convert bandwidth to lags: bandwidth = lags + 1 (matching ivreg2 convention)
    // If bw is missing, use default based on Newey-West (1994)
    if (missing(bw)) {
        bw = reghdfe_dkraay_default_bw(T)
        if (S.verbose > 0) printf("{txt}   - Using default bandwidth: %g\n", bw)
    }
    lags = bw - 1
    
    // Validate lags (must be >= 0 and < T)
    if (lags < 0) {
        printf("{txt}Warning: bandwidth must be >= 1; setting to 1 (0 lags)\n")
        lags = 0
    }
    if (lags >= T) {
        printf("{txt}Warning: bandwidth %g too large for %g time periods; reducing to %g\n", bw, T, T)
        lags = T - 1
    }

    // 2. Dimensions
    K = sol.K + sol.report_constant

    // 3. Prepare weighted residuals
    if (S.weight_type != "") {
        resid = sol.resid :* w
    }
    else {
        resid = sol.resid
    }

    // 4. Sort data by time factor and set up panel
    F_time.panelsetup()
    sorted_resid = F_time.sort(resid)
    sorted_X = sol.report_constant ? F_time.sort((X, J(rows(X), 1, 1))) : F_time.sort(X)

    // 5. Compute h_t = sum of X_i * e_i for each time period t
    // h_t is a T x K matrix where each row is the cross-sectional sum for period t
    h_t = J(T, K, 0)
    for (t = 1; t <= T; t++) {
        h_t[t, .] = quadcross(1, 0, 
            panelsubmatrix(sorted_resid, t, F_time.info),
            panelsubmatrix(sorted_X, t, F_time.info), 0)
    }

    // 6. Compute the "Meat" with Bartlett/Newey-West kernel
    // S_0: Lag 0 component (cross-sectional robust part)
    M = quadcross(h_t, h_t)

    // S_j: Autocovariance components with Bartlett declining weights
    for (j = 1; j <= lags; j++) {
        weight_j = 1 - j / (lags + 1)  // Bartlett kernel weight
        
        // Omega_j = sum over t of h_t' * h_{t-j}
        // h_t[j+1..T, .] are the "current" periods, h_t[1..T-j, .] are the "lagged"
        Omega_j = quadcross(h_t[|j+1, 1 \ T, .|], h_t[|1, 1 \ T-j, .|])
        
        // Add symmetric contribution
        M = M + weight_j * (Omega_j + Omega_j')
    }

    // 7. Degrees of freedom adjustment
    // For Driscoll-Kraay, the relevant "cluster" count is T (number of time periods)
    // We use T/(T-1) as the cluster adjustment, similar to standard clustering
    N_clust = T
    
    // ivreg2-style adjustment:
    dof_adj = (sol.N - 1) / (sol.N - sol.df_m - S.df_a) * N_clust / (N_clust - 1)
    // Alternative: xtscc-style adjustment (does not subtract df_a for absorbed FEs)
    // dof_adj = (sol.N / (sol.N - sol.df_m)) * (N_clust / (N_clust - 1))
    
    if (vce_mode == "vce_asymptotic") dof_adj = N_clust / (N_clust - 1)

    if (S.verbose > 0) {
        printf("{txt}# Estimating Driscoll-Kraay Variance-Covariance Matrix\n\n")
        printf("{txt}   - Time variable: {res}%s{txt}\n", S.timevar)
        printf("{txt}   - Number of time periods: {res}%g{txt}\n", T)
        printf("{txt}   - Bandwidth: {res}%g{txt} (lags = %g)\n", lags + 1, lags)
        printf("{txt}   - Small-sample-adjustment: q = (%g - 1) / (%g - %g) * %g / (%g - 1) = %g\n", 
               sol.N, sol.N, sol.df_m + S.df_a, N_clust, N_clust, dof_adj)
    }

    // 8. Sandwich assembly: V = D * M * D
    sol.V = D * M * D * dof_adj
    
    // 9. Ensure the matrix is positive semi-definite
    if (min(Re(eigenvalues(sol.V))) < 0) {
        if (S.verbose > 0) printf("{txt}   - Applying PSD correction\n")
        (void) reghdfe_fix_psd(sol.V, sol.report_constant)
    }

    // 10. Update degrees of freedom (use T-1 as for clustering on time)
    if (sol.df_r > N_clust - 1) {
        sol.df_r = N_clust - 1
    }

    // 11. Store results for e() posting
    sol.N_clust = N_clust
    sol.N_clust_list = N_clust  // Required for display (as a scalar, will be treated as 1x1)
    sol.num_clusters = 1
    sol.clustervars = (S.timevar)
    sol.dkraay_bw = lags + 1  // bandwidth = lags + 1

    if (S.verbose > 0) printf("\n")
}


// Compute default bandwidth for Driscoll-Kraay
// Based on Newey-West (1994) plug-in formula: floor(4 * (T/100)^(2/9)) + 1
// Note: bandwidth = lags + 1, so we add 1 to the lag formula
// Alternative (used by fixest): floor(T^0.25) + 1
`Integer' reghdfe_dkraay_default_bw(`Integer' T)
{
    `Integer' lags
    lags = floor(4 * (T / 100)^(2/9))
    return(lags + 1)  // bandwidth = lags + 1
}

end
