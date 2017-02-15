// Code that partials out (demean) a specific fixed effect
mata:

`Variables' panelmean(`Variables' y,
                      `Factor' f)
{
    pointer(`Variable')              Pw, Pcounts
    `Boolean' has_weights
    has_weights = asarray(f.extra, "has_weights")
    assert(has_weights==0 | has_weights==1) // bugbug: remove

    if (has_weights) {
        Pw = &asarray(f.extra, "weights")
        Pcounts = &asarray(f.extra, "weighted_counts")
        return(editmissing(`panelsum'(y, *Pw, f.info) :/ *Pcounts, 0))

    }
    else {
        return(`panelsum'(y, f.info) :/ f.counts)
    }
}


`Variables' panelsolve_invsym(`Variables' y,
                              `Factor' f,
                              `Boolean' has_intercept)
{
    `Integer'               i, L, K
    `Variables'             x, tmp_x, tmp_y, xbd, tmp_xbd
    `Variable'              w, tmp_w
    `Matrix'                xmeans
    `RowVector'             tmp_xmeans, tmp_ymeans
    `Matrix'                tmp_xx, tmp_xy
     `Boolean' has_weights
     has_weights = asarray(f.extra, "has_weights")
     assert(has_weights==0 | has_weights==1) // bugbug: remove
   
    // x, y and w must be already sorted by the factor f
    L = f.num_levels
    xbd = J(rows(y), cols(y), .)
    x = asarray(f.extra, "x")
    K = cols(x)
    //xx = asarray(f.extra, "xx")

    if (has_weights) w = asarray(f.extra, "weights")
    if (has_intercept) xmeans = asarray(f.extra, "xmeans")
    
    for (i = 1; i <= L; i++) {
            tmp_y = panelsubmatrix(y, i, f.info)
            tmp_x = panelsubmatrix(x, i, f.info)
            tmp_w = has_weights ? panelsubmatrix(w, i, f.info) : 1
            if (has_intercept) {
                tmp_ymeans = mean(tmp_y, tmp_w)
                tmp_xmeans = K > 1 ? xmeans[i, .] : xmeans[i]
                tmp_xx = quadcrossdev(tmp_x, tmp_xmeans, tmp_w, tmp_x, tmp_xmeans)
                // Bugbug: PRECOMPUTE! invsym(tmp_xx) !!!!!!!!!!!
                // crossdev or quadcrossdev?
                tmp_xy = quadcrossdev(tmp_x, tmp_xmeans, tmp_w, tmp_y, tmp_ymeans)
                tmp_xbd = (tmp_x :- tmp_xmeans) * (invsym(tmp_xx) * tmp_xy) :+ tmp_ymeans
            }
            else {
                tmp_xx = quadcross(tmp_x, tmp_w, tmp_x)
                // Bugbug: PRECOMPUTE! invsym(tmp_xx) !!!!!!!!!!!
                tmp_xy = quadcross(tmp_x, tmp_w, tmp_y)
                tmp_xbd = tmp_x * (invsym(tmp_xx) * tmp_xy)
            }
            xbd[|f.info[i,1], 1 \ f.info[i,2], .|] = tmp_xbd
    }
    return(f.invsort(xbd))
}

/*
`Variables' panelsolve_qrsolve(`Variables' Y, `Variables' X, `Factor' f)
{
    `Integer'               i
    `Variables'             x, y, betas

    betas = J(f.num_levels, 1 + cols(X), .)
    
    for (i = 1; i <= f.num_levels; i++) {
            y = panelsubmatrix(Y, i, F.info)
            x = panelsubmatrix(X, i, F.info) , J(rows(y), 1, 1)
            betas[i, .] = qrsolve(x, y)'
    }
    return(betas)
}

*/
end

/*

`Group' function map_projection(`Problem' S, `Integer' g, `Group' y) {
    `Integer'   K, L, Q // Q is the number of depvars
    `Integer'   j, i_lower, i_upper // j loops over levels, i loops over observations
    `Boolean'   has_weights, sortedby, has_intercept, storing_betas
    `Series'    sorted_w
    `Group'     ans
    `Vector'    tmp_w, tmp_count
    real rowvector b
    real rowvector ymean, alpha // 1*Q
    real rowvector zero // 1*K
    `Matrix'    tmp_y, tmp_x
    pointer(`Series') scalar p_sorted_w
    pragma unset sorted_w // If we just set the pointer, what happens to the underlying data?

    // PROFILE TO SEE IF THIS HELPS OR NOT AT ALL
    //pointer(`Vector') scalar p_offset
    //p_offset = &(S.fes[g].offsets)

    has_weights = S.weightvar !=""
    sortedby = S.fes[g].is_sortedby
    has_intercept = S.fes[g].has_intercept
    K = S.fes[g].num_slopes
    Q = cols(y)
    L = S.fes[g].levels
    tmp_w = 1 // Harmless value for when there are no weights

    // Minimize copy+order operations on y
    if (has_weights) p_sorted_w = sortedby ? &(S.w) : &(sorted_w = S.w[S.fes[g].p, .])
    if (K>0) zero = J(1,K,0)

    ans = sortedby ? y : y[S.fes[g].p, .]

    i_lower = 1
    storing_betas = S.storing_betas & length(S.fes[g].target)>0
    for (j=1; j<=L; j++) {
        i_upper = S.fes[g].offsets[j]
        tmp_count = S.fes[g].counts[j]
        
        if (has_weights) tmp_w = (*p_sorted_w)[| i_lower \ i_upper |]
        tmp_y = ans[| i_lower , 1 \ i_upper , . |]
        // BUGBUG: quadcolsum or colsum ? Depends if there are dense FEs. Maybe condition it on L??
        if (has_weights) {
            ymean = has_intercept ? (quadcolsum(tmp_y :* tmp_w) / tmp_count) : 0
        }
        else {
            ymean = has_intercept ? (quadcolsum(tmp_y) / tmp_count) : 0
        }

        if (K>0) {
            tmp_x = S.fes[g].x[| i_lower , 1 \ i_upper , . |]
            // BUGBUG crossdev/cross or their quad version?
            if (has_intercept) {
                b = S.fes[g].inv_xx[| 1+(j-1)*K , 1 \ j*K , . |] * quadcrossdev(tmp_x, zero, tmp_w, tmp_y, ymean)
                alpha = ymean - S.fes[g].xmeans[j, .] * b
            }
            else {
                b = S.fes[g].inv_xx[| 1+(j-1)*K , 1 \ j*K , . |] * quadcross(tmp_x, tmp_w, tmp_y)
            }
        }
        
        if (storing_betas) {
            if (has_intercept) S.fes[g].tmp_alphas[j, 1] = K==0 ? ymean : alpha
            if (K>0) S.fes[g].tmp_alphas[j, (has_intercept+1)..(has_intercept+K) ] = b'
        }

        // BUGBUG if we split this ternary will it be faster?
        //ans[| i_lower , 1 \ i_upper , . |] = K>0 ? (ymean :+ tmp_x*b) : (ymean :+ J(i_upper-i_lower+1,Q,0))
        if (K==0) {
            ans[| i_lower , 1 \ i_upper , . |] = ymean :+ J(i_upper-i_lower+1,Q,0)
        }
        else if (has_intercept) {
            ans[| i_lower , 1 \ i_upper , . |] = ymean :+ tmp_x*b
        }
        else {
            ans[| i_lower , 1 \ i_upper , . |] = tmp_x*b
        }


        i_lower = i_upper + 1
    }
        
    return(sortedby ? ans : ans[S.fes[g].inv_p, .])
}
end
*/
