// Code that partials out (demean) a specific fixed effect
mata:

`Variables' panelmean(`Variables' y,
                      `Factor' f)
{
    pointer(`Variable')              Pw, Pcounts
    `Boolean' has_weights
    has_weights = asarray(f.extra, "has_weights") == J(0,0,.) ? 0 : asarray(f.extra, "has_weights")
    assert(has_weights==0 | has_weights==1)

    if (has_weights) {
        Pw = &asarray(f.extra, "weights")
        Pcounts = &asarray(f.extra, "weighted_counts")
        return(editmissing(`panelsum'(y, *Pw, f.info) :/ *Pcounts, 0))
    }
    else {
        return(`panelsum'(y, f.info) :/ f.counts)
    }
}


`Matrix' precompute_inv_xx(`Factor' f,
                           `Boolean' has_intercept)
{
    `Integer'               i, L, K, offset
    `Variables'             x, tmp_x
    `Variable'              w, tmp_w
    `Matrix'                xmeans, inv_xx
    `RowVector'             tmp_xmeans
    `Matrix'                tmp_inv_xx
    `Boolean'               has_weights

    has_weights = asarray(f.extra, "has_weights")

    // x and w must be already sorted by the factor f
    x = asarray(f.extra, "x")
    L = f.num_levels
    K = cols(x)
    inv_xx = J(L * K, K, .)

    if (has_weights) w = asarray(f.extra, "weights")
    if (has_intercept) xmeans = asarray(f.extra, "xmeans")
    
    for (i = 1; i <= L; i++) {
            tmp_x = panelsubmatrix(x, i, f.info)
            tmp_w = has_weights ? panelsubmatrix(w, i, f.info) : 1
            if (has_intercept) {
                tmp_xmeans = K > 1 ? xmeans[i, .] : xmeans[i]
                tmp_inv_xx = invsym(quadcrossdev(tmp_x, tmp_xmeans, tmp_w, tmp_x, tmp_xmeans))
            }
            else {
                tmp_inv_xx = invsym(quadcross(tmp_x, tmp_w, tmp_x))
            }
            offset = K * (i - 1)
            inv_xx[|offset + 1, 1 \ offset + K , . |] = tmp_inv_xx
    }
    return(inv_xx)
}


`Variables' panelsolve_invsym(`Variables' y,
                              `Factor' f,
                              `Boolean' has_intercept,
                            | `Matrix' alphas)
{
    `Integer'               i, L, K, offset
    `Variables'             x, tmp_x, tmp_y, xbd, tmp_xbd
    `Variable'              w, tmp_w
    `Matrix'                xmeans, inv_xx
    `RowVector'             tmp_xmeans, tmp_ymeans
    `Matrix'                tmp_xy, tmp_inv_xx
    `Boolean'               has_weights
    `Boolean'               save_alphas
    `Vector'                b

    has_weights = asarray(f.extra, "has_weights")
    save_alphas = args()>=4 & alphas!=J(0,0,.)
    // assert(has_weights==0 | has_weights==1)
    if (save_alphas) assert(cols(y)==1)
   
    // x, y and w must be already sorted by the factor f
    L = f.num_levels
    xbd = J(rows(y), cols(y), .)
    x = asarray(f.extra, "x")
    inv_xx = asarray(f.extra, "inv_xx")
    K = cols(x)

    if (has_weights) w = asarray(f.extra, "weights")
    if (has_intercept) xmeans = asarray(f.extra, "xmeans")
    
    for (i = 1; i <= L; i++) {
            tmp_y = panelsubmatrix(y, i, f.info)
            tmp_x = panelsubmatrix(x, i, f.info)
            tmp_w = has_weights ? panelsubmatrix(w, i, f.info) : 1
            offset = K * (i - 1)
            tmp_inv_xx = inv_xx[|offset + 1, 1 \ offset + K , . |]

            if (has_intercept) {
                tmp_ymeans = mean(tmp_y, tmp_w)
                tmp_xmeans = K > 1 ? xmeans[i, .] : xmeans[i]
                tmp_xy = quadcrossdev(tmp_x, tmp_xmeans, tmp_w, tmp_y, tmp_ymeans)
                if (save_alphas) {
                    b = tmp_inv_xx * tmp_xy
                    alphas[i, .] = tmp_ymeans - tmp_xmeans * b, b'
                    tmp_xbd = (tmp_x :- tmp_xmeans) * b :+ tmp_ymeans
                }
                else {
                    tmp_xbd = (tmp_x :- tmp_xmeans) * (tmp_inv_xx * tmp_xy) :+ tmp_ymeans
                }
            }
            else {
                tmp_xy = quadcross(tmp_x, tmp_w, tmp_y)
                if (save_alphas) {
                    b = tmp_inv_xx * tmp_xy
                    alphas[i, .] = b'
                    tmp_xbd = tmp_x * b
                }
                else {
                    tmp_xbd = tmp_x * (tmp_inv_xx * tmp_xy)
                }
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

// used with lsmr if we have fixed slopes
`Variables' reghdfe_panel_precondition(`Variables' y, `Factor' f)
{
    `Vector' ans
    pointer(`Variable')              Pw
    `Boolean' has_weights

    has_weights = asarray(f.extra, "has_weights")
    if (has_weights) {
        Pw = &asarray(f.extra, "weights")
        ans = `panelsum'(y:^2, *Pw, f.info)
    }
    else {
        ans = `panelsum'(y, f.info)
    }

    ans = y :/ sqrt(ans)[f.levels]
    return(ans)
}
end

