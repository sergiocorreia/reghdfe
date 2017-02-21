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
                              `Boolean' has_intercept,
                            | `Matrix' alphas)
{
    `Integer'               i, L, K
    `Variables'             x, tmp_x, tmp_y, xbd, tmp_xbd
    `Variable'              w, tmp_w
    `Matrix'                xmeans
    `RowVector'             tmp_xmeans, tmp_ymeans
    `Matrix'                tmp_xx, tmp_xy
    `Boolean'               has_weights
    `Boolean'               save_alphas
    `Vector'                b

    has_weights = asarray(f.extra, "has_weights")
    save_alphas = args()>=4 & alphas!=J(0,0,.)
    assert(has_weights==0 | has_weights==1) // bugbug: remove
    if (save_alphas) assert(cols(y)==1)
   
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
                if (save_alphas) {
                    b = invsym(tmp_xx) * tmp_xy
                    alphas[i, .] = tmp_ymeans - tmp_xmeans * b, b'
                    tmp_xbd = (tmp_x :- tmp_xmeans) * b :+ tmp_ymeans
                }
                else {
                    tmp_xbd = (tmp_x :- tmp_xmeans) * (invsym(tmp_xx) * tmp_xy) :+ tmp_ymeans
                }
            }
            else {
                tmp_xx = quadcross(tmp_x, tmp_w, tmp_x)
                // Bugbug: PRECOMPUTE! invsym(tmp_xx) !!!!!!!!!!!
                tmp_xy = quadcross(tmp_x, tmp_w, tmp_y)
                if (save_alphas) {
                    b = invsym(tmp_xx) * tmp_xy
                    alphas[i, .] = b'
                    tmp_xbd = tmp_x * b
                }
                else {
                    tmp_xbd = tmp_x * (invsym(tmp_xx) * tmp_xy)
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
end

