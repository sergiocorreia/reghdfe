// --------------------------------------------------------------------------
// FixedEffects main class
// --------------------------------------------------------------------------

mata:

class FixedEffects
{
    // Factors
    `Integer'               G                   // Number of sets of FEs
    `Integer'               N                   // number of obs
    `Factors'               factors
    `Vector'                sample
    `Varlist'               absvars
    `Varlist'               ivars
    `Varlist'               cvars
    `Boolean'               has_intercept
    `RowVector'             intercepts
    `RowVector'             num_slopes
    `Integer'               num_singletons
    `Boolean'               save_any_fe
    `Boolean'               save_all_fe
    `Varlist'               targets
    `RowVector'             save_fe

    // Optimization options
    `Real'                  tolerance
    `Integer'               maxiter
    `String'                transform           // Kaczmarz Cimmino Symmetric_kaczmarz (k c s)
    `String'                acceleration        // Acceleration method. None/No/Empty is none\
    `Integer'               accel_start         // Iteration where we start to accelerate // set it at 6? 2?3?
    `string'                slope_method
    `Boolean'               prune               // Whether to recursively prune degree-1 edges
    `Boolean'               abort               // Raise error if convergence failed?
    `Integer'               accel_freq          // Specific to Aitken's acceleration
    `Boolean'               storing_alphas      // 1 if we should compute the alphas/fes

    // Optimization objects
    `BipartiteGraph'        bg                  // Used when pruning 1-core vertices
    `Vector'                pruned_weight       // temp. weight for the factors that were pruned
    `Integer'               prune_g1            // Factor 1/2 in the bipartite subgraph that gets pruned
    `Integer'               prune_g2            // Factor 2/2 in the bipartite subgraph that gets pruned

    // Misc
    `Integer'               verbose
    `Boolean'               timeit

    // Weight-specific
    `Boolean'               has_weights
    `Variabe'               weight              // unsorted weight
    `String'                weight_var          // Weighting variable
    `String'                weight_type         // Weight type (pw, fw, etc)

    // Absorbed degrees-of-freedom computations
    `Integer'               G_extended          // Number of intercepts plus slopes
    `Integer'               df_a_redundant      // e(mobility)
    `Integer'               df_a_initial
    `Integer'               df_a                // df_a_inital - df_a_redundant
    `Vector'                doflist_M
    `Vector'                doflist_K
    `Vector'                doflist_M_is_exact
    `Vector'                doflist_M_is_nested
    `Vector'                is_slope
    `Integer'               df_a_nested // Redundant due to bein nested; used for: r2_a r2_a_within rmse

    // VCE and cluster variables
    `String'                vcetype
    `Integer'               num_clusters
    `Varlist'               clustervars
    `Varlist'               base_clustervars
    `String'                vceextra

    // Regression-specific
    `String'                original_varlist    // y x1 x2 (x3 x4 = z1 z2 z3)
    `String'                varlist             // y x1 x2 x3 x4 z1 z2 z3
    `String'                original_depvar
    `String'                original_indepvars
    `String'                original_endogvars
    `String'                original_instruments
    `String'                depvar              // y
    `String'                indepvars           // x1 x2
    `String'                endogvars           // x3 x4
    `String'                instruments         // z1 z2 z3
    
    `Boolean'               drop_singletons
    `String'                absorb              // contents of absorb()
    `String'                select_if           // If condition
    `String'                select_in           // In condition
    `String'                model               // ols, iv
    `String'                summarize_stats
    `Boolean'               summarize_quietly
    `StringRowVector'       dofadjustments // firstpair pairwise cluster continuous
    `Varname'               groupvar
    `String'                residuals
    `String'                diopts

    // Output
    `String'                cmdline
    `String'                subcmd
    `String'                title
    `Boolean'               converged
    `Integer'               iteration_count // e(ic)
    `Varlist'               extended_absvars
    `String'                notes
    `String'                equation_d
    `Integer'               df_r
    `Integer'               df_m
    `Integer'               N_clust
    `Integer'               N_clust_list
    `Real'                  rss
    `Real'                  rmse
    `Real'                  F
    `Real'                  tss
    `Real'                  tss_within
    `Real'                  sumweights
    `Real'                  r2
    `Real'                  r2_within
    `Real'                  r2_a
    `Real'                  r2_a_within
    `Real'                  ll
    `Real'                  ll_0

    // Methods
    `Void'                  new()
    `Matrix'                partial_out()
    `Void'                  _partial_out()
    `Variables'             project_one_fe()
    `Void'                  prune_1core()
    `Void'                  _expand_1core()
    `Void'                  estimate_dof()
    `Void'                  save_touse()
    `Void'                  store_alphas()
    `Void'                  save_variable()
    `Void'                  post_footnote()
    `Void'                  post()
}    


// Set default value of properties
`Void' FixedEffects::new()
{
    num_singletons = .
    sample = J(0, 1, .)
    weight = 1 // set to 1 so cross(x, S.weight, y)==cross(x, y)

    verbose = 0

    // Optimization defaults
    slope_method = "invsym"
    maxiter = 1e4
    tolerance = 1e-8
    transform = "symmetric_kaczmarz"
    acceleration = "conjugate_gradient"
    accel_start = 6
    
    prune = 1
    converged = 0
    abort = 1
    storing_alphas = 0

    // Specific to Aitken:
    accel_freq = 3
}


`Variables' FixedEffects::partial_out(`Anything' data, | `Boolean' save_tss)
{
    // -data- is either a varlist or a matrix
    `Variables'             y
    `Varlist'               vars
    `RowVector'             idx
    `Integer'               i

    if (eltype(data) == "string") {
        vars = tokens(invtokens(data))
        if (verbose > 0) printf("\n{txt} ## Partialling out %g variables: {res}%s{txt}\n\n", cols(vars), invtokens(vars))
        if (verbose > 0) printf("{txt}    - Loading variables into Mata\n")
        y = st_data(sample, invtokens(vars))

        if (cols(y) < cols(vars)) {
            printf("{err}(some columns were dropped due to a bug/quirk in st_data() (%g cols created instead of %g for %s); running slower workaround)\n", cols(y), cols(vars), invtokens(vars))
            assert(0) // when does this happen?
        }
        else if (cols(y) > cols(vars)) {
            printf("{err}(some empty columns were added due to a bug/quirk in {bf:st_data()}; %g cols created instead of %g for {it:%s}; running slower workaround)\n", cols(y), cols(vars), invtokens(vars))
            y = J(rows(y), 0, .)
            for (i=1; i<=cols(vars); i++) {
                y = y, st_data(sample, vars[i])
            }            
        }

        assert_msg(cols(y)==cols(vars), "st_data() constructed more columns than expected")

        idx = `selectindex'(colmax(abs(colminmax(y))) :== 0)
        if (cols(idx) & verbose>-1) {
            printf("{err}WARNING: after demeaning, some variables are just zeros: %s\n", invtokens(vars[idx]))
        }

        _partial_out(y, save_tss)
    }
    else {
        if (verbose > 0) printf("\n{txt} ## Partialling out %g variables\n\n", cols(data))
        _partial_out(y=data, save_tss)
    }
    //if (verbose==0) printf("\n")
    if (verbose==0) printf("{txt}(converged in %s iterations)\n", strofreal(iteration_count))
    return(y)
}


`Void' FixedEffects::store_alphas(`Varname' d_varname)
{
    `Integer'               g, i, j
    `StringRowVector'       varlabel
    `Variable'              d
    `RowVector'             tmp_stdev

    if (verbose > 0) printf("\n{txt} ## Storing estimated fixed effects\n\n")

    // Load -d- variable
    if (verbose > 0) printf("{txt}    - Loading d = e(depvar) - xb - e(resid)\n")
    d = st_data(sample, d_varname)
    assert(!missing(d))

    // Create empty alphas
    if (verbose > 0) printf("{txt}    - Initializing alphas\n")
    for (g=j=1; g<=G; g++) {
        if (!save_fe[g]) continue
        asarray(factors[g].extra, "alphas", J(factors[g].num_levels, intercepts[g] + num_slopes[g], 0))
        asarray(factors[g].extra, "tmp_alphas", J(factors[g].num_levels, intercepts[g] + num_slopes[g], 0))
    }

    // Fill out alphas
    if (verbose > 0) printf("{txt}    - Computing alphas\n")
    storing_alphas = 1
    d = accelerate_sd(this, d, &transform_kaczmarz())
    storing_alphas = 0

    if (verbose > 0) printf("{txt}    - SSR of d wrt FEs: %g\n", quadcross(d,d))

    // Store alphas in dataset
    if (verbose > 0) printf("{txt}    - Creating varlabels\n")
    for (g=j=1; g<=G; g++) {
        if (!save_fe[g]) {
            j = j + intercepts[g] + num_slopes[g]
            continue
        }
        varlabel = J(1, intercepts[g] + num_slopes[g], "")
        for (i=1; i<=cols(varlabel); i++) {
            varlabel[i] = sprintf("[FE] %s", extended_absvars[j])
            j++
        }

        if (num_slopes[g]) {
            if (verbose > 0) printf("{txt}    - Recovering unstandardized variables\n")
            tmp_stdev = asarray(factors[g].extra, "x_stdevs")
            if (intercepts[g]) tmp_stdev = 1, tmp_stdev

            // We need to *divide* the coefs by the stdev, not multiply!
            asarray(factors[g].extra, "alphas",
                asarray(factors[g].extra, "alphas") :/ tmp_stdev
            )
        }

        if (verbose > 0) printf("{txt}    - Storing alphas in dataset\n")
        save_variable(targets[g], asarray(factors[g].extra, "alphas")[factors[g].levels, .], varlabel)
        asarray(factors[g].extra, "alphas", .)
        asarray(factors[g].extra, "tmp_alphas", .)
    }
}


`Void' FixedEffects::_partial_out(`Variables' y, | `Boolean' save_tss)
{
    `RowVector'             stdevs, needs_zeroing
    `FunctionP'             funct_transform, func_accel // transform
    `Real'                  y_mean

    if (args()<2 | save_tss==.) save_tss = 0

    // Solver Warnings
    if (transform=="kaczmarz" & acceleration=="conjugate_gradient") {
        printf("{err}(WARNING: convergence is {bf:unlikely} with transform=kaczmarz and accel=CG)\n")
    }

    // Load transform pointer
    if (transform=="cimmino") funct_transform = &transform_cimmino()
    if (transform=="kaczmarz") funct_transform = &transform_kaczmarz()
    if (transform=="symmetric_kaczmarz") funct_transform = &transform_sym_kaczmarz()
    if (transform=="random_kaczmarz") funct_transform = &transform_rand_kaczmarz()

    // Pointer to acceleration routine
    if (acceleration=="test") func_accel = &accelerate_test()
    if (acceleration=="none") func_accel = &accelerate_none()
    if (acceleration=="conjugate_gradient") func_accel = &accelerate_cg()
    if (acceleration=="steepest_descent") func_accel = &accelerate_sd()
    if (acceleration=="aitken") func_accel = &accelerate_aitken()
    if (acceleration=="hybrid") func_accel = &accelerate_hybrid()

    // Shortcut for trivial case (1 FE)
    if (G==1) func_accel = &accelerate_none()

    // Compute TSS of depvar
    if (save_tss & tss==.) {
        if (has_intercept) {
            y_mean = mean(y[., 1], weight)
            tss = crossdev(y[., 1], y_mean, weight, y[., 1], y_mean) // Sum of w[i] * (y[i]-y_mean) ^ 2
        }
        else {
            tss = cross(y[., 1], weight, y[., 1]) // Sum of w[i] * y[i] ^ 2
        }
        if (weight_type=="aweight" | weight_type=="pweight") tss = tss * rows(y) / sum(weight)
    }

    // Standardize variables
    if (verbose > 0) printf("{txt}    - Standardizing variables\n")
    stdevs = reghdfe_standardize(y)

    // Solve
    if (verbose>0) printf("{txt}    - Running solver (acceleration={res}%s{txt}, transform={res}%s{txt} tol={res}%-1.0e{txt})\n", acceleration, transform, tolerance)
    if (verbose==1) printf("{txt}    - Iterating:")
    if (verbose>1) printf("{txt}      ")
    converged = 0
    iteration_count = 0
    y = (*func_accel)(this, y, funct_transform) :* stdevs
    // converged gets updated by check_convergence()
    
    if (prune) {
        _expand_1core(y)
    }

    // Standardizing makes it hard to detect values that are perfectly collinear with the absvars
    // in which case they should be 0.00 but they end up as 1e-16
    // EG: reghdfe price ibn.foreign , absorb(foreign)

    // This will edit to zero entire columns where *ALL* values are very close to zero
    needs_zeroing = colmax(abs(colminmax(y))) :< 1e-8 // chose something relatively close to epsilon() ~ 1e-16
    needs_zeroing = `selectindex'(needs_zeroing)
    if (cols(needs_zeroing)) {
        y[., needs_zeroing] = J(rows(y), cols(needs_zeroing), 0)
    }
}


`Variables' FixedEffects::project_one_fe(`Variables' y, `Integer' g)
{
    `Factor'                f
    `Boolean'               store_these_alphas
    `Matrix'                alphas, proj_y

    // Cons+K+W, Cons+K, K+W, K, Cons+W, Cons = 6 variants

    f = factors[g]
    store_these_alphas = storing_alphas & save_fe[g]
    if (store_these_alphas) assert(cols(y)==1)

    if (num_slopes[g]==0) {
        if (store_these_alphas) {
            alphas = panelmean(f.sort(y), f)
            asarray(factors[g].extra, "tmp_alphas", alphas)
            return(alphas[f.levels, .])
        }
        else {
            return(panelmean(f.sort(y), f)[f.levels, .])
        }
    }
    else {
        // This includes both cases, with and w/out intercept (## and #)
        if (store_these_alphas) {
            alphas = J(f.num_levels, intercepts[g] + num_slopes[g], .)
            proj_y = panelsolve_invsym(f.sort(y), f, intercepts[g], alphas)
            asarray(factors[g].extra, "tmp_alphas", alphas)
            return(proj_y)
        }
        else {
            return(panelsolve_invsym(f.sort(y), f, intercepts[g]))
        }
    }
}


`Void' FixedEffects::estimate_dof()
{
    `Boolean'               has_int
    `Integer'               g, h                        // index FEs (1..G)
    `Integer'               num_intercepts              // Number of absvars with an intercept term
    `Integer'               i_cluster, i_intercept, j_intercept
    `Integer'               i                           // index 1..G_extended
    `Integer'               j
    `Integer'               bg_verbose                  // verbose level when calling BipartiteGraph()
    `Integer'               m                           // Mobility groups between a specific pair of FEs
    `RowVector'             SubGs
    `RowVector'             offsets, idx, zeros, results
    `Matrix'                tmp
    `Variables'             data
    `DataCol'               cluster_data
    `String'                absvar, clustervar
    `Factor'                F
    `BipartiteGraph'        BG
    `Integer'               pair_count
    
    if (verbose > 0) printf("\n{txt} ## Estimating degrees-of-freedom absorbed by the fixed effects\n\n")

    // Count all FE intercepts and slopes
    SubGs = intercepts + num_slopes
    G_extended = sum(SubGs)
    num_intercepts = sum(intercepts)
    offsets = runningsum(SubGs) - SubGs :+ 1 // start of each FE within the extended list
    idx = `selectindex'(intercepts) // Select all FEs with intercepts
    if (verbose > 0) printf("{txt}    - there are %f fixed intercepts and slopes in the %f absvars\n", G_extended, G)

    // Initialize result vectors and scalars
    doflist_M_is_exact = J(1, G_extended, 0)
    doflist_M_is_nested = J(1, G_extended, 0)
    df_a_nested = 0

    // (1) M will hold the redundant coefs for each extended absvar (G_extended, not G)
    doflist_M = J(1, G_extended, 0)
    assert(0 <= num_clusters & num_clusters <= 10)
    if (num_clusters > 0 & anyof(dofadjustments, "clusters")) {
        
        // (2) (Intercept-Only) Look for absvars that are clustervars
        for (i_intercept=1; i_intercept<=length(idx); i_intercept++) {
            g = idx[i_intercept]
            i = offsets[g]
            absvar = invtokens(tokens(ivars[g]), "#")
            if (anyof(clustervars, absvar)) {
                doflist_M[i] = factors[g].num_levels
                df_a_nested = df_a_nested + doflist_M[i]
                doflist_M_is_exact[i] = doflist_M_is_nested[i] = 1
                idx[i_intercept] = 0
                if (verbose > 0) printf("{txt} - categorical variable {res}%s{txt} is also a cluster variable, so it doesn't reduce DoF\n", absvar)
            }
        }
        idx = select(idx, idx)

        // (3) (Intercept-Only) Look for absvars that are nested within a clustervar
        for (i_cluster=1; i_cluster<= num_clusters; i_cluster++) {
            cluster_data = .
            if (!length(idx)) break // no more absvars to process
            for (i_intercept=1; i_intercept<=length(idx); i_intercept++) {

                g = idx[i_intercept]
                i = offsets[g]
                absvar = invtokens(tokens(ivars[g]), "#")
                clustervar = clustervars[i_cluster]
                if (doflist_M_is_exact[i]) continue // nothing to do

                if (cluster_data == .) {
                    if (strpos(clustervar, "#")) {
                        clustervar = base_clustervars[i_cluster]
                        F = factor(clustervar, sample, ., "", 0, 0, ., 0)
                        cluster_data = F.levels
                        F = Factor() // clear
                    }
                    else {
                        cluster_data = __fload_data(clustervar, sample)
                    }
                }

                if (factors[g].nested_within(cluster_data)) {
                    doflist_M[i] = factors[g].num_levels
                    doflist_M_is_exact[i] = doflist_M_is_nested[i] = 1
                    df_a_nested = df_a_nested + doflist_M[i]
                    idx[i_intercept] = 0
                    if (verbose > 0) printf("{txt} - categorical variable {res}%s{txt} is nested within a cluster variable, so it doesn't reduce DoF\n", absvar)
                }
            }
            idx = select(idx, idx)
        }
        cluster_data = . // save memory
    } // end of the two cluster checks (absvar is clustervar; absvar is nested within clustervar)


    // (4) (Intercept-Only) Every intercept but the first has at least one redundant coef.
    if (length(idx) > 1) {
        if (verbose > 0) printf("{txt}    - there is at least one redundant coef. for every set of FE intercepts after the first one\n")
        doflist_M[offsets[idx[2..num_intercepts]]] = J(1, num_intercepts-1, 1) // Set DoF loss of all intercepts but the first one to 1
    }

    // (5) (Intercept-only) Mobility group algorithm
    // Excluding those already solved, the first absvar is exact

    if (length(idx)) {
        i = idx[1]
        doflist_M_is_exact[i] = 1
    }

    // Compute number of dijsoint subgraphs / mobility groups for each pair of remaining FEs
    if (anyof(dofadjustments, "firstpair") | anyof(dofadjustments, "pairwise")) {
        BG = BipartiteGraph()
        bg_verbose = max(( verbose - 1 , 0 ))
        pair_count = 0

        for (i_intercept=1; i_intercept<=length(idx)-1; i_intercept++) {
            for (j_intercept=i_intercept+1; j_intercept<=length(idx); j_intercept++) {
                g = idx[i_intercept]
                h = idx[j_intercept]
                i = offsets[h]
                BG.init(&factors[g], &factors[h], bg_verbose)
                m = BG.init_zigzag()
                ++pair_count
                if (verbose > 0) printf("{txt}    - mobility groups between FE intercepts #%f and #%f: {res}%f\n", g, h, m)
                doflist_M[i] = max(( doflist_M[i] , m ))
                if (j_intercept==2) doflist_M_is_exact[i] = 1
                if (pair_count & anyof(dofadjustments, "firstpair")) break
            }
            if (pair_count & anyof(dofadjustments, "firstpair")) break
        }
        BG = BipartiteGraph() // clear
    }
    // TODO: add group3hdfe

    // (6) See if cvars are zero (w/out intercept) or just constant (w/intercept)
    if (anyof(dofadjustments, "continuous")) {
        for (i=g=1; g<=G; g++) {
            // If model has intercept, redundant cvars are those that are CONSTANT
            // Without intercept, a cvar has to be zero within a FE for it to be redundant
            // Since S.fes[g].x are already demeaned IF they have intercept, we don't have to worry about the two cases
            has_int = intercepts[g]
            if (has_int) i++
            if (!num_slopes[g]) continue

            data = asarray(factors[g].extra, "x")
            assert(num_slopes[g]==cols(data))
            results = J(1, cols(data), 0)
            // float(1.1) - 1 == 2.384e-08 , so let's pick something bigger, 1e-6
            zeros = J(1, cols(data), 1e-6)
            // This can be speed up by moving the -if- outside the -for-
            for (j = 1; j <= factors[g].num_levels; j++) {
                tmp = colminmax(panelsubmatrix(data, j, factors[g].info))
                if (has_int) {
                    results = results + ((tmp[2, .] - tmp[1, .]) :<= zeros)
                }
                else {
                    results = results + (colsum(abs(tmp)) :<= zeros)
                }
            }
            data = .
            if (sum(results)) {
                if (has_int  & verbose) printf("{txt}    - the slopes in the FE #%f are constant for {res}%f{txt} levels, which don't reduce DoF\n", g, sum(results))
                if (!has_int & verbose) printf("{txt}    - the slopes in the FE #%f are zero for {res}%f{txt} levels, which don't reduce DoF\n", g, sum(results))
                doflist_M[i..i+num_slopes[g]-1] = results
            }
            i = i + num_slopes[g]
        }
    }

    // Store results (besides doflist_..., etc.)
    doflist_K = J(1, G_extended, .)
    for (g=1; g<=G; g++) {
        i = offsets[g]
        j = g==G ? G_extended : offsets[g+1]
        doflist_K[i..j] = J(1, j-i+1, factors[g].num_levels)
    }
    df_a_initial = sum(doflist_K)
    df_a_redundant = sum(doflist_M)
    df_a = df_a_initial - df_a_redundant
}



`Void' FixedEffects::prune_1core()
{
    // Note that we can't prune degree-2 nodes, to keep the graph bipartite
    `Integer'               i, j, g
    `Vector'                subgraph_id
    
    `Vector'                idx
    `RowVector'             i_prune

    assert(G==2) // bugbug remove?

    // Abort if the user set HDFE.prune = 0
    if (!prune) return
    
    // Pick two factors, and check if we really benefit from pruning
    prune = 0
    i_prune = J(1, 2, 0)
    for (g=i=1; g<=2; g++) {
        //if (intercepts[g] & !num_slopes[g] & factors[g].num_levels>100) {
        if (intercepts[g] & !num_slopes[g]) {
            i_prune[i++] = g // increments at the end
            if (i > 2) { // success!
                prune = 1
                break
            }
        }
    }

    if (!prune) return

    // for speed, the factor with more levels goes first
    i = i_prune[1]
    j = i_prune[2]
    //if (factors[i].num_levels < factors[j].num_levels) swap(i, j) // bugbug uncomment it!
    prune_g1 = i
    prune_g2 = j

    bg = BipartiteGraph()
    bg.init(&factors[prune_g1], &factors[prune_g2], verbose)
    bg.init_zigzag(1) // 1 => save subgraphs into bg.subgraph_id
    bg.compute_cores()
    bg.prune_1core(weight)
}

// bugbug fix or remove this fn altogether
`Void' FixedEffects::_expand_1core(`Variables' y)
{
    y = bg.expand_1core(y)
}


`Void' FixedEffects::save_touse(`Varname' touse, | `Boolean' replace)
{
    `Integer'               idx
    `Vector'                mask


    if (verbose > 0) printf("\n{txt} ## Saving e(sample)\n")
    mask = J(st_nobs(), 1, 0)
    mask[sample] = J(rows(sample), 1, 1)
    
    if (args()<2 | replace==.) {
        // Generate
        idx = st_addvar("byte", touse)
        st_store(., idx, mask)
    }
    else {
        // Replace
        st_store(., touse, mask)
    }
}


`Void' FixedEffects::save_variable(`Varname' varname,
                                   `Variables' data,
                                 | `Varlist' varlabel)
{
    `RowVector'               idx
    `Integer'               i
    `Vector'                mask
    idx = st_addvar("double", tokens(varname))
    st_store(sample, idx, data)
    if (args()>=3 & varlabel!="") {
        for (i=1; i<=cols(data); i++) {
            st_varlabel(idx[i], varlabel[i])
        }
    }

}



`Void' FixedEffects::post_footnote()
{
    `Matrix'                table
    `StringVector'          rowstripe
    `StringRowVector'       colstripe
    `String'                text

    assert(!missing(G))
    st_numscalar("e(N_hdfe)", G)
    st_numscalar("e(N_hdfe_extended)", G_extended)
    st_numscalar("e(df_a)", df_a)
    st_numscalar("e(df_a_initial)", df_a_initial)
    st_numscalar("e(df_a_redundant)", df_a_redundant)
    st_numscalar("e(df_a_nested)", df_a_nested)

    st_global("e(absvars)", invtokens(absvars))
    text = invtokens(extended_absvars)
    text = subinstr(text, "1.", "")
    st_global("e(extended_absvars)", text)

    // Absorbed degrees-of-freedom table
    table = (doflist_K \ doflist_M \ (doflist_K-doflist_M) \ !doflist_M_is_exact \ doflist_M_is_nested)'
    rowstripe = extended_absvars'
    rowstripe = J(rows(table), 1, "") , extended_absvars' // add equation col
    colstripe = "Categories" \ "Redundant" \ "Num Coefs" \ "Exact?" \ "Nested?" // colstripe cannot have dots on Stata 12 or earlier
    colstripe = J(cols(table), 1, "") , colstripe // add equation col
    st_matrix("e(dof_table)", table)
    st_matrixrowstripe("e(dof_table)", rowstripe)
    st_matrixcolstripe("e(dof_table)", colstripe)

    st_numscalar("e(ic)", iteration_count)
    st_numscalar("e(drop_singletons)", drop_singletons)
    st_numscalar("e(num_singletons)", num_singletons)
}


`Void' FixedEffects::post(`Options' options)
{
    `String'        text
    `Integer'       i

    post_footnote()

    // ---- constants -------------------------------------------------------

    st_global("e(predict)", "reghdfe_p")
    st_global("e(estat_cmd)", "reghdfe_estat")
    st_global("e(footnote)", "reghdfe_footnote")
    st_global("e(marginsok)", "")
    st_numscalar("e(df_m)", df_m)


    assert(title != "")
    text = sprintf("HDFE %s", title)
    st_global("e(title)", text)
    
    text = sprintf("Absorbing %g HDFE %s", G, plural(G, "group"))
    st_global("e(title2)", text)
    
    st_global("e(model)", model)
    st_global("e(cmdline)", cmdline)

    st_numscalar("e(tss)", tss)
    st_numscalar("e(tss_within)", tss_within)
    st_numscalar("e(rss)", rss)
    st_numscalar("e(mss)", tss - rss)
    st_numscalar("e(rmse)", rmse)
    st_numscalar("e(F)", F)

    st_numscalar("e(ll)", ll)
    st_numscalar("e(ll_0)", ll_0)

    st_numscalar("e(r2)", r2)
    st_numscalar("e(r2_within)", r2_within)
    st_numscalar("e(r2_a)", r2_a)
    st_numscalar("e(r2_a_within)", r2_a_within)
    
    if (!missing(N_clust)) {
        st_numscalar("e(N_clust)", N_clust)
        for (i=1; i<=num_clusters; i++) {
            text = sprintf("e(N_clust%g)", i)
            st_numscalar(text, N_clust_list[i])
        }
        text = "Statistics robust to heteroskedasticity"
        st_global("e(title3)", text)
    }

    if (!missing(sumweights)) st_numscalar("e(sumweights)", sumweights)     


    // ---- .options properties ---------------------------------------------

    st_global("e(depvar)", depvar)
    st_global("e(indepvars)", indepvars)
    st_global("e(endogvars)", endogvars)
    st_global("e(instruments)", instruments)

    if (!missing(N_clust)) {
        st_numscalar("e(N_clustervars)", num_clusters)
        st_global("e(clustvar)", invtokens(clustervars))
        for (i=1; i<=num_clusters; i++) {
            text = sprintf("e(clustvar%g)", i)
            st_global(text, clustervars[i])
        }
    }

    if (residuals != "") {
        st_global("e(resid)", residuals)
    }

    // Stata uses e(vcetype) for the SE column headers
    // In the default option, leave it empty.
    // In the cluster and robust options, set it as "Robust"
    text = strproper(vcetype)
    if (text=="Cluster") text = "Robust"
    if (text=="Unadjusted") text = ""
    assert(anyof( ("", "Robust", "Jackknife", "Bootstrap") , text))
    if (text!="") st_global("e(vcetype)", text)

    text = vcetype
    if (text=="unadjusted") text = "ols"
    st_global("e(vce)", text)

    // Weights
    if (weight_type != "") {
        st_global("e(wexp)", "= " + weight_var)
        st_global("e(wtype)", weight_type)
    }
}

end
