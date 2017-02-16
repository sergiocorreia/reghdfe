// --------------------------------------------------------------------------
// FixedEffects main class
// --------------------------------------------------------------------------

mata:

class FixedEffects
{
    // Factors
    `Integer'               G // number of FE levels (2 for two-way-fe, etc.)
    `Integer'               N // number of obs
    `Boolean'               has_intercept
    `Boolean'               save_any_fe
    `Boolean'               save_all_fe
    `Varlist'               absvars
    `Varlist'               ivars
    `Varlist'               cvars
    `Varlist'               targets
    `RowVector'             intercepts
    `RowVector'             num_slopes
    `Factors'               factors
    `Vector'                sample
    `Integer'               num_singletons

    // Weight-specific
    `Boolean'               has_weights
    `String'                weighttype
    `Varname'               weightvar
    `Variabe'               weight // unsorted weight

    // Misc
    `Integer'               verbose
    `Boolean'               timeit

    // Optimization options
    `Real'                  tolerance
    `Integer'               maxiter
    `String'                transform           // Kaczmarz Cimmino Symmetric_kaczmarz (k c s)
    `String'                acceleration        // Acceleration method. None/No/Empty is none\
    `Integer'               accel_start         // Iteration where we start to accelerate // set it at 6? 2?3?
    `string'                slope_method
    `Boolean'               prune               // Whether to recursively prune degree-1 edges
    `Boolean'               abort               // Raise error if convergence failed?

    // Specific to Aitken's acceleration
    `Integer'               accel_freq

    // Used when pruning 1-core vertices
    `BipartiteGraph'        bg


    `Vector'                pruned_weight       // temp. weight for the factors that were pruned
    `Integer'               prune_g1            // Factor 1/2 in the bipartite subgraph that gets pruned
    `Integer'               prune_g2            // Factor 2/2 in the bipartite subgraph that gets pruned


    // Additional options not used here
    `Options'               options

    // Results object
    `Output'                output

    // Methods
    `Void'                  new()
    `Matrix'                partial_out()
    `Void'                  _partial_out()
    `Variables'             project_one_fe()
    `Void'                  prune_1core()
    `Void'                  _expand_1core()
    `Void'                  estimate_dof()
    `Void'                  save_touse()
    `Void'                  save_variable()
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
    output.converged = 0
    prune = 1
    abort = 1

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
    if (verbose==0) printf("{txt}(converged in %s iterations)\n", strofreal(output.iteration_count))
    return(y)
}


`Void' FixedEffects::_partial_out(`Variables' y, | `Boolean' save_tss)
{
    `RowVector'             stdevs, needs_zeroing
    `FunctionP'             funct_transform, func_accel // transform
    `Real'                y_mean

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
    if (save_tss & output.tss==.) {
        if (has_intercept) {
            y_mean = mean(y[., 1], weight)
            output.tss = crossdev(y[., 1], y_mean, weight, y[., 1], y_mean) // Sum of w[i] * (y[i]-y_mean) ^ 2
        }
        else {
            output.tss = cross(y[., 1], weight, y[., 1]) // Sum of w[i] * y[i] ^ 2
        }
        if (weighttype=="aweight" | weighttype=="pweight") output.tss = output.tss * rows(y) / sum(weight)
    }

    // Standardize variables
    if (verbose > 0) printf("{txt}    - Standardizing variables\n")
    stdevs = reghdfe_standardize(y)

    // Solve
    if (verbose>0) printf("{txt}    - Running solver (acceleration={res}%s{txt}, transform={res}%s{txt} tol={res}%-1.0e{txt})\n", acceleration, transform, tolerance)
    if (verbose==1) printf("{txt}    - Iterating:")
    if (verbose>1) printf("{txt}      ")
    output.converged = 0
    output.iteration_count = 0
    y = (*func_accel)(this, y, funct_transform) :* stdevs
    // output.converged gets updated by check_convergence()
    
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

    if (weight==1) weight = J(rows(y), 1, .) // bugbug remove this!??!
}


`Variables' FixedEffects::project_one_fe(`Variables' y, `Integer' g)
{
    `Factor'                f

    // Cons+K+W, Cons+K, K+W, K, Cons+W, Cons = 6 variants

    f = factors[g]

    if (num_slopes[g]==0) {
        // bugbug see if we can remove the .]
        return(panelmean(f.sort(y), f)[f.levels, .])
    }
    else if (intercepts[g]) {
        return(panelsolve_invsym(f.sort(y), f, intercepts[g]))
    }
    else {
        return(panelsolve_invsym(f.sort(y), f, intercepts[g]))
    }
}


`Void' FixedEffects::estimate_dof()
{
    `Boolean'               has_int
    `Integer'               SuperG                      // Number of "extended absvars" (intercepts and slopes in the G absvars)
    `Integer'               g, h                        // index FEs (1..G)
    `Integer'               num_intercepts              // Number of absvars with an intercept term
    `Integer'               i_cluster, i_intercept, j_intercept
    `Integer'               i                           // index 1..SuperG
    `Integer'               j
    `Integer'               bg_verbose                  // verbose level when calling BipartiteGraph()
    `Integer'               m                           // Mobility groups between a specific pair of FEs
    `Integer'               M_due_to_nested             // Used for r2_a r2_a_within rmse
    `RowVector'             SubGs, M, M_is_exact, M_is_nested
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
    SuperG = sum(SubGs)
    num_intercepts = sum(intercepts)
    offsets = runningsum(SubGs) - SubGs :+ 1 // start of each FE within the extended list
    idx = `selectindex'(intercepts) // Select all FEs with intercepts
    if (verbose > 0) printf("{txt}    - there are %f fixed intercepts and slopes in the %f absvars\n", SuperG, G)

    // Initialize result vectors and scalars
    M_is_exact = J(1, SuperG, 0)
    M_is_nested = J(1, SuperG, 0)
    M_due_to_nested = 0

    // (1) M will hold the redundant coefs for each extended absvar (G_extended, not G)
    M = J(1, SuperG, 0)

    assert(0 <= options.num_clusters & options.num_clusters <= 10)
    if (options.num_clusters > 0 & anyof(options.dofadjustments, "clusters")) {
        
        // (2) (Intercept-Only) Look for absvars that are clustervars
        for (i_intercept=1; i_intercept<=length(idx); i_intercept++) {
            g = idx[i_intercept]
            i = offsets[g]
            absvar = invtokens(tokens(ivars[g]), "#")
            if (anyof(options.clustervars, absvar)) {
                M[i] = factors[g].num_levels
                M_due_to_nested = M_due_to_nested + M[i]
                M_is_exact[i] = M_is_nested[i] = 1
                idx[i_intercept] = 0
                if (verbose > 0) printf("{txt} - categorical variable {res}%s{txt} is also a cluster variable, so it doesn't reduce DoF\n", absvar)
            }
        }
        idx = select(idx, idx)

        // (3) (Intercept-Only) Look for absvars that are nested within a clustervar
        for (i_cluster=1; i_cluster<= options.num_clusters; i_cluster++) {
            cluster_data = .
            if (!length(idx)) break // no more absvars to process
            for (i_intercept=1; i_intercept<=length(idx); i_intercept++) {

                g = idx[i_intercept]
                i = offsets[g]
                absvar = invtokens(tokens(ivars[g]), "#")
                clustervar = options.clustervars[i_cluster]
                if (M_is_exact[i]) continue // nothing to do

                if (cluster_data == .) {
                    if (strpos(clustervar, "#")) {
                        clustervar = options.base_clustervars[i_cluster]
                        F = factor(clustervar, sample, ., "", 0, 0, ., 0)
                        cluster_data = F.levels
                        F = Factor() // clear
                    }
                    else {
                        cluster_data = __fload_data(clustervar, sample)
                    }
                }

                if (factors[g].nested_within(cluster_data)) {
                    M[i] = factors[g].num_levels
                    M_is_exact[i] = M_is_nested[i] = 1
                    M_due_to_nested = M_due_to_nested + M[i]
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
        M[offsets[idx[2..num_intercepts]]] = J(1, num_intercepts-1, 1) // Set DoF loss of all intercepts but the first one to 1
    }


    // (5) (Intercept-only) Mobility group algorithm
    // Excluding those already solved, the first absvar is exact

    if (length(idx)) {
        i = idx[1]
        M_is_exact[i] = 1
    }

    // Compute number of dijsoint subgraphs / mobility groups for each pair of remaining FEs
    if (anyof(options.dofadjustments, "firstpair") | anyof(options.dofadjustments, "pairwise")) {
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
                M[i] = max(( M[i] , m ))
                if (j_intercept==2) M_is_exact[i] = 1
                if (pair_count & anyof(options.dofadjustments, "firstpair")) break
            }
            if (pair_count & anyof(options.dofadjustments, "firstpair")) break
        }
        BG = BipartiteGraph() // clear
    }
    // TODO: add group3hdfe

    // (6) See if cvars are zero (w/out intercept) or just constant (w/intercept)
    if (anyof(options.dofadjustments, "continuous")) {
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
                M[i..i+num_slopes[g]-1] = results
            }
            i = i + num_slopes[g]
        }
    }

    // Store results
    output.G_extended = SuperG
    output.doflist_M = M
    output.doflist_M_is_exact = M_is_exact
    output.doflist_M_is_nested = M_is_nested
    output.doflist_K = J(1, SuperG, .)
    for (g=1; g<=G; g++) {
        i = offsets[g]
        j = g==G ? SuperG : offsets[g+1]
        output.doflist_K[i..j] = J(1, j-i+1, factors[g].num_levels)
    }
    output.dof_M = sum(M)
    output.df_a = sum(output.doflist_K) - output.dof_M
    output.M_due_to_nested = M_due_to_nested
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
                                   `Variable' data,
                                 | `String' varlabel)
{
    `Integer'               idx
    `Vector'                mask
    idx = st_addvar("double", varname)
    st_store(sample, idx, data)
    if (args()>=3 & varlabel!="") st_varlabel(idx, varlabel)

}


end
