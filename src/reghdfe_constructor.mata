// --------------------------------------------------------------------------
// FixedEffects constructor (also precomputes factors)
// --------------------------------------------------------------------------

mata:

`FixedEffects' fixed_effects(`Varlist' absvars,
                           | `Varname' touse,
                             `String' weighttype,
                             `Varname' weightvar,
                             `Boolean' drop_singletons,
                             `Boolean' verbose,
                             class FixedEffects matrix S // S must be of orgtype matrix, else it cannot be optional
                             )
{
    `Varname'               absvar, cvars
    `Integer'               i, j, g, gg, remaining
    `Boolean'               use_sample
    `Vector'                idx
    `Integer'               spaces
    `Integer'               num_singletons_i
    `Variables'             cvar_data
    `Variable'              w
    `FactorPointer'         pf

    // Set default value of arguments
    if (args()<2) touse = ""
    if (args()<3) weighttype = ""
    if (args()<4) weightvar = ""
    if (args()<5 | drop_singletons==.) drop_singletons = 1
    if (args()<6 | verbose==.) verbose = 0
    if (args()<7) S = FixedEffects()

    S.verbose = verbose

    // Parse absvars
    if (S.verbose > 0) printf("\n{txt} ## Parsing absvars\n")
    stata(`"ms_parse_absvars "' + absvars)
    if (S.verbose > 0) stata("sreturn list")
    S.G = strtoreal(st_global("s(G)"))
    S.absvars = tokens(st_global("s(absvars)"))
    S.has_intercept = strtoreal(st_global("s(has_intercept)"))
    S.save_any_fe = strtoreal(st_global("s(save_any_fe)"))
    S.save_all_fe = strtoreal(st_global("s(save_all_fe)"))
    S.absvars = tokens(st_global("s(absvars)"))
    S.ivars = tokens(st_global("s(ivars)"))
    S.cvars = tokens(st_global("s(cvars)"))
    S.targets = tokens(st_global("s(targets)"))
    S.intercepts = strtoreal(tokens(st_global("s(intercepts)")))
    S.num_slopes = strtoreal(tokens(st_global("s(num_slopes)")))

    // TODO:
    // Allow verbose, timeit, all optimization options, keepsingleton
    // "s(options)"

    if (S.verbose > -1 & !S.has_intercept) printf("{txt}(warning: no intercepts terms in absorb(); regression lacks constant term)\n")


    // Store results (HDFE.options cannot access HDFE easily because its class is defined before the HDFE class)
    S.output.G = S.G
    S.output.absvars = S.absvars
    S.output.extended_absvars = tokens(st_global("s(extended_absvars)"))
    S.output.equation_d = st_global("s(equation_d)")
    S.output.tss = .

    assert(1<=S.G & S.G<=10)
    assert(S.G == cols(S.ivars))
    assert(S.G == cols(S.cvars))
    assert(S.G == cols(S.targets))
    assert(S.G == cols(S.intercepts))
    assert(S.G == cols(S.num_slopes))


    // Fill out object
    S.G = cols(S.absvars)
    S.factors = Factor(S.G)
    assert_msg(anyof(("", "fweight", "pweight", "aweight"), weighttype), "wrong weight type")
    S.weighttype = weighttype
    S.weightvar = weightvar

    S.sample = (touse=="") ? 1::st_nobs() : `selectindex'(st_data(., touse))

    if (drop_singletons) {
        S.num_singletons = 0
        num_singletons_i = 0
        if (weighttype=="fweight") {
            S.weight = st_data(S.sample, S.weightvar) // just to use it in F.drop_singletons()
        }
    }


    // (1) create the factors and remove singletons
    remaining = S.G
    i = 0
    if (S.verbose > 0) {
        printf("\n{txt} ## Initializing Mata object for %g fixed effects\n\n", S.G)
        spaces = max((0, max(strlen(S.absvars))-4))
        printf("{txt}   |  i | g | %s Name | Int? | #Slopes |    Obs.   |   Levels   | Sorted? | #Drop Singl. |\n", " " * spaces)
        printf("{txt}   |----|---|-%s------|------|---------|-----------|------------|---------|--------------|\n", "-" * spaces)
        displayflush()
    }

    while (remaining) {
        ++i
        g = 1 + mod(i-1, S.G)
        absvar = S.absvars[g]
        use_sample = (touse!="") | (i>1 & drop_singletons)
        
        if (S.verbose > 0) {
            printf("{txt}   | %2.0f | %1.0f | {res}%s{txt} | ", i, g, (spaces+5-strlen(absvar)) * " " + absvar)
            printf("{txt}  {%s}%1.0f{txt}  |    %1.0f    |", S.intercepts[g] ? "txt" : "err", S.intercepts[g], S.num_slopes[g])
            displayflush()
        }
        if (i<=S.G) {
            if (use_sample) assert_msg(rows(S.sample), "empty sample")
            S.factors[g] = factor(S.ivars[g], use_sample ? S.sample : .)
        }
        if (S.verbose > 0) {
            printf("{res}%10.0g{txt} | {res}%10.0g{txt} | %7.0f |", S.factors[g].num_obs, S.factors[g].num_levels, S.factors[g].is_sorted)
            displayflush()
        }
 
        if (drop_singletons) {
            
            if (weighttype=="fweight") {
                idx = S.factors[g].drop_singletons(S.weight)
            }
            else {
                idx = S.factors[g].drop_singletons()
            }

            num_singletons_i = rows(idx)
            S.num_singletons = S.num_singletons + num_singletons_i
            if (S.verbose > 0) {
                printf(" %10.0g   |", num_singletons_i)
                displayflush()
            }

            if (num_singletons_i==0) {
                --remaining
            }
            else {
                remaining = S.G - 1
                
                // sample[idx] = . // not allowed in Mata; instead, make 0 and then select()
                S.sample[idx] = J(rows(idx), 1, 0)
                S.sample = select(S.sample, S.sample)

                for (j=i-1; j>=max((1, i-remaining)); j--) {
                    gg = 1 + mod(j-1, S.G)
                    S.factors[gg].drop_obs(idx)
                    if (S.verbose > 0) printf("{res} .")
                }
            }
        }
        else {
            if (S.verbose > 0) printf("      n/a     |")
            --remaining
        }
        if (S.verbose > 0) printf("\n")
    }
    if (S.verbose > 0) printf("\n")

    if (S.factors[1].num_obs == 0) {
        if (S.verbose > -1) printf("{err}(no observations remaining after dropping singletons!)\n")
        return(S)
    }

    if ( drop_singletons & (S.num_singletons > 0) & (S.verbose > -1) ) {
        printf("{txt}(dropped %s singleton observations)\n", strofreal(S.num_singletons))
    }

    S.N = S.factors[1].num_obs // store number of obs.
    assert(S.N = S.factors[S.G].num_obs)


    // (2) run *.panelsetup() after the sample is defined
    if (S.verbose > 0) printf("\n{txt} ## Initializing panelsetup() for each fixed effect\n\n")
    for (g=1; g<=S.G; g++) {
        absvar = S.absvars[g]
        if (S.verbose > 0) printf("{txt}    - panelsetup({res}%s{txt})\n", absvar)
        S.factors[g].panelsetup()
    }


    // (3) load weight
    S.has_weights = (S.weighttype !="" & S.weightvar!="")
    for (g=1; g<=S.G; g++) {
        asarray(S.factors[g].extra, "has_weights", S.has_weights)
    }
    if (S.has_weights) {
        if (S.verbose > 0) printf("\n{txt} ## Loading and sorting weights\n\n")
        S.weight = st_data(S.sample, S.weightvar)
        if (S.verbose > 0) printf("{txt}    - loading %s weight from variable %s)\n", S.weighttype, S.weightvar)
        
        for (g=1; g<=S.G; g++) {
            if (S.verbose > 0) printf("{txt}    - sorting weight for factor {res}%s{txt}\n", S.absvars[g])
            pf = &(S.factors[g])
            w = (*pf).sort(S.weight)
            asarray((*pf).extra, "weights", w)
            asarray((*pf).extra, "weighted_counts", `panelsum'(w, (*pf).info))
        }
        w = .
    }


    // (4) prune edges of degree-1
    S.prune = 0 // bugbug
    if (S.prune) S.prune_1core()


    // (5) load cvars
    if (sum(S.num_slopes)) {
        if (S.verbose > 0) printf("\n{txt} ## Loading slope variables\n\n")
        for (g=1; g<=S.G; g++) {
            cvars = tokens(S.cvars[g])
            if (cols(cvars) > 0) {
                // Load, standardize, sort by factor, precompute (TODO), and store
                if (S.verbose > 0) printf("{txt}    - cvars({res}%s{txt})\n", invtokens(cvars))
                pf = &(S.factors[g])
                cvar_data = (*pf).sort(st_data(S.sample, cvars))
                (void) reghdfe_standardize(cvar_data)
                asarray((*pf).extra, "x", cvar_data)
                if (S.intercepts[g]) {
                    asarray((*pf).extra, "xmeans", panelmean(cvar_data, *pf))
                }
            }
        }
        cvar_data = .
    }

    return(S)
}

end
