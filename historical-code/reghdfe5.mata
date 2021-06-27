// --------------------------------------------------------------------------
// Mata Code: FE Estimator (REGHDFE)
// --------------------------------------------------------------------------
// - Project URL:   https://github.com/sergiocorreia/reghdfe
// - Dependency:    https://github.com/sergiocorreia/ftools

    *mata: mata set matastrict on
    mata: mata set mataoptimize on
    *mata: mata set matadebug off
    *mata: mata set matalnum off

// Include ftools -----------------------------------------------------------
    cap findfile "ftools.mata"
    if (_rc) {
        di as error "reghdfe requires the {bf:ftools} package, which is not installed"
        di as error `"    - install from {stata ssc install ftools:SSC}"'
        di as error `"    - install from {stata `"net install ftools, from("https://github.com/sergiocorreia/ftools/raw/master/src/")"':Github}"'
        exit 9
    }
    include "`r(fn)'"


// Custom types -------------------------------------------------------------
    loc FixedEffects        class FixedEffects scalar
    loc Factors             class Factor rowvector
    loc BipartiteGraph      class BipartiteGraph scalar
    loc FactorPointer       pointer(`Factor') scalar


// Versioning ---------------------------------------------------------------
	*ms_get_version reghdfe // from parsetools package
	*assert("`package_version'" != "")
    *mata: string scalar reghdfe_version() return("`package_version'")
    *mata: string scalar reghdfe_stata_version() return("`c(stata_version)'")
    *mata: string scalar reghdfe_joint_version() return("`package_version'|`c(stata_version)'")


// Includes -----------------------------------------------------------------
// Bipartite Graphs ---------------------------------------------------------
// - For simplicity, assume the graph represent (firm, ceo) pairs
// - TODO: Check when we don't need all these objects anymore and clean them up!

mata:

class BipartiteGraph
{
	// Computed by init()
	`Boolean'				verbose
	`Integer'				N		// Num. obs
	`Integer'				N1		// Num. levels of FE 1
	`Integer'				N2		// Num. levels of FE 2
	`Integer'				N12		// N1 + N2
	`FactorPointer'			PF1
	`FactorPointer'			PF2
	`Factor'				F12
	`Factor'				F12_1
	`Factor'				F12_2

	// Computed by init_zigzag()
	`Vector' 				queue
	`Vector' 				stack
	`Vector' 				keys1_by_2
	`Vector' 				keys2_by_1
	`Integer'				num_subgraphs
	`Variable'				subgraph_id		// (optional)

	// Computed by compute_cores()
	`Vector' 				cores
	`Vector' 				drop_order

	// Computed after prune_1core()
	`Integer'               N_drop
	`Variable'              mask            // mask (0|1) of obs that are dropped after prunning of degree-1 edges
	`Boolean'               prune           // Whether to recursively prune degree-1 edges
	`Vector'                drop2idx
	`Matrix'                drop2info
	`Variable'              sorted_w
	`Boolean'				has_weights
	`Variable'              sorted_true_weight


	// Methods
	`Void'					init()
	`Real'					init_zigzag()
	`Void'					compute_cores()
	`Void'					prune_1core()
	`Variables'				expand_1core()
	`Variables'				partial_out()
	`Variables'				__partial_out_map()
	`Variables'				__partial_out_laplacian()
}


`Void' BipartiteGraph::init(`FactorPointer' PF1,
                            `FactorPointer' PF2,
                            `Boolean' verbose)
{
	if (verbose) {
		printf("\n{txt}## Initializing bipartite graph\n\n")
		printf("    - FE #1: {res}%s{txt}\n", invtokens((*PF1).varlist))
		printf("    - FE #2: {res}%s{txt}\n", invtokens((*PF2).varlist))
	}
	this.verbose = verbose
	this.PF1 = PF1
	this.PF2 = PF2

	N = (*PF1).num_obs
	N1 = (*PF1).num_levels
	N2 = (*PF2).num_levels
	N12 = N1 + N2
	(*PF1).panelsetup() // Just in case
	(*PF2).panelsetup() // Just in case

	// F12 must be created from F1.levels and F2.levels (not from the original keys)
	// This is set automatically by join_factors() with the correct flag:
	//			F12 = join_factors(F1, F2, ., ., 1)
	// But you can also run (slower)
	//			F12 = _factor( (F1.levels, F2.levels) )
	//			asarray(F12.extra, "levels_as_keys", 1)
	if (verbose) printf("{txt}   - computing F12:  ")
	// join_factors(F1, (*PF2) [, count_levels, save_keys, levels_as_keys])
	F12 = join_factors((*PF1), (*PF2), ., ., 1)
	if (verbose) printf("{txt} edges found: {res}%-10.0gc{txt}\n", F12.num_levels)
	F12.panelsetup()
	
	if (verbose) printf("{txt}   - computing F12_1:")
	// _factor(data [, integers_only, verbose, method, sort_levels, count_levels, hash_ratio, save_keys])
	F12_1 = _factor(F12.keys[., 1], 1, 0, "", 1, 1, ., 0)
	if (verbose) printf("{txt} edges found: {res}%-10.0gc{txt}\n", F12_1.num_levels)
	F12_1.panelsetup()
	
	if (verbose) printf("{txt}   - computing F12_2:")
	F12_2 = _factor(F12.keys[., 2], 1, 0, "", 1, 1, ., 0)
	if (verbose) printf("{txt} edges found: {res}%-10.0gc{txt}\n", F12_2.num_levels)
	F12_2.panelsetup()
}


// --------------------------------------------------------------------------
// init_zigzag()
// --------------------------------------------------------------------------
//   Construct -queue- and -stack- vectors that allow zigzag iteration
//   
//   queue: firm and CEOs that will be processed, in the required order
//   	   note: negative values indicate CEOs
//   
//   stack: for each firm/CEO, the list of nodes it connects to
//          note: stacks are zero-separated
//   
// --------------------------------------------------------------------------
//   As a byproduct, also computes the number of disjoint subgraphs.
//   See the algorithm from on Abowd, Creecy and Kramarz (WP 2002) p4. Sketch:
//
//   		g = 0
//   		While there are firms w/out a group:
//   		     g++
//   		     Assign the first firm w/out a group to group g
//   		     Repeat until no further changes:
//   		         Add all persons employed by a firm in g to g
//   		         Add all firms that employ persons in g to g
//   		 return(g)
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
`Real' BipartiteGraph::init_zigzag(| `Boolean' save_subgraphs)
{
	`Vector'				counter1
	`Vector'				counter2
	`Vector' 				done1
	`Vector' 				done2

	`Integer'				i_stack		// use to process the queue
	`Integer'				last_i		// use to fill out the queue
	`Integer'				start_j		// use to search for firms to start graph enumeration
	`Integer'				i_queue
	`Integer'				id 			// firm number if id>0; error if id=0; ceo number if id<0
	`Integer'				j 			// firm # (or viceversa)
	`Integer'				k 			// ceo # (or viceversa)
	`Integer'				c 			// temporary counter
	`Integer'				i 			// temporary iterator
	
	`Matrix'				matches // list of CEOs that matched with firm j (or viceversa)

	if (verbose) printf("\n{txt}## Initializing zigzag iterator for bipartite graph\n\n")
	assert(F12_1.panel_is_setup)
	assert(F12_2.panel_is_setup)

	// If subgraph_id (mobility groups) is anything BUT zero, we will save them
	if (args()==0 | save_subgraphs==.) save_subgraphs = 0
	if (save_subgraphs) {
		subgraph_id = J(N2, 1, .)
	}

	queue = J(N12, 1, 0)
	stack = J(F12.num_levels + N12, 1, .) // there are N12 zeros
	counter1 = J(N1, 1, 0)
	counter2 = J(N2, 1, 0)

	keys1_by_2 = F12_2.sort(F12.keys[., 1])
	keys2_by_1 = F12_1.sort(F12.keys[., 2])
	done1 = J(N1, 1, 0) // if a firm is already on the queue
	done2 = J(N2, 1, 0) // if a CEO is already on the queue

	// Use -j- for only for firms and -k- only for CEOs
	// Use -i_queue- to iterate over the queue and -i_stack- over the stack
	// Use -last_i- to fill out the queue (so its the last filled value)
	// Use -i- to iterate arbitrary vectors
	// Use -id- to indicate a possible j or k (negative for k)
	// Use -start_j- to remember where to start searching for new subgraphs

	i_stack = 0
	last_i = 0
	start_j = 1
	num_subgraphs = 0

	for (i_queue=1; i_queue<=N12; i_queue++) {
		id = queue[i_queue] // >0 if firm ; <0 if CEO; ==0 if nothing yet
		j = k = . // just to avoid bugs
		
		// Pick starting point (useful if the graph is disjoint!)
		if (id == 0) {
			assert(last_i + 1 == i_queue)
			for (j=start_j; j<=N1; j++) {
				if (!done1[j]) {
					queue[i_queue] = id = j
					start_j = j + 1
					++last_i
					break
				}
			}
			// printf("{txt} - starting subgraph with firm %g\n", j)
			++num_subgraphs
			assert(id != 0) // Sanity check
		}

		if (id > 0) {
			// It's a firm
			j = id
			done1[j] = 1
			matches = panelsubmatrix(keys2_by_1, j, F12_1.info)
			for (i=1; i<=rows(matches); i++) {
				k = matches[i]
				c = counter2[k]
				counter2[k] = c + 1
				if (!done2[k]) {
					if (!c) {
						queue[++last_i] = -k
					}
				 	stack[++i_stack] = k
				}
			}
			stack[++i_stack] = 0
		}
		else {
			// It's a CEO
			k = -id
			done2[k] = 1
			matches = panelsubmatrix(keys1_by_2, k, F12_2.info)
			for (i=1; i<=rows(matches); i++) {
				j = matches[i]
				c = counter1[j]
				counter1[j] = c + 1
				if (!done1[j]) {
					if (!c) {
						queue[++last_i] = j
					}
					stack[++i_stack] = j
				}
			}
			stack[++i_stack] = 0
			if (save_subgraphs) subgraph_id[k] = num_subgraphs
		}
	}

	// Sanity checks
	assert(counter1 == F12_1.counts)
	assert(counter2 == F12_2.counts)
	assert(!anyof(queue, 0)) // queue can't have zeros at the end
	assert(allof(done1, 1))
	assert(allof(done2, 1))
	assert(!missing(queue))
	assert(!missing(stack))

	if (save_subgraphs) subgraph_id = subgraph_id[(*PF2).levels]
	
	if (verbose) printf("{txt}   - disjoint subgraphs found: {res}%g{txt}\n", num_subgraphs)
	return(num_subgraphs)
}


// --------------------------------------------------------------------------
// compute_cores()
// --------------------------------------------------------------------------
// 	  Computes vertex core numbers, which allows k-core pruning
//	  Algorithm used is listed here: https://arxiv.org/abs/cs/0310049
// --------------------------------------------------------------------------
// Note:
//    maybe use the k-cores for something useful? eg:
//    we might want to weight the core numbers by the strength (# of obs together)
//    https://arxiv.org/pdf/1611.02756.pdf --> # of butterflies in bipartite graph
//    this paper also has useful data sources for benchmarks
//    # of primary and secondary vertices, edges
// --------------------------------------------------------------------------

`Void' BipartiteGraph::compute_cores()
{
	`Factor'				Fbin
	`Boolean'				is_firm
	`Integer'				M, ND, j, jj
	`Integer'				i_v, i_u, i_w
	`Integer'				pv, pu, pw
	`Integer'				v, u, w
	`Integer'				dv, du
	`Vector'				bin, deg, pos, invpos, vert, neighbors

	if (verbose) printf("\n{txt}## Computing vertex core numbers\n\n")

	// v, u, w are vertices; <0 for CEOs and >0 for firms
	// vert is sorted by degree; deg is unsorted
	// pos[i] goes from sorted to unsorted, so:
	// 		vert[i] === original_vert[ pos[i] ]
	// invpos goes from unsorted to sorted, so:
	//		vert[invpos[j]] === original_vert[j]

	// i_u represents the pos. of u in the sorted tables
	// pu represents the pos. of u in the unsorted/original tables

	assert(F12_1.panel_is_setup==1)
	assert(F12_2.panel_is_setup==1)
	assert(rows(queue)==N12)
	assert(rows(keys1_by_2)==F12.num_levels)
	assert(rows(keys2_by_1)==F12.num_levels)

	deg = F12_1.counts \ F12_2.counts
	ND = max(deg) // number of degrees

	Fbin = _factor(deg, 1, 0)
	Fbin.panelsetup()

	bin = J(ND, 1, 0)
	bin[Fbin.keys] = Fbin.counts
	bin = rows(bin) > 1 ? runningsum(1 \ bin[1..ND-1]) : 1
	
	pos = Fbin.p
	invpos = invorder(Fbin.p)

	vert = Fbin.sort(( (1::N1) \ (-1::-N2) ))

	for (i_v=1; i_v<=N12; i_v++) {
		v = vert[i_v]
		is_firm = (v > 0)

		neighbors = is_firm ? panelsubmatrix(keys2_by_1, v, F12_1.info) : panelsubmatrix(keys1_by_2, -v, F12_2.info)
		M = rows(neighbors)
		
		for (j=1; j<=M; j++) {	
			pv = pos[i_v]
			jj = neighbors[j]
			pu = is_firm ? N1 + jj : jj // is_firm is *not* for the neighbor
			dv = deg[pv]
			du = deg[pu]
		
			if (dv < du) {
				i_w = bin[du]
				w = vert[i_w]
				u = is_firm ? -jj : jj // is_firm is *not* for the neighbor
				if (u != w) {
					pw = pos[i_w]
					i_u = invpos[pu]
					pos[i_u] = pw
					pos[i_w] = pu
					vert[i_u] = w
					vert[i_w] = u
					invpos[pu] = i_w
					invpos[pw] = i_u
				}
				bin[du] = bin[du] + 1
				deg[pu] = deg[pu] - 1
			}
		} // end for neighbor u (u ~ v)
	} // end for each node v
	
	if (verbose) {
		//printf("{txt}      Table: core numbers and vertex count\n")
		Fbin = _factor(deg, 1, 0)
		//printf("\n")
		mm_matlist(Fbin.counts, "%-8.0gc", 2, strofreal(Fbin.keys), "Freq.", "Core #")
		printf("\n")
	}

	// ((F1.keys \ F2.keys), (F12_1.keys \ -F12_2.keys))[selectindex(deg:==1), .]		
	
	// Store the values in the class
	swap(drop_order, vert)
	swap(cores, deg)
}

// --------------------------------------------------------------------------
// prune_1core()
// --------------------------------------------------------------------------
//    Prune edges with degree-1
//    That is, recursively remove CEOs that only worked at one firm,
//    and firms that only had one CEO in the sample, until every agent
//    in the dataset has at least two matches
// --------------------------------------------------------------------------
`Void' BipartiteGraph::prune_1core(| `Variable' weight)
{
	`Integer'               N_drop2, i, j, i1, i2, j1, j2, K_drop2
	`Vector'                drop1, drop2
	`Vector'                tmp_mask
	`Vector'                proj1, proj2
	`Variable'              w, tmp_weight

	has_weights = (args()>0 & rows(weight) > 1)
	if (has_weights) sorted_true_weight = (*PF1).sort(weight)
	tmp_weight = has_weights ? weight : J(N, 1, 1)

	N_drop = sum(cores :== 1)
	if (!N_drop) {
	    if (verbose) printf("{txt}   - no 1-core vertices found\n")
	    prune = 0
	    return
	}
	if (verbose) printf("{txt}   - 1-core vertices found: {res}%g{txt}\n", N_drop)

	drop_order = drop_order[1..N_drop]
	drop1 = `selectindex'(cores[1..N1] :== 1)
	cores = .
	drop1 = (1::N1)[drop1]
	drop2 = -select(drop_order, drop_order:<0)
	
	K_drop2 = rows(drop2)
	N_drop2 = K_drop2 ? sum((*PF2).info[drop2, 2] :- (*PF2).info[drop2, 1] :+ 1) : 0

	tmp_mask = J(N1, 1, 0)
	if (rows(drop1)) tmp_mask[drop1] = J(rows(drop1), 1, 1)
	mask = tmp_mask[(*PF1).levels, 1]
	tmp_mask = J(N2, 1, 0)
	if (K_drop2) tmp_mask[drop2] = J(K_drop2, 1, 1)
	mask = mask :| tmp_mask[(*PF2).levels, 1]
	tmp_mask = .

	drop2idx = J(N_drop2, 1, .)
	drop2info = J(N2, 2, .)

	j1 = 1
	for (i=1; i<=K_drop2; i++) {
	    j = drop2[i]
	    i1 = (*PF2).info[j, 1]
	    i2 = (*PF2).info[j, 2]

	    j2 = j1 + i2 - i1
	    drop2idx[j1::j2] = i1::i2
	    drop2info[j, .] = (j1, j2)
	    j1 = j2 + 1
	}

	if (!(*PF2).is_sorted) {
	    assert(((*PF2).p != J(0, 1, .)))
	    drop2idx = (*PF2).p[drop2idx, .]
	}

	if (!(*PF1).is_sorted) {
	    assert(((*PF1).inv_p != J(0, 1, .)))
	    drop2idx = invorder((*PF1).p)[drop2idx, .]
	}

	// To undo pruning, I need (*PF1).info[drop1, .] & drop2info & drop2idx

	// Set weights of pruned obs. to zero
	tmp_weight[`selectindex'(mask)] = J(sum(mask), 1, 0)

	// Update sorted weights for g=1,2
	w = (*PF1).sort(tmp_weight)
	asarray((*PF1).extra, "has_weights", 1)
	asarray((*PF1).extra, "weights", w)
	asarray((*PF1).extra, "weighted_counts", `panelsum'(w, (*PF1).info))
	w = .

	w = (*PF2).sort(tmp_weight)
	tmp_weight = . // cleanup
	asarray((*PF2).extra, "has_weights", 1)
	asarray((*PF2).extra, "weights", w)
	asarray((*PF2).extra, "weighted_counts", `panelsum'(w, (*PF2).info))
	w = .
	
	// Select obs where both FEs are degree-1 (and thus omitted)
	sorted_w = J(N, 1, 1)
	
	proj1 = panelmean((*PF1).sort(sorted_w), *PF1)[(*PF1).levels, .]
	proj2 = panelmean((*PF2).sort(sorted_w), *PF2)[(*PF2).levels, .]
	sorted_w = ((sorted_w - proj1) :!= 1) :| ((sorted_w - proj2) :!= 1)
	proj1 = proj2 = .
	sorted_w = (*PF1).sort(sorted_w)

	prune = 1
}

// --------------------------------------------------------------------------
// prune_1core()
// --------------------------------------------------------------------------
//    Prune edges with degree-1
//    That is, recursively remove CEOs that only worked at one firm,
//    and firms that only had one CEO in the sample, until every agent
//    in the dataset has at least two matches
// --------------------------------------------------------------------------
`Variables' BipartiteGraph::expand_1core(`Variables' y)
{
    `Boolean'               zero_weights
    `Variable'              sorted_y
    `Integer'               i, j, j1, j2, i2, k1, k2, nk
    `Matrix'                tmp_y
    `Vector'                tmp_w, tmp_idx, new_w
    `RowVector'             tmp_mean

    if (prune==0) return(y)
    if (verbose) printf("\n{txt}## Expanding 2-core into original dataset\n\n")
    assert(N_drop == rows(drop_order))

    sorted_y = (*PF1).sort(y)

    i2 = 0
    for (i=N_drop; i>=1; i--) {
        j = drop_order[i]
        if (j > 0) {
            j1 = (*PF1).info[j, 1]
            j2 = (*PF1).info[j, 2]

            tmp_y = sorted_y[| j1 , 1 \ j2 , . |] // panelsubmatrix(sorted_y, j, (*PF1).info)
            tmp_w = sorted_w[|j1, 1 \ j2, .|] // panelsubmatrix(sorted_w, j, (*PF1).info)
            new_w = has_weights ? sorted_true_weight[|j1, 1 \ j2, .|] : J(j2-j1+1, 1, 1)
            zero_weights = !sum(tmp_w)
            if (!zero_weights) {
                tmp_mean = mean(tmp_y, tmp_w)
                assert(!missing(tmp_mean)) // bugbug remove later
                sorted_y[| j1 , 1 \ j2 , . |] = tmp_y :- tmp_mean
            }
            sorted_w[| j1 , 1 \ j2 , 1 |] = new_w
        }
        else {
            ++i2
            j1 = drop2info[-j, 1]
            j2 = drop2info[-j, 2]
            tmp_idx = drop2idx[| j1 , 1 \ j2 , 1 |]
            tmp_y = sorted_y[tmp_idx, .]
            tmp_w = sorted_w[tmp_idx]
            zero_weights = !sum(tmp_w)
            if (zero_weights) {
                tmp_w = has_weights ? sorted_true_weight[tmp_idx] : J(j2-j1+1, 1, 1)
            }
            tmp_mean = mean(tmp_y, tmp_w)
            assert(!missing(tmp_mean)) // bugbug remove later
            nk = rows(tmp_idx)
            for (k1=1; k1<=nk; k1++) {
               k2 = tmp_idx[k1]
               sorted_y[k2, .] = sorted_y[k2, .] - tmp_mean
               sorted_w[k2] = has_weights ? sorted_true_weight[k2] : 1
            }
        }
    }

    if (verbose) printf("{txt}   - number of coefficients solved triangularly: {res}%s{txt}\n", strofreal(rows(drop_order)))
    return((*PF1).invsort(sorted_y))
}


`Variables' BipartiteGraph::partial_out(`Variables' y)
{

}


`Variables' BipartiteGraph::__partial_out_map(`Variables' y)
{

}


`Variables' BipartiteGraph::__partial_out_laplacian(`Variables' y)
{

}

end

// --------------------------------------------------------------------------
// FixedEffects main class
// --------------------------------------------------------------------------

mata:

class FixedEffects
{
	// Factors
	`Integer'               G                   // Number of sets of FEs
	`Integer'               N                   // number of obs
	`Integer'               M                   // Sum of all possible FE coefs
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

	// Constant-related (also see -has_intercept-)
	`Boolean'               report_constant
	`Boolean'               compute_constant

	// Optimization options
	`Real'                  tolerance
	`Real'                  extra_tolerance		// Try to achieve this tol if it only takes a few more iters: ceil(10%)
	`Integer'               maxiter
	`String'                transform           // Kaczmarz Cimmino Symmetric_kaczmarz (k c s)
	`String'                acceleration        // Acceleration method. None/No/Empty is none\
	`Integer'               accel_start         // Iteration where we start to accelerate // set it at 6? 2?3?
	`string'                slope_method
	`Boolean'               prune               // Whether to recursively prune degree-1 edges
	`Boolean'               abort               // Raise error if convergence failed?
	`Integer'               accel_freq          // Specific to Aitken's acceleration
	`Boolean'               storing_alphas      // 1 if we should compute the alphas/fes
	`Real'                  conlim              // specific to LSMR
	`Real'                  btol                // specific to LSMR
	`Boolean'				always_run_lsmr_preconditioner
	`Integer'				min_ok

	// Optimization objects
	`BipartiteGraph'        bg                  // Used when pruning 1-core vertices
	`Vector'                pruned_weight       // temp. weight for the factors that were pruned
	`Integer'               prune_g1            // Factor 1/2 in the bipartite subgraph that gets pruned
	`Integer'               prune_g2            // Factor 2/2 in the bipartite subgraph that gets pruned
	`Integer'               num_pruned          // Number of vertices (levels) that were pruned

	// Misc
	`Integer'               verbose
	`Boolean'               timeit
	`Boolean'               compact
	`Integer'               poolsize
	`Boolean'               store_sample
	`Real'                  finite_condition
	`Real'                  compute_rre         // Relative residual error: || e_k - e || / || e ||
	`Real'                  rre_depvar_norm
	`Vector'                rre_varname
	`Vector'                rre_true_residual
	`String'                panelvar
	`String'                timevar

	`RowVector'             not_basevar         // Boolean vector indicating whether each regressor is or not a basevar
	`String'                fullindepvars       // indepvars including basevars

	// Weight-specific
	`Boolean'               has_weights
	`Variable'              weight              // unsorted weight
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
	`String'                varlist             // y x1 x2 x3
	`String'                depvar              // y
	`String'                indepvars           // x1 x2 x3
	`String'                tousevar
	
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
	`Variable'              residuals_vector
	`RowVector'             kept // 1 if the regressors are not deemed as omitted (by partial_out+cholsolve+invsym)
	`String'                diopts

	// Output
	`String'                cmdline
	`String'                subcmd
	`String'                title
	`Boolean'               converged
	`Integer'               iteration_count // e(ic)
	`Varlist'               extended_absvars
	`String'                notes
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
	`Real'                  accuracy
	`RowVector'             means
	`RowVector'				all_stdevs

	// Methods
	`Void'                  new()
	`Void'                  destroy()
	`Void'                  load_weights() // calls update_sorted_weights, etc.
	`Void'                  update_sorted_weights()
	`Void'                  update_cvar_objects()
	`Matrix'                partial_out()
	`Matrix'                partial_out_pool()
	`Void'                  _partial_out()
	`Variables'             project_one_fe()
	`Void'                  prune_1core()
	`Void'                  _expand_1core()
	`Void'                  estimate_dof()
	`Void'                  estimate_cond()
	`Void'                  save_touse()
	`Void'                  store_alphas()
	`Void'                  save_variable()
	`Void'                  post_footnote()
	`Void'                  post()
	`FixedEffects'          reload() // create new instance of object

	// LSMR-Specific Methods
	`Real'                  lsmr_norm()
	`Vector'                lsmr_A_mult()
	`Vector'                lsmr_At_mult()
}    


// Set default value of properties
`Void' FixedEffects::new()
{
	num_singletons = .
	sample = J(0, 1, .)
	weight = 1 // set to 1 so cross(x, S.weight, y)==cross(x, y)

	verbose = 0
	timeit = 0
	compact = 0
	poolsize = .
	finite_condition = .
	residuals = ""
	residuals_vector = .
	panelvar = timevar = ""
	iteration_count = 0
	accuracy = -1 // Epsilon at the time of convergence

	// Optimization defaults
	slope_method = "invsym"
	maxiter = 1e4
	tolerance = 1e-8
	transform = "symmetric_kaczmarz"
	acceleration = "conjugate_gradient"
	accel_start = 6
	conlim = 1e+8 // lsmr only
	btol = 1e-8 // lsmr only (note: atol is just tolerance)
	always_run_lsmr_preconditioner = 0
	min_ok = 1

	prune = 0
	converged = 0
	abort = 1
	storing_alphas = 0
	report_constant = compute_constant = 1

	// Specific to Aitken:
	accel_freq = 3

	not_basevar = J(1, 0, .)

	means = all_stdevs = J(1, 0, .) // necessary with pool() because we append to it
	kept = J(1, 0, .) // necessary with pool() because we append to it
}


`Void' FixedEffects::destroy()
{
	// stata(sprintf("cap drop %s", tousevar))
}


// This adds/removes weights or changes their type
`Void' FixedEffects::load_weights(`String' weighttype, `String' weightvar, `Variable' weight, `Boolean' verbose)
{
	`Integer'				g
	`FactorPointer'         pf
	`Matrix'                precond // used for lsmr
	`Varname'               cvars_g
   
	this.has_weights = (weighttype != "" & weightvar != "")
	if (this.verbose > 0 & verbose > 0 & this.has_weights) printf("{txt}## Loading weights [%s=%s]\n", weighttype, weightvar)

	// Update main properties
	this.weight_var = weightvar
	this.weight_type = weighttype

	// Update booleans
	for (g=1; g<=this.G; g++) {
		asarray(this.factors[g].extra, "has_weights", this.has_weights)
	}

	// Optionally load weight from dataset
	if (this.has_weights & weight==J(0,1,.)) {
		weight = st_data(this.sample, this.weight_var)
	}

	// Update weight vectors
	if (this.has_weights) {
		if (this.verbose > 0 & verbose > 0) printf("{txt}## Sorting weights for each absvar\n")
		this.update_sorted_weights(weight)
	}
	else {
		// If no weights, clear this up
		this.weight = 1 // same as defined by new()
		for (g=1; g<=this.G; g++) {
			asarray(this.factors[g].extra, "weights", .)
			asarray(this.factors[g].extra, "weighted_counts", .)
		}
	}
	
	// Update cvar objects (do AFTER updating weights!)
	// (this is meaningless with iweights)
	if (weighttype != "iweight") this.update_cvar_objects()

	// Preconditioners for LSMR
	if (acceleration=="lsmr" | always_run_lsmr_preconditioner) {

		// Compute M
		M = 0
		for (g=1; g<=G; g++) {
			M = M + factors[g].num_levels * (intercepts[g] + num_slopes[g])
		}

		// Preconditioner
		for (g=1; g<=G; g++) {
			pf = &(factors[g])
			if (intercepts[g]) {
				precond = has_weights ? asarray((*pf).extra, "weighted_counts") : (*pf).counts
				asarray((*pf).extra, "precond_intercept", sqrt(1 :/ precond))
			}

			if (num_slopes[g]) {
				cvars_g = tokens(this.cvars[g])
				precond = st_data(this.sample, cvars_g)
				precond = reghdfe_panel_precondition(precond, (*pf))
				asarray((*pf).extra, "precond_slopes", precond)
			}

			precond = .
		}
	}

}


// This just updates the weight but doesn't change the type or variable of the weight
`Void' FixedEffects::update_sorted_weights(`Variable' weight)
{
	`Integer'               g
	`Real'                  min_w
	`Variable'              w
	`FactorPointer'         pf

	assert_msg(!hasmissing(weight), "weights can't be missing")
	this.weight = weight
	assert(rows(weight)==rows(sample))
	if (verbose > 0) printf("{txt}   - loading %s weight from variable %s\n", weight_type, weight_var)
	for (g=1; g<=G; g++) {
		if (verbose > 0) printf("{txt}   - sorting weight for factor {res}%s{txt}\n", absvars[g])
		pf = &(factors[g])
		w = (*pf).sort(weight)

		// Rescale weights so there are no weights below 1
		if (weight_type != "fweight") {
			min_w = colmin(w)
			if (min_w < 1e-6) min_w = 1e-6 // Prevent bugs if a weight is very close to zero
			//assert_msg(min_w > 0, "weights must be positive")
			//if (min_w <= 0) printf("{err} not all weights are positive\n")
			if (0 < min_w & min_w < 1) {
				w = w :/ min_w
			}
		}

		asarray((*pf).extra, "weights", w)
		asarray((*pf).extra, "weighted_counts", `panelsum'(w, (*pf).info))
	}
}


`Void' FixedEffects::update_cvar_objects()
{
	`Integer'               g
	`FactorPointer'         pf

	for (g=1; g<=G; g++) {
		pf = &(factors[g])
		// Update mean(z; w) and inv(z'z; w) where z is a slope variable and w is the weight
		if (num_slopes[g]) {
			if (verbose > 0) printf("{txt}   - precomputing cvar objects for factor {res}%s{txt}\n", absvars[g])
			if (intercepts[g]) {
			    asarray((*pf).extra, "xmeans",
			            panelmean(asarray((*pf).extra, "x"), *pf))
			}
			asarray((*pf).extra, "inv_xx", precompute_inv_xx(*pf, intercepts[g]))
		}
	}
}


`Variables' FixedEffects::partial_out(`Anything' data,
									| `Boolean' save_tss,
									  `Boolean' standardize_data,
									  `Boolean' first_is_depvar)
{
	// -data- is either a varlist or a matrix
	`Variables'             y
	`Varlist'               vars
	`Integer'               i
	`Integer'               k

	if (args()<2 | save_tss==.) save_tss = 0
	if (args()<3 | standardize_data==.) standardize_data = 1
	if (args()<4 | first_is_depvar==.) first_is_depvar = 1

	if (eltype(data) == "string") {
		vars = tokens(invtokens(data)) // tweak to allow string scalars and string vectors
		k = cols(vars)

		if (poolsize < k) {
			if (verbose > 0) printf("\n{txt}## Loading and partialling out %g variables in blocks of %g\n\n", k, poolsize)
			if (timeit) timer_on(50)
			partial_out_pool(vars, save_tss, standardize_data, first_is_depvar, poolsize, y=.)
			if (timeit) timer_off(50)
		}
		else {
			if (verbose > 0) printf("\n{txt}## Partialling out %g variables: {res}%s{txt}\n\n", cols(vars), invtokens(vars))
			if (verbose > 0) printf("{txt}   - Loading variables into Mata\n")
			if (timeit) timer_on(50)
			_st_data_wrapper(sample, invtokens(vars), y=., verbose)
			if (timeit) timer_off(50)

			// Workaround to odd Stata quirk
			if (timeit) timer_on(51)
			if (cols(y) > cols(vars)) {
				printf("{err}(some empty columns were added due to a bug/quirk in {bf:st_data()}; %g cols created instead of %g for {it:%s}; running slower workaround)\n", cols(y), cols(vars), invtokens(vars))
				partial_out_pool(vars, save_tss, standardize_data, first_is_depvar, 1, y=.)
			}
			else {
				_partial_out(y, save_tss, standardize_data, first_is_depvar)
			}
			if (timeit) timer_off(51)
			
		}
	}
	else {
		if (verbose > 0) printf("\n{txt}## Partialling out %g variables\n\n", cols(data))
		if (timeit) timer_on(54)
		_partial_out(y=data, save_tss, standardize_data, first_is_depvar)
		if (timeit) timer_off(54)
	}

	if (verbose==0) printf(`"{txt}({browse "http://scorreia.com/research/hdfe.pdf":MWFE estimator} converged in %s iteration%s)\n"', strofreal(iteration_count), iteration_count > 1 ? "s" : "s")
	return(y)
}



`Variables' FixedEffects::partial_out_pool(`Anything' vars,
										   `Boolean' save_tss,
										   `Boolean' standardize_data,
										   `Boolean' first_is_depvar,
										   `Integer' step,
										   `Variables' y)
{
	`Variables'             part_y
	`Integer'               i, j, ii
	`Integer'               k
	`StringRowVector'       keepvars

	k = cols(vars)
	assert(step > 0)
	assert(step < k)
	y = J(rows(sample), 0, .)

	for (i=1; i<=k; i=i+step) {
		
		j = i + step - 1
		if (j>k) j = k

		// Load data
		_st_data_wrapper(sample, vars[i..j], part_y=., verbose)

		if (cols(part_y) > j - i + 1) {
			printf("{err}(some empty columns were added due to a bug/quirk in {bf:st_data()}; running slower workaround)\n")
			if (timeit) timer_on(51)
			part_y = J(rows(sample), 0, .)
			for (ii=i; ii<=j; ii++) {
				part_y = part_y, st_data(sample, vars[ii])
			}
			if (timeit) timer_off(51)
		}

		// Drop loaded vars as quickly as possible
		if (compact) {
			// st_dropvar(vars[i..j]) // bugbug what if repeated??
			keepvars = base_clustervars , timevar, panelvar, (j == k ? "" : vars[j+1..k])
			keepvars = tokens(invtokens(keepvars))
			if (cols(keepvars)) {
				stata(sprintf("fvrevar %s, list", invtokens(keepvars)))
				stata(sprintf("keep %s", st_global("r(varlist)")))
			}
			else {
				stata("clear")
			}
		}

		_partial_out(part_y, save_tss, standardize_data, first_is_depvar)
		y = y, part_y
		part_y = .
	}
}


`Void' FixedEffects::store_alphas(`Anything' d_varname)
{
	`Integer'               g, i, j
	`StringRowVector'       varlabel
	`Variable'              d
	`RowVector'             tmp_stdev

	if (verbose > 0) printf("\n{txt}## Storing estimated fixed effects\n\n")

	// -d- can be either the data or the variable name

	// Load -d- variable
	if (eltype(d_varname) == "string") {
		if (verbose > 0) printf("{txt}   - Loading d = e(depvar) - xb - e(resid)\n")
		d = st_data(sample, d_varname)
	}
	else {
		d = d_varname
	}
	assert(!missing(d))

	// Create empty alphas
	if (verbose > 0) printf("{txt}   - Initializing alphas\n")
	for (g=j=1; g<=G; g++) {
		if (!save_fe[g]) continue
		asarray(factors[g].extra, "alphas", J(factors[g].num_levels, intercepts[g] + num_slopes[g], 0))
		asarray(factors[g].extra, "tmp_alphas", J(factors[g].num_levels, intercepts[g] + num_slopes[g], 0))
	}

	// Fill out alphas
	if (verbose > 0) printf("{txt}   - Computing alphas\n")
	storing_alphas = 1
	converged = 0
	d = accelerate_sd(this, d, &transform_kaczmarz())
	storing_alphas = 0

	if (verbose > 0) printf("{txt}   - SSR of d wrt FEs: %g\n", quadcross(d,d))

	// Store alphas in dataset
	if (verbose > 0) printf("{txt}   - Creating varlabels\n")
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
			if (verbose > 0) printf("{txt}   - Recovering unstandardized variables\n")
			tmp_stdev = asarray(factors[g].extra, "x_stdevs")
			if (intercepts[g]) tmp_stdev = 1, tmp_stdev

			// We need to *divide* the coefs by the stdev, not multiply!
			asarray(factors[g].extra, "alphas",
				asarray(factors[g].extra, "alphas") :/ tmp_stdev
			)
		}

		if (verbose > 0) printf("{txt}   - Storing alphas in dataset\n")
		save_variable(targets[g], asarray(factors[g].extra, "alphas")[factors[g].levels, .], varlabel)
		asarray(factors[g].extra, "alphas", .)
		asarray(factors[g].extra, "tmp_alphas", .)
	}
}


`Void' FixedEffects::_partial_out(`Variables' y,
								| `Boolean' save_tss,
								  `Boolean' standardize_data,
								  `Boolean' first_is_depvar,
								  `Boolean' flush)
{
	`RowVector'             stdevs, needs_zeroing, kept2
	`FunctionP'             funct_transform, func_accel // transform
	`Real'                  y_mean, collinear_tol
	`Vector'                lhs
	`Vector'                alphas
	`StringRowVector'       vars
	`Integer'               i

	if (args()<2 | save_tss==.) save_tss = 0
	if (args()<3 | standardize_data==.) standardize_data = 1
	if (args()<4 | first_is_depvar==.) first_is_depvar = 1
	if (args()<5 | flush==.) flush = 0 // whether or not to flush the values of means & kept

	assert(anyof((0, 1, 2), standardize_data)) // 0=Don't standardize; 1=Std. and REVERT after partial; 2=Std., partial, and KEEP STANDARDIZED

	if (flush) {
		iteration_count = 0
		accuracy = -1
		means = stdevs = J(1, 0, .)
		kept = J(1, 0, .)
	}

	// Shortcut for trivial case (1 FE)
	if (G==1) acceleration = "none"

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

	// Compute TSS of depvar
	if (timeit) timer_on(60)
	if (save_tss & tss==.) {
		lhs = y[., 1]
		if (has_intercept) {
			y_mean = mean(lhs, weight)
			tss = crossdev(lhs, y_mean, weight, lhs, y_mean) // Sum of w[i] * (y[i]-y_mean) ^ 2
		}
		else {
			tss = cross(lhs, weight, lhs) // Sum of w[i] * y[i] ^ 2
		}
		lhs = .
		if (weight_type=="aweight" | weight_type=="pweight") tss = tss * rows(y) / sum(weight)
	}
	if (timeit) timer_off(60)


	// Compute 2-norm of each var, to see if we need to drop as regressors
	kept2 = diagonal(cross(y, y))'

	// Compute and save means of each var
	means = means , ( compute_constant ? mean(y, weight) : J(1, cols(y), 1) )

	// Intercept LSMR case
	if (acceleration=="lsmr") {
		// RRE benchmarking
		if (compute_rre) rre_depvar_norm = norm(y[., 1])
		if (cols(y)==1) {
			y = lsmr(this, y, alphas=.)
			alphas = . // or return them!
		}
		else {
			for (i=1; i<=cols(y); i++) {
				y[., i] = lsmr(this, y[., i], alphas=.)
			}
			alphas = .
		}
	}
	else {

		// Standardize variables
		if (timeit) timer_on(61)
		if (standardize_data) {
			if (verbose > 0) printf("{txt}   - Standardizing variables\n")
			stdevs = reghdfe_standardize(y)
			all_stdevs = all_stdevs, stdevs
			kept2 = kept2 :/ stdevs :^ 2
		}
		if (timeit) timer_off(61)

		// RRE benchmarking
		if (compute_rre) {
			rre_true_residual = rre_true_residual / (standardize_data ? stdevs[1] : 1)
			rre_depvar_norm = norm(y[., 1])
		}

		// Solve
		if (verbose>0) printf("{txt}   - Running solver (acceleration={res}%s{txt}, transform={res}%s{txt} tol={res}%-1.0e{txt})\n", acceleration, transform, tolerance)
		if (verbose==1) printf("{txt}   - Iterating:")
		if (verbose>1) printf("{txt}      ")
		converged = 0 // converged will get updated by check_convergence()

		if (timeit) timer_on(62)
		if (G==1 & factors[1].method=="none" & num_slopes[1]==0 & !(storing_alphas & save_fe[1])) {
			// Speedup for constant-only case (no fixed effects)
			assert(factors[1].num_levels == 1)
			iteration_count = 1
			accuracy = 0
			if (standardize_data == 1) {
				y = stdevs :* y :- stdevs :* mean(y, has_weights ? asarray(factors[1].extra, "weights") : 1) // Undoing standardization
			}
			else {
				y = y :- mean(y, has_weights ? asarray(factors[1].extra, "weights") : 1)
			}
		}
		else {
			if (standardize_data == 1) {
				y = (*func_accel)(this, y, funct_transform) :* stdevs // Undoing standardization
			}
			else {
				y = (*func_accel)(this, y, funct_transform) // 'this' is like python's self
			}
		}
		if (timeit) timer_off(62)
		
		if (prune) {
			assert_msg(G==2, "prune option requires only two FEs")
			if (timeit) timer_on(63)
			_expand_1core(y)
			if (timeit) timer_off(63)
		}
	}

	assert_msg(!hasmissing(y), "error partialling out; missing values found")

	// Standardizing makes it hard to detect values that are perfectly collinear with the absvars
	// in which case they should be 0.00 but they end up as e.g. 1e-16
	// EG: reghdfe price ibn.foreign , absorb(foreign)

	// This will edit to zero entire columns where *ALL* values are very close to zero
	if (timeit) timer_on(64)
	vars = cols(varlist) > 1 ? varlist : tokens(varlist)
	if (cols(vars)!=cols(y)) vars ="variable #" :+ strofreal(1..cols(y))
	collinear_tol = min(( 1e-6 , tolerance / 10))

	kept2 = (diagonal(cross(y, y))' :/ kept2) :> (collinear_tol)
	if (first_is_depvar & kept2[1]==0) {
		kept2[1] = 1
		if (verbose > -1) printf("{txt}warning: %s might be perfectly explained by fixed effects (tol =%3.1e)\n", vars[1], collinear_tol)
	}
	needs_zeroing = `selectindex'(!kept2)
	if (cols(needs_zeroing)) {
		y[., needs_zeroing] = J(rows(y), cols(needs_zeroing), 0)
		for (i=1; i<=cols(vars); i++) {
			if (!kept2[i] & verbose>-1 & (i > 1 | !first_is_depvar)) {
				printf("{txt}note: {res}%s{txt} is probably collinear with the fixed effects (all partialled-out values are close to zero; tol =%3.1e)\n", vars[i], collinear_tol)
			}
		}
	}

	kept = kept, kept2
	if (timeit) timer_off(64)
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
			if (cols(y)==1 & f.num_levels > 1) {
				return(panelmean(f.sort(y), f)[f.levels])
			}
			else {
				return(panelmean(f.sort(y), f)[f.levels, .])
			}
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
	`Integer'               i_cluster, i_intercept, j_intercept, i_start
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
	`Boolean'				save_subgraph
	`String'				grouplabel
	
	if (verbose > 0) printf("\n{txt}## Estimating degrees-of-freedom absorbed by the fixed effects\n\n")

	// Count all FE intercepts and slopes
	SubGs = intercepts + num_slopes
	G_extended = sum(SubGs)
	num_intercepts = sum(intercepts)
	offsets = runningsum(SubGs) - SubGs :+ 1 // start of each FE within the extended list
	idx = `selectindex'(intercepts) // Select all FEs with intercepts
	if (verbose > 0) printf("{txt}   - there are %f fixed intercepts and slopes in the %f absvars\n", G_extended, G)

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
						clustervar = subinstr(clustervars[i_cluster], "#", " ", .)
						F = factor(clustervar, sample, ., "", 0, 0, ., 0)
						cluster_data = F.levels
						F = Factor() // clear
					}
					else {
						cluster_data = __fload_data(clustervar, sample, 0)
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
	// Note that this excludes the ones nested within clusters
	// If we DO have FEs nested within clusters, we should also include the first intercept
	if ((length(idx) > 1) | (length(idx) >= 1 & df_a_nested > 0)) {
		if (verbose > 0) printf("{txt}   - there is at least one redundant coef. for every set of FE intercepts after the first one\n")
		i_start = df_a_nested ? 1 : 2
		doflist_M[offsets[idx[i_start..length(idx)]]] = J(1, length(idx)-i_start+1, 1) // Set DoF loss of all intercepts but the first one to 1
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
				++pair_count
				save_subgraph = (pair_count == 1) & (groupvar != "")
				m = BG.init_zigzag(save_subgraph)
				if (verbose > 0) printf("{txt}   - mobility groups between FE intercepts #%f and #%f: {res}%f\n", g, h, m)
				if (save_subgraph) {
					if (verbose > 2) printf("{txt}   - Saving identifier for the first mobility group: {res}%s\n", groupvar)
					st_store(sample, st_addvar("long", groupvar), BG.subgraph_id)
					grouplabel = sprintf("Mobility group between %s and %s", invtokens(factors[g].varlist, "#"), invtokens(factors[h].varlist, "#"))
					st_varlabel(groupvar, grouplabel)
				}
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
				if (has_int  & verbose) printf("{txt}   - the slopes in the FE #%f are constant for {res}%f{txt} levels, which don't reduce DoF\n", g, sum(results))
				if (!has_int & verbose) printf("{txt}   - the slopes in the FE #%f are zero for {res}%f{txt} levels, which don't reduce DoF\n", g, sum(results))
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
	// Note that we can't prune degree-2 nodes, or the graph stops being bipartite
	`Integer'               i, j, g
	`Vector'                subgraph_id
	
	`Vector'                idx
	`RowVector'             i_prune

	// For now; too costly to use prune for G=3 and higher
	// (unless there are *a lot* of degree-1 vertices)
	if (G!=2) return //assert_msg(G==2, "G==2") // bugbug remove?

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
	(void) bg.init_zigzag(1) // 1 => save subgraphs into bg.subgraph_id
	bg.compute_cores()
	bg.prune_1core(weight)
	num_pruned = bg.N_drop
}

// bugbug fix or remove this fn altogether
`Void' FixedEffects::_expand_1core(`Variables' y)
{
	y = bg.expand_1core(y)
}


`Void' FixedEffects::save_touse(| `Varname' touse, `Boolean' replace)
{
	`Integer'               idx
	`Vector'                mask

	// Set default arguments
	if (args()<1 | touse=="") {
		assert(tousevar != "")
		touse = tousevar
	}
	// Note that args()==0 implies replace==1 (else how would you find the name)
	if (args()==0) replace = 1
	if (args()==1 | replace==.) replace = 0

	if (verbose > 0) printf("\n{txt}## Saving e(sample)\n")

	// Compute dummy vector
	mask = create_mask(st_nobs(), 0, sample, 1)

	// Save vector as variable
	if (replace) {
		st_store(., touse, mask)
	}
	else {
		idx = st_addvar("byte", touse, 1)
		st_store(., idx, mask)
	}
}


`Void' FixedEffects::save_variable(`Varname' varname,
								   `Variables' data,
								 | `Varlist' varlabel)
{
	`RowVector'               idx
	`Integer'               i
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
	st_global("e(dofmethod)", invtokens(dofadjustments))

	if (absvars == "") {
		absvars = extended_absvars = "_cons"
	}

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
	st_numscalar("e(N_full)", st_numscalar("e(N)") + num_singletons)

	// Prune specific
	if (prune==1) {
		st_numscalar("e(pruned)", 1)
		st_numscalar("e(num_pruned)", num_pruned)
	}

	if (!missing(finite_condition)) st_numscalar("e(finite_condition)", finite_condition)
}


`Void' FixedEffects::post()
{
	`String'        text
	`Integer'       i

	post_footnote()

	// ---- constants -------------------------------------------------------

	st_global("e(predict)", "reghdfe5_p")
	st_global("e(estat_cmd)", "reghdfe5_estat")
	st_global("e(footnote)", "reghdfe5_footnote")
	//st_global("e(marginsok)", "")
	st_global("e(marginsnotok)", "Residuals SCore")
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

	st_numscalar("e(report_constant)", compute_constant & report_constant)


	// ---- .options properties ---------------------------------------------

	st_global("e(depvar)", depvar)
	st_global("e(indepvars)", invtokens(indepvars))

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


// --------------------------------------------------------------------------
// Recreate HDFE object
// --------------------------------------------------------------------------
`FixedEffects' FixedEffects::reload(`Boolean' copy)
{
	`FixedEffects' ans
	assert(copy==0 | copy==1)
	
	// Trim down current object as much as possible
	// this. is optional but useful for clarity
	if (copy==0) {
		this.factors = Factor()
		this.sample = .
		this.bg = BipartiteGraph()
		this.pruned_weight = .
		this.rre_varname = .
		this.rre_true_residual = .
	}

	// Initialize new object
	ans = fixed_effects(this.absorb, this.tousevar, this.weight_type, this.weight_var, this.drop_singletons, this.verbose)

	// Fill out new object with values of current one
	ans.depvar = this.depvar
	ans.indepvars = this.indepvars
	ans.varlist = this.varlist
	ans.model = this.model
	ans.vcetype = this.vcetype
	ans.num_clusters = this.num_clusters
	ans.clustervars = this.clustervars
	ans.base_clustervars = this.base_clustervars
	ans.vceextra = this.vceextra
	ans.summarize_stats = this.summarize_stats
	ans.summarize_quietly = this.summarize_quietly
	ans.notes = this.notes
	ans.store_sample = this.store_sample
	ans.timeit = this.timeit
	ans.compact = this.compact
	ans.poolsize = this.poolsize
	ans.diopts = this.diopts

	ans.fullindepvars = this.fullindepvars
	ans.not_basevar = this.not_basevar

	ans.compute_constant = this.compute_constant
	ans.report_constant = this.report_constant
	ans.tolerance = this.tolerance
	ans.save_any_fe = this.save_any_fe

	ans.slope_method = this.slope_method
	ans.maxiter = this.maxiter
	ans.transform = this.transform
	ans.acceleration = this.acceleration
	ans.accel_start = this.accel_start
	ans.conlim = this.conlim
	ans.btol = this.btol
	ans.min_ok = this.min_ok
	ans.prune = this.prune
	ans.always_run_lsmr_preconditioner = this.always_run_lsmr_preconditioner

	return(ans)
}


// --------------------------------------------------------------------------
// Estimate finite condition number of the graph Laplacian
// --------------------------------------------------------------------------
`Void' FixedEffects::estimate_cond()
{
	`Matrix'                D1, D2, L
	`Vector'                lambda
	`RowVector'             tmp
	`Factor'                F12

	if (finite_condition!=-1) return

	if (verbose > 0) printf("\n{txt}## Computing finite condition number of the Laplacian\n\n")

	if (verbose > 0) printf("{txt}   - Constructing vectors of levels\n")
	F12 = join_factors(factors[1], factors[2], ., ., 1)
	
	// Non-sparse (lots of memory usage!)
	if (verbose > 0) printf("{txt}   - Constructing design matrices\n")
	D1 = designmatrix(F12.keys[., 1])
	D2 = designmatrix(F12.keys[., 2])
	assert_msg(rows(D1)<=2000, "System is too big!")
	assert_msg(rows(D2)<=2000, "System is too big!")

	if (verbose > 0) printf("{txt}   - Constructing graph Laplacian\n")
	L =   D1' * D1 , - D1' * D2 \
		- D2' * D1 ,   D2' * D2
	if (verbose > 0) printf("{txt}   - L is %g x %g \n", rows(L), rows(L))
		
	if (verbose > 0) printf("{txt}   - Computing eigenvalues\n")
	assert_msg(rows(L)<=2000, "System is too big!")
	eigensystem(L, ., lambda=.)
	lambda = Re(lambda')

	if (verbose > 0) printf("{txt}   - Selecting positive eigenvalues\n")
	lambda = edittozerotol(lambda, 1e-8)
	tmp = select(lambda,  edittozero(lambda, 1))
	tmp = minmax(tmp)
	finite_condition = tmp[2] / tmp[1]

	if (verbose > 0) printf("{txt}   - Finite condition number: {res}%s{txt}\n", strofreal(finite_condition))
}


`Real' FixedEffects::lsmr_norm(`Matrix' x)
{
	assert(cols(x)==1 | rows(x)==1)
	// BUGBUG: what if we have a corner case where there are as many obs as params?
	if (has_weights & cols(x)==1 & rows(x)==rows(weight)) {
		return(sqrt(quadcross(x, weight, x)))
	}
	else if (cols(x)==1) {
		return(sqrt(quadcross(x, x)))
	}
	else {
		return(sqrt(quadcross(x', x')))
	}
}


// Ax: given the coefs 'x', return the projection 'Ax'
`Vector' FixedEffects::lsmr_A_mult(`Vector' x)
{
	`Integer' g, k, idx_start, idx_end, i
	`Vector' ans
	`FactorPointer'         pf

	ans = J(N, 1, 0)
	idx_start = 1

	for (g=1; g<=G; g++) {
		pf = &(factors[g])
		k = (*pf).num_levels

		if (intercepts[g]) {
			idx_end = idx_start + k - 1
			ans = ans + (x[|idx_start, 1 \ idx_end , 1 |] :* asarray((*pf).extra, "precond_intercept") )[(*pf).levels, .]
			idx_start = idx_end + 1
		}

		if (num_slopes[g]) {
			for (i=1; i<=num_slopes[g]; i++) {
				idx_end = idx_start + k - 1
				ans = ans + x[|idx_start, 1 \ idx_end , 1 |][(*pf).levels] :* asarray((*pf).extra, "precond_slopes")
				idx_start = idx_end + 1
			}
		}

	}
	//assert(!missing(ans))
	return(ans)
}


// A'x: Compute the FEs and store them in a big stacked vector
`Vector' FixedEffects::lsmr_At_mult(`Vector' x)
{
	`Integer' m, g, i, idx_start, idx_end, k
	`Vector' ans
	`FactorPointer'         pf
	`Vector' alphas
	`Matrix' tmp_alphas

	alphas = J(M, 1, .)
	idx_start = 1

	for (g=1; g<=G; g++) {
		pf = &(factors[g])
		k = (*pf).num_levels

		if (intercepts[g]) {
			idx_end = idx_start + k - 1
			if (has_weights) {
				alphas[| idx_start , 1 \ idx_end , 1 |] = `panelsum'((*pf).sort(x :* weight), (*pf).info) :* asarray((*pf).extra, "precond_intercept")
			}
			else {
				alphas[| idx_start , 1 \ idx_end , 1 |] = `panelsum'((*pf).sort(x), (*pf).info) :* asarray((*pf).extra, "precond_intercept")
			}
			idx_start = idx_end + 1
		}

		if (num_slopes[g]) {
			idx_end = idx_start + k * num_slopes[g] - 1
			if (has_weights) {
				tmp_alphas = `panelsum'((*pf).sort(x :* weight :* asarray((*pf).extra, "precond_slopes")), (*pf).info)
			}
			else {
				tmp_alphas = `panelsum'((*pf).sort(x :* asarray((*pf).extra, "precond_slopes")), (*pf).info)
			}
			alphas[| idx_start , 1 \ idx_end , 1 |] = vec(tmp_alphas)
			idx_start = idx_end + 1
		}
	}
	//assert(!missing(alphas))
	return(alphas)
}

end

// --------------------------------------------------------------------------
// FixedEffects constructor (also precomputes factors)
// --------------------------------------------------------------------------

mata:

`FixedEffects' fixed_effects(`Varlist' absvars,
						   | `Varname' touse,
							 `String' weighttype,
							 `Varname' weightvar,
							 `Boolean' drop_singletons,
							 `Boolean' verbose)
{
	`FixedEffects'          S
	`Varname'               absvar, cvars
	`Integer'               i, j, g, gg, remaining
	`Vector'                idx
	`Integer'               spaces
	`Integer'               num_singletons_i
	`Variables'             cvar_data
	`FactorPointer'         pf

	// Set default value of arguments
	if (args()<2) touse = ""
	if (args()<3) weighttype = ""
	if (args()<4) weightvar = ""
	if (args()<5 | drop_singletons==.) drop_singletons = 1
	if (args()<6 | verbose==.) verbose = 0
	
	S = FixedEffects()
	S.verbose = verbose
	S.drop_singletons = drop_singletons

	// Parse absvars
	if (S.verbose > 0) printf("\n{txt}## Parsing absvars and HDFE options\n")
	
	if (touse == "") touse = st_tempname()
	st_global("reghdfe_touse", touse)
	stata(`"reghdfe5_parse "' + absvars)
	S.sample = `selectindex'(st_data(., touse))
	S.tousevar = touse // useful if later on we want to clone the HDFE object
	st_global("reghdfe_touse", "")

	if (st_global("s(residuals)") != "") S.residuals = st_global("s(residuals)")
	if (st_global("s(verbose)")!="") S.verbose = verbose = strtoreal(st_global("s(verbose)"))
	if (st_global("s(drop_singletons)")!="") S.drop_singletons = drop_singletons = strtoreal(st_global("s(drop_singletons)"))
	assert(S.verbose < .)
	assert(S.drop_singletons==0 | S.drop_singletons==1)

	if (S.verbose > 0) stata("sreturn list")
	S.G = strtoreal(st_global("s(G)"))
	S.absorb = absvars // useful if later on we want to clone the HDFE object
	S.absvars = tokens(st_global("s(absvars)"))
	S.has_intercept = strtoreal(st_global("s(has_intercept)"))
	S.save_any_fe = strtoreal(st_global("s(save_any_fe)"))
	S.save_all_fe = strtoreal(st_global("s(save_all_fe)"))
	S.ivars = tokens(st_global("s(ivars)"))
	S.cvars = tokens(st_global("s(cvars)"))
	S.targets = strtrim(tokens(st_global("s(targets)")))
	S.intercepts = strtoreal(tokens(st_global("s(intercepts)")))
	S.num_slopes = strtoreal(tokens(st_global("s(num_slopes)")))
	S.save_fe = S.targets :!= ""
	S.report_constant = strtoreal(st_global("s(report_constant)"))
	S.always_run_lsmr_preconditioner = strtoreal(st_global("s(precondition)"))

	// Ensure that S.report_constant and S.has_intercept are 0/1
	assert(anyof((0,1), S.has_intercept))
	assert(anyof((0,1), S.report_constant))
	S.compute_constant = S.has_intercept & S.report_constant

	if (st_global("s(tolerance)") != "") S.tolerance = strtoreal(st_global("s(tolerance)"))
	if (st_global("s(maxiter)") != "") S.maxiter = strtoreal(st_global("s(maxiter)"))
	if (st_global("s(prune)") != "") S.prune = strtoreal(st_global("s(prune)"))
	if (st_global("s(transform)") != "") S.transform = st_global("s(transform)")
	if (st_global("s(acceleration)") != "") S.acceleration = st_global("s(acceleration)")

	// Override LSMR if G=1
	if (S.G==1 & S.acceleration=="lsmr") S.acceleration = "conjugate_gradient"

	S.dofadjustments = tokens(st_global("s(dofadjustments)"))
	S.groupvar = st_global("s(groupvar)")
	if (st_global("s(finite_condition)")=="1") S.finite_condition = -1 // signal to compute it
	S.compute_rre = (st_global("s(compute_rre)")=="1")
	if (S.compute_rre) S.rre_varname = st_global("s(rre)")
	
	S.poolsize = strtoreal(st_global("s(poolsize)"))

	if (S.verbose > -1 & !S.has_intercept) printf("{txt}(warning: no intercepts terms in absorb(); regression lacks constant term)\n")

	S.extended_absvars = tokens(st_global("s(extended_absvars)"))
	S.tss = .

	assert(1<=S.G)
	if (S.G>10) printf("{txt}(warning: absorbing %2.0f dimensions of fixed effects; check that you really want that)\n", S.G)
	assert(S.G == cols(S.ivars))
	assert(S.G == cols(S.cvars))
	assert(S.G == cols(S.targets))
	assert(S.G == cols(S.intercepts))
	assert(S.G == cols(S.num_slopes))

	// Fill out object
	S.G = cols(S.absvars)
	S.factors = Factor(S.G)

	assert_msg(anyof(("", "fweight", "pweight", "aweight", "iweight"), weighttype), "wrong weight type")
	S.weight_type = weighttype
	S.weight_var = weightvar

	S.num_singletons = 0
	if (drop_singletons) {
		num_singletons_i = 0
		if (weighttype=="fweight" | weighttype=="iweight") {
			S.weight = st_data(S.sample, weightvar) // just to use it in F.drop_singletons()
		}
	}


	// (1) create the factors and remove singletons
	remaining = S.G
	i = 0
	if (S.verbose > 0) {
		printf("\n{txt}## Initializing Mata object for %g fixed effects\n\n", S.G)
		spaces = max((0, max(strlen(S.absvars))-4))
		printf("{txt}   {c TLC}{hline 4}{c TT}{hline 3}{c TT}{hline 1}%s{hline 6}{c TT}{hline 6}{c TT}{hline 9}{c TT}{hline 11}{c TT}{hline 12}{c TT}{hline 9}{c TT}{hline 14}{c TRC}\n", "{hline 1}" * spaces)
		printf("{txt}   {c |}  i {c |} g {c |} %s Name {c |} Int? {c |} #Slopes {c |}    Obs.   {c |}   Levels   {c |} Sorted? {c |} #Drop Singl. {c |}\n", " " * spaces)
		printf("{txt}   {c LT}{hline 4}{c +}{hline 3}{c +}{hline 1}%s{hline 6}{c +}{hline 6}{c +}{hline 9}{c +}{hline 11}{c +}{hline 12}{c +}{hline 9}{c +}{hline 14}{c RT}\n", "{hline 1}" * spaces)
		displayflush()
	}

	while (remaining) {
		++i
		g = 1 + mod(i-1, S.G)
		absvar = S.absvars[g]
		
		if (S.verbose > 0) {
			printf("{txt}   {c |} %2.0f {c |} %1.0f {c |} {res}%s{txt} {c |} ", i, g, (spaces+5-strlen(absvar)) * " " + absvar)
			printf("{txt}{%s}%3s{txt}  {c |}    %1.0f    {c |}", S.intercepts[g] ? "txt" : "err", S.intercepts[g] ? "Yes" : "No", S.num_slopes[g])
			displayflush()
		}

		if (S.verbose > 0) {
			printf("{res}%10.0g{txt} {c |}", rows(S.sample))
			displayflush()
		}

		if (rows(S.sample) < 2) {
			if (S.verbose > 0) printf("\n")
			exit(error(2001))
		}

		if (i<=S.G) {
			if (S.ivars[g] == "_cons" & S.G == 1) {
				// Special case without any fixed effects

				S.factors[g] = Factor()
				pf = &(S.factors[g])
				(*pf).num_obs = (*pf).counts = rows(S.sample)
				(*pf).num_levels = 1
				//(*pf).levels = . // Not filled to save space
				(*pf).levels = J(rows(S.sample), 1, 1)
				(*pf).is_sorted = 1
				(*pf).method = "none"

				// The code below is equivalent but 3x slower
				// S.factors[g] = _factor(J(rows(S.sample),1,1), 1, ., "hash0", ., 1, ., 0)
			}
			else {
				// We don't need to save keys (or sort levels but that might change estimates of FEs)
				S.factors[g] = factor(S.ivars[g], S.sample, ., "", ., 1, ., 0)
			}
		}

		if (S.verbose > 0) {
			printf(" {res}%10.0g{txt} {c |} %7s {c |}", S.factors[g].num_levels, S.factors[g].is_sorted ? "Yes" : "No")
			displayflush()
		}
 
		if (drop_singletons) {
			
			if (weighttype=="fweight") {
				idx = S.factors[g].drop_singletons(S.weight)
			}
			else if (weighttype=="iweight") {
				idx = S.factors[g].drop_singletons(S.weight, 1) // zero_threshold==1
			}
			else {
				idx = S.factors[g].drop_singletons()
			}

			num_singletons_i = rows(idx)
			S.num_singletons = S.num_singletons + num_singletons_i
			if (S.verbose > 0) {
				printf(" %10.0g   {c |}", num_singletons_i)
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
			if (S.verbose > 0) printf("      n/a     {c |}")
			--remaining
		}
		if (S.verbose > 0) printf("\n")
	}
	if (S.verbose > 0) {
		printf("{txt}   {c BLC}{hline 4}{c BT}{hline 3}{c BT}{hline 1}%s{hline 6}{c BT}{hline 6}{c BT}{hline 9}{c BT}{hline 11}{c BT}{hline 12}{c BT}{hline 9}{c BT}{hline 14}{c BRC}\n", "{hline 1}" * spaces)
	}

	if ( drop_singletons & S.num_singletons>0 & S.verbose>-1 | S.factors[1].num_obs<2) {
		if (weighttype=="iweight") {
			// PPML-specific
			printf(`"{txt}(dropped %s observations that are either {browse "http://scorreia.com/research/singletons.pdf":singletons} or {browse "http://scorreia.com/research/separation.pdf":separated} by a fixed effect)\n"', strofreal(S.num_singletons))
		}
		else {
			printf(`"{txt}(dropped %s {browse "http://scorreia.com/research/singletons.pdf":singleton observations})\n"', strofreal(S.num_singletons))
		}
	}

	if (S.factors[1].num_obs < 2) {
		exit(error(2001))
	}

	S.N = S.factors[1].num_obs // store number of obs.
	assert(S.N = S.factors[S.G].num_obs)
	assert(S.N > 1)


	// (2) run *.panelsetup() after the sample is defined
	if (S.verbose > 0) printf("\n{txt}## Initializing panelsetup() for each fixed effect\n\n")
	for (g=1; g<=S.G; g++) {
		absvar = S.absvars[g]
		if (S.verbose > 0) printf("{txt}   - panelsetup({res}%s{txt})\n", absvar)
		S.factors[g].panelsetup()
	}

	// (3) load cvars
	if (sum(S.num_slopes)) {
		if (S.verbose > 0) printf("\n{txt}## Loading slope variables\n\n")
		for (g=1; g<=S.G; g++) {
			cvars = tokens(S.cvars[g])
			if (S.num_slopes[g]) {
				// Load, standardize, sort by factor and store
				// Don't precompute aux objects (xmeans, inv_xx) as they depend on the weights
				// and will be computed on step (5)
				if (S.verbose > 0) printf("{txt}   - cvars({res}%s{txt})\n", invtokens(cvars))
				pf = &(S.factors[g])
				cvar_data = (*pf).sort(st_data(S.sample, cvars))
				asarray((*pf).extra, "x_stdevs", reghdfe_standardize(cvar_data))
				asarray((*pf).extra, "x", cvar_data)
			}
		}
		cvar_data = .
	}

	// (4) prune edges of degree-1
	// S.prune = 0 // bugbug
	if (S.prune) S.prune_1core()

	// (5) load weight
	S.load_weights(weighttype, weightvar, J(0,1,.), 1) // update S.has_weights, S.factors, etc.

	// Save "true" residuals for RRE
	if (S.compute_rre) {
		assert_msg(S.rre_varname != "")
		S.rre_true_residual = st_data(S.sample, S.rre_varname)
	}

	return(S)
}

end

// Common functions ---------------------------------------------------------
mata:

// --------------------------------------------------------------------------
// BUGBUG: not sure if this is still used...
// --------------------------------------------------------------------------
`StringRowVector' clean_tokens(`String' vars)
{
	`StringRowVector' 		ans
	`Integer'				i
	ans = tokens(vars)
	for (i=1; i<=cols(ans); i++) {
		ans[i] = invtokens(tokens(ans[i]))
	}
	return(ans)
}


// --------------------------------------------------------------------------
// Workaround to st_data's odd behavior
// --------------------------------------------------------------------------
// This does three things:
// 1) Wrap up interactions in parens (up to four) to avoid Stata's quirk/bug
// 2) If issue persists, load columns one-by-one
// 1) Instead of returning it reuses existing matrices (might use less mem?)
//
// Example of the issue:
// 	sysuse auto, clear
//	mata: cols(st_data(., "1.rep78 2.rep78 3.rep78#1.foreign")) // expected 3, got 6
// Happens b/c st_data doesn't work variable by variable but expands the interactions
// We can partly fix it by surrounding interactions with parens
// But st_data() only supports up to 4 parens


`Void' _st_data_wrapper(`Variables' index, `StringRowVector' vars, `Variables' data, `Boolean' verbose)
{
	`RowVector'			is_interaction
	`StringRowVector'	fixed_vars
	`Integer'			i, k

	vars = tokens(invtokens(vars))

	// Add parenthesis only for Stata 11-14, as on Stata 15+ they are i) not needed and ii) corrupt output
	// For i) see "help set fvtrack"
	// For ii) see "test/stdata3.do" on Github
	if (st_numscalar("c(stata_version)") < 15) {
		is_interaction = strpos(vars, "#") :> 0
		is_interaction = is_interaction :& (runningsum(is_interaction) :<= 4) // Only up to 4 parenthesis supported
		fixed_vars = subinstr(strofreal(is_interaction), "0", "")
		fixed_vars = subinstr(fixed_vars, "1", "(") :+ vars :+ subinstr(fixed_vars, "1", ")")
	}
	else {
		fixed_vars = vars
	}

	// Override code above, to minimize any risk of incorrect data
	// Since this is an undocumented feature, it might or might not work on some older versions of Stata
	// (See also email from jpitblado@stata.com)
	fixed_vars = vars

	data = st_data(index, fixed_vars)
	k = cols(vars)

	if (cols(data) > k) {
	    if (verbose > 0) printf("{err}(some empty columns were added due to a bug/quirk in {bf:st_data()}; %g cols created instead of %g for {it:%s}; running slower workaround)\n", cols(data), k, invtokens(vars))
	    data = J(rows(data), 0, .)
	    for (i=1; i<=k; i++) {
	        data = data, st_data(index, vars[i])
	    }
	}
	assert(cols(data)==k)
}
	

// --------------------------------------------------------------------------
// Each col of A will have stdev of 1 unless stdev is quite close to 0
// --------------------------------------------------------------------------
`RowVector' function reghdfe_standardize(`Matrix' A)
{
	`RowVector'				stdevs, means
	`Integer'				K, N // i, 

	// We don't need to good accuracy for the stdevs, so we have a few alternatives:
	// Note: cross(1,A) is the same as colsum(A), but faster
	// Note: cross(A, A) is very fast, but we only need the main diagonals
	// [A: 1sec] stdevs = sqrt( (colsum(A:*A) - (cross(1, A) :^ 2 / N)) / (N-1) )
	// [B: .61s] stdevs = sqrt( (diagonal(cross(A, A))' - (cross(1, A) :^ 2 / N)) / (N-1) )
	// [C: .80s] stdevs = diagonal(sqrt(variance(A)))'
	// [D: .67s] means = cross(1, A) / N; stdevs =  sqrt(diagonal(crossdev(A, means, A, means))' / (N-1))

	assert_msg(!isfleeting(A), "input cannot be fleeting")
	N = rows(A)
	K = cols(A)

	stdevs = J(1, K, .)

	// (A) Very precise

	// (B) Precise
	// means = cross(1, A) / N
	// stdevs =  sqrt(diagonal(quadcrossdev(A, means, A, means))' / (N-1))

	// (C) 20% faster; don't use it if you care about accuracy
	stdevs = sqrt( (diagonal(cross(A, A))' - (cross(1, A) :^ 2 / N)) / (N-1) )
	assert_msg(!missing(stdevs), "stdevs are missing; is N==1?") // Shouldn't happen as we don't expect N==1
	stdevs = colmax(( stdevs \ J(1, K, 1e-3) ))
	A = A :/ stdevs

	// (D) Equilibrate matrix columns instead of standardize (i.e. just divide by column max)
	// _perhapsequilc(A, stdevs=.)
	// stdevs = 1 :/ stdevs
	// assert_msg(!missing(stdevs), "stdevs are missing; is N==1?")

	// (E) Don't do anything
	// stdevs = J(1, cols(A), 1)

	return(stdevs)
}


// --------------------------------------------------------------------------
// Divide two row vectors but adjust the denominator if it's too small
// --------------------------------------------------------------------------
`RowVector' safe_divide(`RowVector' numerator, `RowVector' denominator, | `Real' epsi) {
	 // If the numerator goes below machine precision, we lose accuracy
	 // If the denominator goes below machine precision, the division explodes
	 if (args()<3 | epsi==.) epsi = epsilon(1)
	return( numerator :/ colmax(denominator \ J(1,cols(denominator),epsi)) )
}


// If X is not square...
// `Matrix' R
// real colvector tau, p

// _hqrdp(A, tau, R, p=.)
// B = hqrdmultq1t(A, tau, B)
// rank = _solveupper(R, B, tol)
// B = B[invorder(p),.]
// +- +-

// invsym(makesymmetric(..))





// --------------------------------------------------------------------------
// Robust solver for Ax=b
// --------------------------------------------------------------------------
// Mata utility for sequential use of solvers
// Default is cholesky;
// if that fails, use QR;
// if overridden, use QR.
// Author: Schaffer, Mark E <M.E.Schaffer@hw.ac.uk>
// --------------------------------------------------------------------------
// Warning:
// cholqrsolve calls qrsolve which calls _qrsolve which calls ...
// Does all the indirection makes it too slow to use within a panel?
// --------------------------------------------------------------------------
`Matrix' function reghdfe_cholqrsolve(`Matrix' A,
                                      `Matrix' B,
                                    | `Boolean' useqr)
{
	`Matrix' C
	if (args()<3 | useqr==.) useqr = 0
	
	if (!useqr) {
		C = cholsolve(A, B)
		if (hasmissing(C)) useqr = 1
	}

	if (useqr) {
		C = qrsolve(A, B)
	}

	return(C)
}


// --------------------------------------------------------------------------
// OLS Regression
// --------------------------------------------------------------------------
`Void' function reghdfe_post_ols(`FixedEffects' S,
                                 `Variables' X,
                                 `String' bname,
                                 `String' Vname,
                                 `String' nname,
                                 `String' rname,
                                 `String' dfrname)
{
	`Integer'				N
	`Integer'				rank
	`Integer'				df_r
	`Vector'				b
	`Matrix'				V
	`Variable'				resid
	`Real'					eps
	`Integer'				i
	`RowVector'				kept
	`Vector'				not_basevar


	`Vector'				idx
	`Vector'				temp_b
	`Matrix'				temp_V
	`Integer'				k

	if (S.timeit) timer_on(90)
	reghdfe_solve_ols(S, X, b=., V=., N=., rank=., df_r=., resid=., kept=., "vce_small")
	assert(cols(X) - 1 == rows(b) - S.compute_constant) // The 1st column of X is actually Y
	assert((rows(b) == rows(V)) & (rows(b) == cols(V)))
	if (S.timeit) timer_off(90)

	// Add base vars
	if (S.compute_constant) {
		if (S.verbose > 1) printf("\n{txt}## Adding _cons to varlist\n")
		assert_msg(rows(S.not_basevar) == 1, "rows(S.not_basevar) == 1")
		S.not_basevar = S.not_basevar, 1
		S.fullindepvars = S.fullindepvars + " _cons"
		S.indepvars = S.indepvars + " _cons"
	}
	if (S.not_basevar != J(1, 0, .)) {
		if (S.verbose > 1) printf("\n{txt}## Adding base variables to varlist\n")
		k = cols(S.not_basevar)
		assert_msg(cols(S.not_basevar) == k, "cols(S.not_basevar) == k")
		idx = `selectindex'(S.not_basevar)
		swap(b, temp_b)
		swap(V, temp_V)
		b = J(k, 1, 0)
		V = J(k, k, 0)
		b[idx, 1] = temp_b
		V[idx, idx] = temp_V
	}

	st_matrix(bname, b')

	if (S.verbose > 1) printf("\n{txt}## Reporting omitted variables\n")
	// Add "o." prefix to omitted regressors
	eps = sqrt(epsilon(1))
	for (i=1; i<=rows(b); i++) {
		if (b[i]==0 & S.not_basevar[i] & S.verbose > -1) {
			printf("{txt}note: %s omitted because of collinearity\n", tokens(S.fullindepvars)[i])
			//stata(sprintf("_ms_put_omit %s", indepvars[i]))
			//indepvars[i] = st_global("s(ospec)")
			// This is now one in reghdfe.ado with -_ms_findomitted-
		}
	}

	st_matrix(Vname, V)
	st_numscalar(nname, N)
	st_numscalar(rname, rank)
	st_numscalar(dfrname, df_r)

	// Need to save resids if saving FEs, even if temporarily
	if (S.residuals == "" & S.save_any_fe) {
		S.residuals = "__temp_reghdfe_resid__"
	}

	if (S.residuals != "") {
		if (S.verbose > 0) printf("\n{txt}## Storing residuals in {res}%s{txt}\n\n", S.residuals)
		if (S.compact == 1) {
			S.residuals_vector = resid
		}
		else {
			S.save_variable(S.residuals, resid, "Residuals")
		}
	}
}


`Void' function reghdfe_solve_ols(`FixedEffects' S,
                                  `Variables' X,
                                  `Vector' b,
                                  `Matrix' V,
                                  `Integer' N,
                                  `Integer' rank,
                                  `Integer' df_r,
                                  `Vector' resid,
                                  `RowVector' kept,
                                  `String' vce_mode,
                                | `Variable' true_w)
{
	// Hack: the first col of X is actually y!
	`Integer'				K, KK, tmp_N
	`Matrix'				xx, inv_xx, W, inv_V, just_X
	`Vector' 				w
	`Integer'				used_df_r
	`Integer'				dof_adj
	
	`Boolean'				is_standardized
	`Real'					stdev_y
	`RowVector'				stdev_x

	if (true_w == . | args() < 11) true_w = J(0, 1, .)
	if (S.vcetype == "unadjusted" & S.weight_type=="pweight") S.vcetype = "robust"
	if (S.verbose > 0) printf("\n{txt}## Solving least-squares regression of partialled-out variables\n\n")
	assert_in(vce_mode, ("vce_none", "vce_small", "vce_asymptotic"))

	is_standardized = S.all_stdevs != J(1, 0, .)
	if (is_standardized) S.means = S.means :/ S.all_stdevs

	// Weight FAQ:
	// - fweight: obs. i represents w[i] duplicate obs. (there is no loss of info wrt to having the "full" dataset)
	// - aweight: obs. i represents w[i] distinct obs. that were mean-collapsed (so there is loss of info and hetero)
	//	 soln: normalize them so they sum to N (the true number of obs in our sample), and then treat them as fweight
	// - pweight: each obs. represents only one obs. from the pop, that was drawn from w[i] individuals
	//	          we want to make inference on the population, so if we interviewed 100% of the men and only 10% of women,
	//            then without weighting we would be over-representing men, which leads to a loss of efficiency +-+-
	// 			  it is the same as aweight + robust
	// We need to pick N and w
	N = rows(X) // Default; will change with fweights
	S.sumweights = S.weight_type != "" ? quadsum(S.weight) : N
	assert(rows(S.means) == 1)
	assert(cols(S.means) == cols(X))

	w = 1
	if (rows(true_w)) {
		// Custom case for IRLS (ppmlhdfe) where S.weight = mu * true_w
		assert_msg(S.weight_type == "aweight")
		N = sum(true_w)
		w = S.weight * sum(true_w) / sum(S.weight)
	}
	else if (S.weight_type=="fweight") {
		N = S.sumweights
		w = S.weight
	}
	else if (S.weight_type=="aweight" | S.weight_type=="pweight") {
		w = S.weight * (N / S.sumweights)
	}

	// Build core matrices
	if (S.timeit) timer_on(91)

	K = cols(X) - 1
	xx = quadcross(X, w, X)
	S.tss_within = xx[1,1]
	xx = K ? xx[| 2 , 2 \ K+1 , K+1 |] : J(0, 0, .)
	if (S.timeit) timer_off(91)

	// This matrix indicates what regressors are not collinear
	assert_msg(cols(S.kept)==K+1, "partial_out() was run with a different set of vars")

	// Bread of the robust VCV matrix
	// Compute this early so we can update the list of collinear regressors
	if (S.timeit) timer_on(95)
	assert_msg( cols(tokens(invtokens(S.indepvars)))==cols(xx) , "HDFE.indepvars is missing or has the wrong number of columns")
	inv_xx = reghdfe_rmcoll(tokens(invtokens(S.indepvars)), xx, kept) // this modifies -kept-

	// // Workaround for case with extremely high weights, where ivnsym loses precision and incorrectly excludes vars
	// if (S.has_weights) {
	// 	if (max(S.weight) > 1e5) {
	// 		kept = (1..K)
	// 	}
	// }

	S.df_m = rank = K - diag0cnt(inv_xx)
	KK = S.df_a + S.df_m
	S.df_r = N - KK // replaced when clustering
	if (S.timeit) timer_off(95)

	// Compute betas
	// - There are two main options
	//	 a) Use cholqrsolve on xx and xy. Faster but numerically inaccurate
	//      See: http://www.stata.com/statalist/archive/2012-02/msg00956.html
	//   b) Use qrsolve. More accurate but doesn't handle weights easily
	// - Ended up doing (b) with a hack for weights
	b = J(K, 1, 0)
	if (cols(kept)) {
		if (S.has_weights) {
			b[kept] = qrsolve(X[., 1:+kept] :* sqrt(S.weight), X[., 1] :* sqrt(S.weight))
		}
		else {
			b[kept] = qrsolve(X[., 1:+kept], X[., 1])
		}
	}

	if (S.timeit) timer_on(92)
	if (!isfleeting(resid) | vce_mode != "vce_none") resid = X * (1 \ -b) // y - X * b
	if (S.timeit) timer_off(92)

	if (S.compute_constant) {
		tmp_N = (S.weight_type=="aweight" | S.weight_type=="pweight") ? N : S.sumweights
		if (rows(true_w)) tmp_N = N
		reghdfe_extend_b_and_inv_xx(S.means, tmp_N, b, inv_xx)
	}

	// Stop if no VCE/R2/RSS needed
	if (vce_mode == "vce_none") {
		assert(!is_standardized)
		return
	}

	if (S.timeit) timer_on(93)
	if (S.vcetype != "unadjusted") {
		if (S.compute_constant) {
			if (isfleeting(X)) {
				// Save some memory... unsure if it helps
				swap(just_X, X)
				just_X = K ? just_X[., 2..K+1] :+ S.means[2..cols(S.means)] : J(rows(just_X), 0, .)
			}
			else {
				just_X = K ? X[., 2..K+1] :+ S.means[2..cols(S.means)] : J(rows(X), 0, .)
			}
		}
		else {
			just_X = K ? X[., 2..K+1] : J(rows(X), 0, .)
		}
	}
	if (S.timeit) timer_off(93)

	if (S.timeit) timer_on(94)
	S.rss = quadcross(resid, w, resid) // do before reghdfe_robust() modifies w
	if (S.timeit) timer_off(94)

	// Compute full VCE
	if (S.timeit) timer_on(96)
	assert_msg(anyof( ("unadjusted", "robust", "cluster") , S.vcetype), "invalid vcetype" + S.vcetype)
	if (S.vcetype == "unadjusted") {
		if (S.verbose > 0) {
			printf("{txt}   - Small-sample-adjustment: q = N / (N-df_m-df_a) = %g / (%g - %g - %g) = %g\n", N, N, rank, S.df_a, N / S.df_r )
		}
		dof_adj = N / S.df_r
		if (vce_mode == "vce_asymptotic") dof_adj = N / (N-1) // 1.0
		V = (S.rss / N) * dof_adj * inv_xx
	}
	else if (S.vcetype == "robust") {
		V = reghdfe_robust(S, just_X, inv_xx, resid, w, N, KK, vce_mode, true_w)
	}
	else {
		V = reghdfe_cluster(S, just_X, inv_xx, resid, w, N, KK, vce_mode)
	}
	if (S.timeit) timer_off(96)

	// Wald test: joint significance
	if (S.timeit) timer_on(97)
	inv_V = invsym(V[kept, kept]) // this might not be of full rank but numerical inaccuracies hide it
	if (diag0cnt(inv_V)) {
		if (S.verbose > -1) printf("{txt}warning: missing F statistic; dropped variables due to collinearity or too few clusters\n")
		W = .
	}
	else if (length(b[kept])==0) {
		W = .
	}
	else {
		// We could probably do this with the simpler formula instead of Wald
		W = b[kept]' * inv_V * b[kept] / S.df_m
		if (missing(W) & S.verbose > -1) printf("{txt}warning: missing F statistic\n")
	}
	if (S.timeit) timer_off(97)

	// V can be missing if b is completely absorbed by the FEs
	if (missing(V)) {
		if (S.verbose > 0) printf("{txt}   - VCE has missing values, setting it to zeroes (are your regressors all collinear?)\n")
		V = J(rows(V), rows(V), 0)
	}

	// Undo standardization
	if (is_standardized) {
		// Sanity checks
		assert(rows(S.all_stdevs)==1)
		assert(cols(S.all_stdevs) - 1 == rows(b) - S.compute_constant) // Subtract "y" on left; subtract "_cons" on right
		
		// Recover stdevs
		stdev_y = S.all_stdevs[1]
		stdev_x = K ? S.all_stdevs[2..cols(S.all_stdevs)] : J(1, 0, .)
		if (S.compute_constant) stdev_x = stdev_x, 1
		stdev_x = stdev_x :/ stdev_y
		
		// Transform output (note that S.tss is already ok)
		S.rss = S.rss * stdev_y ^ 2
		S.tss_within = S.tss_within * stdev_y ^ 2
		resid = resid * stdev_y
		V = V :/ (stdev_x' * stdev_x)
		b = b :/ stdev_x'
	}

	// Results
	S.title = "Linear regression"
	// S.model = "ols"
	used_df_r = N - KK - S.df_a_nested
	S.r2 = 1 - S.rss / S.tss
	S.r2_a = 1 - (S.rss / used_df_r) / (S.tss / (N - S.has_intercept ) )
	S.r2_within = 1 - S.rss / S.tss_within
	S.r2_a_within = 1 - (S.rss / used_df_r) / (S.tss_within / (used_df_r + rank))

	S.ll = - 0.5 * N * (1 + ln(2 * pi()) + ln(S.rss / N))
	S.ll_0 = - 0.5 * N * (1 + ln(2 * pi()) + ln(S.tss_within / N))

	S.rmse = sqrt(S.rss / used_df_r)
	if (used_df_r==0) S.rmse = sqrt(S.rss)
	S.F = W
	df_r = S.df_r // reghdfe_cluster might have updated it (this gets returned to the caller function)
}


// --------------------------------------------------------------------------
// Robust VCE
// --------------------------------------------------------------------------
// Advice: Delegate complicated regressions to -avar- and specialized routines
// BUGBUG: do we standardize X again? so V is well behaved?
// Notes:
// - robust is the same as cluster robust where cluster==_n
// - cluster just "collapses" X_i * e_i for each group, and builds M from that

`Matrix' reghdfe_robust(`FixedEffects' S,
                        `Variables' X,
						`Matrix' D,
						`Variable' resid,
						`Variable' w,
						`Integer' N,
						`Integer' K,
						`String' vce_mode,
					    `Variable' true_w)
{
	`Matrix'				M, V
	`Integer'				dof_adj

	if (S.verbose > 0) printf("\n{txt}## Estimating Robust Variance-Covariance Matrix of the Estimators (VCE)\n\n")
	if (S.verbose > 0) printf("{txt}   - VCE type: {res}%s{txt}\n", S.vcetype)
	if (S.verbose > 0) printf("{txt}   - Weight type: {res}%s{txt}\n", S.weight_type=="" ? "<none>" : S.weight_type)

	if (rows(true_w)) {
		assert(S.weight_type=="aweight")
		w = (resid :* w) :^ 2 :/ true_w // resid^2 * aw^2 * fw
	}
	else if (S.weight_type=="") {
		w = resid :^ 2
	}
	else if (S.weight_type=="fweight") {
		w = resid :^ 2 :* w
	}
	else if (S.weight_type=="aweight" | S.weight_type=="pweight") {
		w = (resid :* w) :^ 2
	}

	dof_adj = N / (N - K)
	if (vce_mode == "vce_asymptotic") dof_adj = N / (N-1) // 1.0
	M = S.compute_constant ? quadcross(X, 1, w, X, 1) : quadcross(X, w, X)
	if (S.verbose > 0) {
		printf("{txt}   - Small-sample-adjustment: q = N / (N-df_m-df_a) = %g / (%g - %g - %g) = %g\n", N, N, K-S.df_a, S.df_a, N / (N-K) )
	}
	V = D * M * D * dof_adj
	return(V)
}

`Matrix' reghdfe_cluster(`FixedEffects' S,
                        `Variables' X,
						`Matrix' D,
						`Variable' resid,
						`Variable' w,
						`Integer' N,
						`Integer' K,
						`String' vce_mode)
{
	`Matrix' 				M, V
	`Integer'				dof_adj, N_clust, df_r, nested_adj
	`Integer'				Q, q, g, sign, i, j
	pointer(`Factor') rowvector FPlist
	`FactorPointer'			FP
	`Varlist'				vars
	`String'				var, var_with_spaces
	`Boolean'				clustervar_is_absvar, required_fix
	`Matrix'				tuples
	`RowVector'				tuple
	`RowVector'				N_clust_list
	`Matrix'				joined_levels
	`Integer'				Msize

	w = resid :* w
	Msize = cols(X) + S.compute_constant

	vars = S.clustervars
	Q = cols(vars)
	if (S.verbose > 0) printf("\n{txt}## Estimating Cluster Robust Variance-Covariance Matrix of the Estimators (VCE)\n\n")
	if (S.verbose > 0) printf("{txt}   - VCE type: {res}%s{txt} (%g-way clustering)\n", S.vcetype, Q)
	if (S.verbose > 0) printf("{txt}   - Cluster variables: {res}%s{txt}\n", invtokens(vars))
	if (S.verbose > 0) printf("{txt}   - Weight type: {res}%s{txt}\n", S.weight_type=="" ? "<none>" : S.weight_type)
	assert_msg(0 < Q & Q < 10)

	// Get or build factors associated with the clustervars
	FPlist = J(1, Q, NULL)
	N_clust_list = J(1, Q, .)
	for (q=1; q<=Q; q++) {
		var = vars[q]
		clustervar_is_absvar = 0
		for (g=1; g<=S.G; g++) {
			if (invtokens(S.factors[g].varlist, "#") == var) {
				clustervar_is_absvar = 1
				FP = &(S.factors[g])
				break
			}
		}
		var_with_spaces = subinstr(var, "#", " ")
		if (!clustervar_is_absvar) FP = &(factor(var_with_spaces, S.sample, ., "", ., ., ., 0))
		N_clust_list[q] = (*FP).num_levels
		if (S.verbose > 0) printf("{txt}   - {res}%s{txt} has {res}%g{txt} levels\n", var, N_clust_list[q])
		FPlist[q] = FP
	}

	// Build the meat part of the V matrix
	if (S.verbose > 0) printf("{txt}   - Computing the 'meat' of the VCE\n")
	M = J(Msize, Msize, 0)
	tuples = .
	for (q=1; q<=Q; q++) {
		tuples = reghdfe_choose_n_k(Q, q, tuples)
		sign = mod(q, 2) ? 1 : -1 // + with odd number of variables, - with even
		for (j=1; j<=rows(tuples); j++) {
			tuple = tuples[j, .]
			if (S.verbose > 0) printf("{txt}      - Level %g/%g; sublevel %g/%g; M = M %s ClusterVCE(%s)\n", q, Q, j, rows(tuples), sign > 0 ? "+" : "-" , invtokens(strofreal(tuple)))
			if (q==1) {
				assert(tuple==j)
				FP =  FPlist[j]
			}
			else if (q==2) {
				FP = &join_factors( *FPlist[tuple[1]] , *FPlist[tuple[2]] , ., ., 1)
			}
			else {
				joined_levels = (*FPlist[tuple[1]]).levels
				for (i=2; i<=cols(tuple); i++) {
					joined_levels = joined_levels, (*FPlist[tuple[i]]).levels
				}
				FP = &_factor(joined_levels, ., ., "", ., ., ., 0)
			}
			M = M + sign * reghdfe_vce_cluster_meat(FP, X, w, Msize, S.compute_constant)
		}
	}

	// Build VCE
	N_clust = min(N_clust_list)
	
	nested_adj = S.df_a_nested > 0 // minor adj. so we match xtreg when the absvar is nested within cluster
	// (when ..nested.., df_a is zero so we divide N-1 by something that can potentially be N (!))
	// so we either add the 1 back, or change the numerator (and the N_clust-1 factor!)
	// addendum: this just ensures we subtract the constant when we have nested FEs

	dof_adj = (N - 1) / (N - nested_adj - K) * N_clust / (N_clust - 1) // adjust for more than 1 cluster
	if (vce_mode == "vce_asymptotic") dof_adj = N_clust / (N_clust - 1)  // 1.0
	if (S.verbose > 0) {
		printf("{txt}   - Small-sample-adjustment: q = (%g - 1) / (%g - %g) * %g / (%g - 1) = %g\n", N, N, K+nested_adj, N_clust, N_clust, dof_adj)
	}
	V = D * M * D * dof_adj
	if (Q > 1) {
		required_fix = reghdfe_fix_psd(V)
		if (required_fix) printf("{txt}Warning: VCV matrix was non-positive semi-definite; adjustment from Cameron, Gelbach & Miller applied.\n")
	}

	// Store e()
	assert(!missing(S.df_r))
	df_r = N_clust - 1
	if (S.df_r > df_r) {
		S.df_r = df_r
	}
	else if (S.verbose > 0) {
		printf("{txt}   - Unclustered df_r (N - df_m - df_a = %g) are {it:lower} than clustered df_r (N_clust-1 = %g)\n", S.df_r, df_r)
		printf("{txt}     Thus, we set e(df_r) as the former.\n")
		printf("{txt}     This breaks consistency with areg but ensures internal consistency\n")
		printf("{txt}     between vce(robust) and vce(cluster _n)\n")
	}
	
	S.N_clust = N_clust
	S.N_clust_list = N_clust_list

	return(V)
}


`Matrix' reghdfe_vce_cluster_meat(`FactorPointer' FP,
                                  `Variables' X,
                                  `Variable' resid,
                                  `Integer' Msize,
                                  `Boolean' compute_constant)
{
	`Integer'				i, N_clust
	`Variables'				X_sorted
	`Variable'				resid_sorted
	`Matrix'				X_tmp
	`Vector'				resid_tmp
	`RowVector'				Xe_tmp
	`Matrix'				M

	if (cols(X)==0 & !compute_constant) return(J(0,0,0))

	N_clust = (*FP).num_levels
	(*FP).panelsetup()
	X_sorted = (*FP).sort(X)
	resid_sorted = (*FP).sort(resid)
	M = J(Msize, Msize, 0)

	if (cols(X)) {
		for (i=1; i<=N_clust; i++) {
			X_tmp = panelsubmatrix(X_sorted, i, (*FP).info)
			resid_tmp = panelsubmatrix(resid_sorted, i, (*FP).info)
			Xe_tmp = quadcross(1, 0, resid_tmp, X_tmp, compute_constant) // Faster than colsum(e_tmp :* X_tmp)
			M = M + quadcross(Xe_tmp, Xe_tmp)
		}		
	}
	else {
		// Workaround for when there are no Xs except for _cons
		assert(compute_constant)
		for (i=1; i<=N_clust; i++) {
			resid_tmp = panelsubmatrix(resid_sorted, i, (*FP).info)
			M = M + quadsum(resid_tmp) ^ 2
		}
	}

	return(M)
}


// Enumerate all combinations of K integers from N integers
// Kroneker approach based on njc's tuples.ado
`Matrix' reghdfe_choose_n_k(`Integer' n, `Integer' k, `Matrix' prev_ans)
{
	`RowVector' v
	`Integer' q
	`Matrix' candidate
	`Matrix' ans
	v = 1::n
	if (k==1) return(v)

	q = rows(prev_ans)
	assert(q==comb(n, k-1))
	assert(cols(prev_ans)==k-1)
	candidate = v # J(q, 1, 1)
	candidate = candidate , J(n, 1, prev_ans)	
	ans = select(candidate, candidate[., 1] :< candidate[., 2])
	return(ans)
}


// --------------------------------------------------------------------------
// Fix non-positive VCV
// --------------------------------------------------------------------------
// If the VCV matrix is not positive-semidefinite, use the fix from
// Cameron, Gelbach & Miller - Robust Inference with Multi-way Clustering (JBES 2011)
// 1) Use eigendecomposition V = U Lambda U' where U are the eigenvectors and Lambda = diag(eigenvalues)
// 2) Replace negative eigenvalues into zero and obtain FixedLambda
// 3) Recover FixedV = U * FixedLambda * U'
`Boolean' function reghdfe_fix_psd(`Matrix' V) {
	`Matrix'				U
	`Matrix'				lambda
	`Boolean' 				required_fix

	if (!issymmetric(V)) _makesymmetric(V)
	if (!issymmetric(V)) exit(error(505))
	symeigensystem(V, U=., lambda=.)
	if (min(lambda)<0) {
		lambda = lambda :* (lambda :>= 0)
		// V = U * diag(lambda) * U'
		V = quadcross(U', lambda, U')
		required_fix = 1
	}
	else {
		required_fix = 0
	}
	return(required_fix)
}


// --------------------------------------------------------------------------
// Remove collinear variables
// --------------------------------------------------------------------------
// Based on ivreg2's s_rmcoll2
`Matrix' reghdfe_rmcoll(`Varlist' varnames,
                        `Matrix' xx,
                        `RowVector' kept)
{
	`Integer'				K, num_dropped
	`Matrix'				inv_xx, smat, alt_inv_xx
	`RowVector'				vl_drop, vl_keep

	assert(rows(xx)==cols(xx))
	K = cols(xx)
	inv_xx = K ? invsym(xx, 1..K) : J(0, 0, .)
	
	// Specifying the sweep order in invsym() can lead to incorrectly dropped regressors
	// (EG: with very VERY high weights)
	// We'll double check in this case
	num_dropped = diag0cnt(inv_xx)
	if (K & num_dropped) {
		alt_inv_xx = invsym(xx)
		if (num_dropped != diag0cnt(alt_inv_xx)) {
			inv_xx = alt_inv_xx
			num_dropped = diag0cnt(alt_inv_xx)
		}
	}

	st_numscalar("r(k_omitted)", num_dropped)
	smat = (diagonal(inv_xx) :== 0)'
	vl_drop = select(varnames, smat)
	vl_keep = select(varnames, !smat)
	if (cols(vl_keep)) st_global("r(varlist)", invtokens(vl_keep))
	if (cols(vl_drop)) st_global("r(omitted)", invtokens(vl_drop))
	kept = `selectindex'(!smat) // Return it, so we can exclude these variables from the joint Wald test
	return(inv_xx)
}


// --------------------------------------------------------------------------
// Use regression-through-mean and block partition formula to enlarge b and inv(XX)
// --------------------------------------------------------------------------
`Void' reghdfe_extend_b_and_inv_xx(
	`RowVector' 	means,
	`Integer'		N,
	`Vector'		b,
	`Matrix' 		inv_xx)
{
	// How to add back _cons:
	// 1) To recover coefficient, apply "regression through means formula":
	//	  b0 = mean(y) - mean(x) * b1

	// 2) To recover variance ("full_inv_xx")
	//	  apply formula for inverse of partitioned symmetric matrix
	//    http://fourier.eng.hmc.edu/e161/lectures/gaussianprocess/node6.html
	//    http://www.cs.nthu.edu.tw/~jang/book/addenda/matinv/matinv/
	//
	//    Given A = [X'X X'1]		B = [B11 B21']		B = inv(A)
	//              [1'X 1'1]			[B21 B22 ]
	//
	//	  B11 is just inv(xx) (because of Frisch-Waugh)
	//	  B21 ("side") = means * B11
	//	  B22 ("corner") = 1 / sumweights * (1 - side * means')
	//
	//	- Note that means is NOT A12, but A12/N or A12 / (sum_weights)

	//  - Note: aw and pw (and unweighted) use normal weights,
	//	  but for fweights we expected S.sumweights

	`RowVector'		means_x, side
	`Real'			corner

	means_x = cols(means) > 1 ? means[2..cols(means)] : J(1, 0, .)
	b = b \ means[1] - means_x * b // means * (1 \ -b)
	corner = (1 / N) + means_x * inv_xx * means_x'
	side = - means_x * inv_xx
	inv_xx = (inv_xx , side' \ side , corner)

}

end

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


mata:

// --------------------------------------------------------------------------
// Transformations: Compute RESIDUALS, not projections
// --------------------------------------------------------------------------

`Void' function transform_cimmino(`FixedEffects' S, `Variables' y, `Variables' ans,| `Boolean' get_proj) {
	`Integer' 	g
	if (args()<4 | get_proj==.) get_proj = 0
	ans = S.project_one_fe(y, 1)
	for (g=2; g<=S.G; g++) {
		ans = ans + S.project_one_fe(y, g)
	}
	ans = get_proj ? ans / S.G : y - ans / S.G
}

// --------------------------------------------------------------------------

`Void' function transform_kaczmarz(`FixedEffects' S, `Variables' y, `Variables' ans,| `Boolean' get_proj) {
	`Integer' 	g

	if (args()<4 | get_proj==.) get_proj = 0

	ans = y - S.project_one_fe(y, 1)
	for (g=2; g<=S.G; g++) {
		ans = ans - S.project_one_fe(ans, g)
	}
	if (get_proj) ans = y - ans
}

// --------------------------------------------------------------------------
// This seems slower than kaczmarz (sym kaczmarz!); not used currently
`Void' function transform_rand_kaczmarz(`FixedEffects' S, `Variables' y, `Variables' ans,| `Boolean' get_proj) {
	`Integer' 	g
	`Vector' rand
	if (args()<4 | get_proj==.) get_proj = 0
	rand = sort( ( (1::S.G) , uniform(S.G,1) ) , 2 )[.,1]
	ans = y - S.project_one_fe(y, rand[1])
	for (g=2; g<=S.G; g++) {
		ans = ans - S.project_one_fe(ans, rand[g])
	}
	for (g=S.G-1; g>=1; g--) {
		ans = ans - S.project_one_fe(ans, rand[g])
	}
	if (get_proj) ans = y - ans
}

// --------------------------------------------------------------------------

 `Void' function transform_sym_kaczmarz(`FixedEffects' S, `Variables' y, `Variables' ans,| `Boolean' get_proj) {
	`Integer' 	g
	if (args()<4 | get_proj==.) get_proj = 0
	if (S.timeit) timer_on(72)
	ans = y - S.project_one_fe(y, 1)
	if (S.timeit) timer_off(72)
	for (g=2; g<=S.G; g++) {
		if (S.timeit) timer_on(72)
		ans = ans - S.project_one_fe(ans, g)
		if (S.timeit) timer_off(72)
	}
	for (g=S.G-1; g>=1; g--) {
		if (S.timeit) timer_on(72)
		ans = ans - S.project_one_fe(ans, g)
		if (S.timeit) timer_off(72)
	}
	if (get_proj) ans = y - ans
}
end

mata:

// --------------------------------------------------------------------------
// Acceleration Schemes
// --------------------------------------------------------------------------

`Variables' function accelerate_test(`FixedEffects' S, `Variables' y, `FunctionP' T) {
	`Integer'	iter, g
	`Variables'		resid
	`Factor' f
	pragma unset resid

	assert(S.converged == 0)

	for (iter=1; iter<=S.maxiter; iter++) {
		for (g=1; g<=S.G; g++) {
			f = S.factors[g]
			if (g==1) resid = y - panelmean(f.sort(y), f)[f.levels, .]
			else resid = resid - panelmean(f.sort(resid), f)[f.levels, .]
		}
		if (check_convergence(S, iter, resid, y)) break
		y = resid
	}
	return(resid)
}

// --------------------------------------------------------------------------

`Variables' function accelerate_none(`FixedEffects' S, `Variables' y, `FunctionP' T) {
	`Integer'	iter
	`Variables'		resid
	pragma unset resid

	assert(S.converged == 0)

	for (iter=1; iter<=S.maxiter; iter++) {
		(*T)(S, y, resid) // Faster version of "resid = S.T(y)"
		if (check_convergence(S, iter, resid, y)) break
		y = resid
	}
	return(resid)
}
// --------------------------------------------------------------------------

// Start w/out acceleration, then switch to CG
`Variables' function accelerate_hybrid(`FixedEffects' S, `Variables' y, `FunctionP' T) {
	`Integer' iter, accel_start
	`Variables' resid
	pragma unset resid

	accel_start = 6
	assert(S.converged == 0)

	for (iter=1; iter<=accel_start; iter++) {
		(*T)(S, y, resid) // Faster version of "resid = S.T(y)"
		if (check_convergence(S, iter, resid, y)) break
		y = resid
	}

	T = &transform_sym_kaczmarz() // Override

	return(accelerate_cg(S, y, T))
}

// --------------------------------------------------------------------------
// Memory cost is approx = 4*size(y) (actually 3 since y is already there)
// But we need to add maybe 1 more due to u:*v
// And I also need to check how much does project and T use..
// Double check with a call to memory

// For discussion on the stopping criteria, see the following presentation:
// Arioli & Gratton, "Least-squares problems, normal equations, and stopping criteria for the conjugate gradient method". URL: https://www.stfc.ac.uk/SCD/resources/talks/Arioli-NAday2008.pdf

// Basically, we will use the Hestenes and Stiefel rule

`Variables' function accelerate_cg(`FixedEffects' S, `Variables' y, `FunctionP' T) {
	// BUGBUG iterate the first 6? without acceleration??
	`Integer'	iter, d, Q
	`Variables'		r, u, v
	`RowVector' alpha, beta, ssr, ssr_old, improvement_potential
	`Matrix' recent_ssr
	pragma unset r
	pragma unset v

	assert(S.converged == 0)
	if (S.timeit) timer_on(70)
	Q = cols(y)
	
	d = 1 // BUGBUG Set it to 2/3 // Number of recent SSR values to use for convergence criteria (lower=faster & riskier)
	// A discussion on the stopping criteria used is described in
	// http://scicomp.stackexchange.com/questions/582/stopping-criteria-for-iterative-linear-solvers-applied-to-nearly-singular-system/585#585

	if (S.timeit) timer_on(73)
	improvement_potential = weighted_quadcolsum(S, y, y)
	recent_ssr = J(d, Q, .)
	if (S.timeit) timer_off(73)
	
	if (S.timeit) timer_on(71)
	(*T)(S, y, r, 1)
	if (S.timeit) timer_off(71)
	if (S.timeit) timer_on(73)
	ssr = weighted_quadcolsum(S, r, r) // cross(r,r) when cols(y)==1 // BUGBUG maybe diag(quadcross()) is faster?
	u = r
	if (S.timeit) timer_off(73)

	for (iter=1; iter<=S.maxiter; iter++) {
		if (S.timeit) timer_on(71)
		(*T)(S, u, v, 1) // This is the hottest loop in the entire program
		if (S.timeit) timer_off(71)
		if (S.timeit) timer_on(73)
		alpha = safe_divide( ssr , weighted_quadcolsum(S, u, v) )
		if (S.timeit) timer_off(73)
		if (S.timeit) timer_on(74)
		recent_ssr[1 + mod(iter-1, d), .] = alpha :* ssr
		improvement_potential = improvement_potential - alpha :* ssr
		y = y - alpha :* u
		if (S.timeit) timer_off(74)
		if (S.timeit) timer_on(75)
		if (S.compute_rre & !S.prune) reghdfe_rre_benchmark(y[., 1], S.rre_true_residual, S.rre_depvar_norm)
		r = r - alpha :* v
		ssr_old = ssr
		if (S.timeit) timer_off(75)
		if (S.timeit) timer_on(73)
		if (S.verbose>=5) r
		ssr = weighted_quadcolsum(S, r, r)
		beta = safe_divide( ssr , ssr_old) // Fletcher-Reeves formula, but it shouldn't matter in our problem
		if (S.timeit) timer_off(73)
		u = r + beta :* u
		// Convergence if sum(recent_ssr) > tol^2 * improvement_potential
		if (S.timeit) timer_on(76)
		if ( check_convergence(S, iter, colsum(recent_ssr), improvement_potential, "hestenes") ) {
			break
			if (S.timeit) timer_off(76)
		}
		if (S.timeit) timer_off(76)
	}
	if (S.timeit) timer_off(70)
	return(y)
}

// --------------------------------------------------------------------------

`Variables' function accelerate_sd(`FixedEffects' S, `Variables' y, `FunctionP' T) {
	`Integer'	iter, g
	`Variables' proj
	`RowVector' t
	pragma unset proj

	assert(S.converged == 0)

	for (iter=1; iter<=S.maxiter; iter++) {
		(*T)(S, y, proj, 1)
		if (check_convergence(S, iter, y-proj, y)) break
		t = safe_divide( weighted_quadcolsum(S, y, proj) , weighted_quadcolsum(S, proj, proj) )
		if (uniform(1,1)<0.1) t = 1 // BUGBUG: Does this REALLY help to randomly unstuck an iteration?

		y = y - t :* proj
		if (S.compute_rre & !S.prune) reghdfe_rre_benchmark(y[., 1], S.rre_true_residual, S.rre_depvar_norm)

		if (S.storing_alphas) {
			for (g=1; g<=S.G; g++) {
				//g, ., ., t
				//asarray(S.factors[g].extra, "alphas"), asarray(S.factors[g].extra, "tmp_alphas")
				if (S.save_fe[g]) {
					asarray(S.factors[g].extra, "alphas",
					    asarray(S.factors[g].extra, "alphas") +
					    t :* asarray(S.factors[g].extra, "tmp_alphas")
					)
				}
			}
		}
	}
	return(y-proj)
}

// --------------------------------------------------------------------------
// This is method 3 of Macleod (1986), a vector generalization of the Aitken-Steffensen method
// Also: "when numerically computing the sequence.. stop..  when rounding errors become too 
// important in the denominator, where the ^2 operation may cancel too many significant digits"
// Note: Sometimes the iteration gets "stuck"; can we unstuck it with adding randomness
// in the accelerate decision? There should be a better way.. (maybe symmetric kacz instead of standard one?)

`Variables' function accelerate_aitken(`FixedEffects' S, `Variables' y, `FunctionP' T) {
	`Integer'	iter
	`Variables'		resid, y_old, delta_sq
	`Boolean'	accelerate
	`RowVector' t
	pragma unset resid

	assert(S.converged == 0)
	y_old = J(rows(y), cols(y), .)

	for (iter=1; iter<=S.maxiter; iter++) {
		
		(*T)(S, y, resid)
		accelerate = iter>=S.accel_start & !mod(iter,S.accel_freq)

		// Accelerate
		if (accelerate) {
			delta_sq = resid - 2 * y + y_old // = (resid - y) - (y - y_old) // Equivalent to D2.resid
			// t is just (d'd2) / (d2'd2)
			t = safe_divide( weighted_quadcolsum(S,  (resid - y) , delta_sq) ,  weighted_quadcolsum(S, delta_sq , delta_sq) )
			resid = resid - t :*  (resid - y)
		}


		// Only check converge on non-accelerated iterations
		// BUGBUG: Do we need to disable the check when accelerating?
		// if (check_convergence(S, iter, accelerate? resid :* .  : resid, y)) break
		if (S.compute_rre & !S.prune) reghdfe_rre_benchmark(resid[., 1], S.rre_true_residual, S.rre_depvar_norm)
		if (check_convergence(S, iter, resid, y)) break
		y_old = y // y_old is resid[iter-2]
		y = resid // y is resid[iter-1]
	}
	return(resid)
}

// --------------------------------------------------------------------------

`Boolean' check_convergence(`FixedEffects' S, `Integer' iter, `Variables' y_new, `Variables' y_old,| `String' method) {
	`Boolean'	is_last_iter
	`Real'		update_error
	`Real'		eps_threshold

	// max() ensures that the result when bunching vars is at least as good as when not bunching
	if (args()<5) method = "vectors"

	if (S.G==1 & !S.storing_alphas) {
		// Shortcut for trivial case (1 FE)
		update_error = 0
	}
	else if (method=="vectors") {
		update_error = max(mean(reldif(y_new, y_old), S.weight))
	}
	else if (method=="hestenes") {
		// If the regressor is perfectly explained by the absvars, we can have SSR very close to zero but negative
		// (so sqrt is missing)

		eps_threshold = 1e-15 // 10 * epsilon(1) ; perhaps too aggressive and should be 1e-14 ?
		if (S.verbose > 0 & all(y_new :< eps_threshold)) {
			printf("{txt} note: eps. is very close to zero (%g), so hestenes assumed convergence to avoid numerical precision errors\n", min(y_new))
		}
		update_error = safe_divide(edittozerotol(y_new, eps_threshold ),
		                           editmissing(y_old, epsilon(1)),
		                           epsilon(1) )
		update_error = sqrt(max(update_error))
	}
	else {
		exit(error(100))
	}

	assert_msg(!missing(update_error), "update error is missing")

	S.converged = S.converged + (update_error <= S.tolerance)
	is_last_iter = iter==S.maxiter
	
	if (S.converged >= S.min_ok) {
		S.iteration_count = max((iter, S.iteration_count))
		S.accuracy = max((S.accuracy, update_error))
		if (S.verbose==1) printf("{txt}   converged in %g iterations (last error =%3.1e)\n", iter, update_error)
		if (S.verbose>1) printf("\n{txt}   - Converged in %g iterations (last error =%3.1e)\n", iter, update_error)
	}
	else if (is_last_iter & S.abort) {
		printf("\n{err}convergence not achieved in %g iterations (last error=%e); try increasing maxiter() or decreasing tol().\n", S.maxiter, update_error)
		exit(430)
	}
	else {
		if ((S.verbose>=2 & S.verbose<=3 & mod(iter,1)==0) | (S.verbose==1 & mod(iter,1)==0)) {
			printf("{res}.{txt}")
			displayflush()
		}
		if ( (S.verbose>=2 & S.verbose<=3 & mod(iter,100)==0) | (S.verbose==1 & mod(iter,100)==0) ) {
			printf("{txt}%9.1f\n      ", update_error/S.tolerance)
		}

		if (S.verbose==4 & method!="hestenes") printf("{txt} iter={res}%4.0f{txt}\tupdate_error={res}%-9.6e\n", iter, update_error)
		if (S.verbose==4 & method=="hestenes") printf("{txt} iter={res}%4.0f{txt}\tupdate_error={res}%-9.6e  {txt}norm(ssr)={res}%g\n", iter, update_error, norm(y_new))
		
		if (S.verbose>=5) {
			printf("\n{txt} iter={res}%4.0f{txt}\tupdate_error={res}%-9.6e{txt}\tmethod={res}%s\n", iter, update_error, method)
			"old:"
			y_old
			"new:"
			y_new
		}
	}
	return(S.converged >= S.min_ok)
}

// --------------------------------------------------------------------------

`Matrix' weighted_quadcolsum(`FixedEffects' S, `Matrix' x, `Matrix' y) {
	// BUGBUG: override S.has_weights with pruning
	// One approach is faster for thin matrices
	// We are using cross instead of quadcross but it should not matter for this use
	if (S.has_weights) {
		if (cols(x) < 14) {
			return(quadcross(x :* y, S.weight)')
		}
		else {
			return(diagonal(quadcross(x, S.weight, y))')
		}
	}
	else {
		if (cols(x) < 25) {
			return(diagonal(quadcross(x, y))')
		}
		else {
			return(colsum(x :* y))
		}
	}
}


// RRE benchmarking
//	|| yk - y || / || y || === || ek - e || / || y ||
`Real' reghdfe_rre_benchmark(`Vector' resid, `Vector' true_resid, `Real' norm_y) {
	`Real' ans
	ans = norm(resid - true_resid) / norm_y
	return(ans)
}

end

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
