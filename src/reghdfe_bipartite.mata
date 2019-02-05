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
	assert(asarray(F12.extra, "levels_as_keys") == 1)

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
