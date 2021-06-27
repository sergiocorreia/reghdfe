// --------------------------------------------------------------------------
// Estimate degrees-of-freedom (DoF) lost due to fixed effects
// --------------------------------------------------------------------------

mata:

`Void' estimate_dof(`FixedEffects' S, `StringRowVector' dofadjustments, `Varname' groupvar)
{
	`Integer'               g, n, i, j
	`Boolean'				skip_next
	`RowVector'				intercept_index
	`Integer'				indiv_fe_i, indiv_fe_g

	
	if (S.verbose > 0) printf("\n{txt}# Estimating degrees-of-freedom absorbed by the fixed effects\n\n")

	// Accessible properties:
	// S.factors[g] 
	//	.has_intercept
	//	.num_slopes
	//  .is_individual_fe

	// Here the focus is not on 'G' but on 'G_extended' which counts slopes and intercepts separately
	// EG: "i.turn i.trunk##c.(price gear)" has G=2 but G_extended=4: i.turn i.trunk i.trunk#c.price i.trunk#c.gear

	S.df_a_nested = 0
	
	// Compute G_extended and initialize table with number of coefficients per extended absvar
	S.doflist_K = J(1, 0, .)
	S.G_extended = 0
	for (g=1; g<=S.G; g++) {
		n = S.factors[g].has_intercept + S.factors[g].num_slopes
		S.G_extended = S.G_extended + n
		S.doflist_K = S.doflist_K , J(1, n, S.factors[g].num_levels)
	}

	// Initialize table with number of redundant coefficients for each extended absvar
	S.doflist_M = S.doflist_M_is_exact = S.doflist_M_is_nested = J(1, S.G_extended, 0)
	
	// intercept index: for each extended absvar, either 0 or 'g' (when the 'g' has intercept)
	intercept_index = J(1, S.G_extended, 0)
	for (g=i=1; g<=S.G; g++) {
		if (S.factors[g].has_intercept & !S.factors[g].is_individual_fe) intercept_index[i] = g
		if (S.factors[g].is_individual_fe) {
			assert(indiv_fe_i == .) // There can only be ONE indiv FE
			indiv_fe_i = i
			indiv_fe_g = g
		}
		n = S.factors[g].has_intercept + S.factors[g].num_slopes
		i = i + n
	}

	// Look for absvars that are clustervars or are nested within clustervars
	assert(0 <= S.num_clusters & S.num_clusters <= 10)
	if (S.num_clusters > 0 & anyof(dofadjustments, "clusters")) {
		dof_update_nested(S, intercept_index)
	}

	// (Intercept-Only) Every intercept but the first has at least one redundant coef.
	// Note that this excludes the ones nested within clusters
	// If we DO have FEs nested within clusters, we should also include the first intercept
	n = cols(selectindex(intercept_index))
	if ((n > 1) | (n > 0 & S.df_a_nested > 0)) {
		if (S.verbose > 0) printf("{txt}   - there is at least one redundant coef. for every set of FE intercepts after the first one\n")
		skip_next = S.df_a_nested == 0
		for (i=1; i<=S.G_extended; i++) {
			g = intercept_index[i]
			if (!g) continue
			if (skip_next) {
				skip_next = 0
				continue
			}
			S.doflist_M[i] = 1
		}
	}

	// First intercept absvar always has an exactly computed DoF
	for (i=1; i<=S.G_extended; i++) {
		g = intercept_index[i]
		if (!g) continue
		S.doflist_M_is_exact[i] = 1
		break
	}

	// Mobility group algorithm
	if (anyof(dofadjustments, "firstpair") | anyof(dofadjustments, "pairwise")) {
		// TODO: add group3hdfe
		dof_update_mobility_group(S, intercept_index, dofadjustments, groupvar)
	}

	if (anyof(dofadjustments, "continuous")) {
		dof_update_cvars(S)
	}

	// Algorithm for individual FEs
	if (S.individual_id != "") {
		dof_update_individual_fe(S, indiv_fe_i, indiv_fe_g)
	}

	// Store results (besides doflist_..., etc.)
	S.df_a_initial = sum(S.doflist_K)
	S.df_a_redundant = sum(S.doflist_M)
	S.df_a = S.df_a_initial - S.df_a_redundant
	S.dofadjustments = dofadjustments
	if (S.verbose > 0) printf("\n")
}


`Void' dof_update_nested(`FixedEffects' S,
						 `RowVector' intercept_index)
{
	`Integer'				i, g, i_cluster
	`String'                absvar, clustervar, vars
	`DataCol'               cluster_data
	`Factor'                F

	// (1) (Intercept-Only) Look for absvars that are clustervars
	for (i=1; i<=S.G_extended; i++) {
		g = intercept_index[i]
		if (!g) continue
		absvar = invtokens(S.factors[g].ivars, "#")
		if (anyof(S.clustervars, absvar)) {
			S.doflist_M[i] = S.doflist_K[i]
			S.df_a_nested = S.df_a_nested + S.doflist_M[i]
			S.doflist_M_is_exact[i] = S.doflist_M_is_nested[i] = 1
			intercept_index[i] = 0
			if (S.verbose > 0) printf("{txt}   - categorical variable {res}%s{txt} is also a cluster variable, so it doesn't reduce DoF\n", absvar)
		}
	}

	// (2) (Intercept-Only) Look for absvars that are nested within a clustervar
	for (i_cluster=1; i_cluster<= S.num_clusters; i_cluster++) {
		cluster_data = .
		if (!any(intercept_index)) break // no more absvars to process
		for (i=1; i<=S.G_extended; i++) {
			g = intercept_index[i]
			if (!g | S.doflist_M_is_exact[i]) continue // nothing to do
			absvar = invtokens(S.factors[g].ivars, "#")
			clustervar = S.clustervars[i_cluster]

			// Load clustervar if needed
			if (cluster_data == .) {
				if (strpos(clustervar, "#")) {
					vars = subinstr(clustervar, "#", " ", .)
					F = factor(vars, S.sample, ., "", 0, 0, ., 0)
					cluster_data = F.levels
					F = Factor() // clear
				}
				else {
					cluster_data = __fload_data(clustervar, S.sample, 0)
				}
			}

			if (S.factors[g].nested_within(cluster_data)) {
				S.doflist_M[i] = S.doflist_K[i]
				S.doflist_M_is_exact[i] = S.doflist_M_is_nested[i] = 1
				S.df_a_nested = S.df_a_nested + S.doflist_M[i]
				intercept_index[i] = 0
				if (S.verbose > 0) printf("{txt}   - categorical variable {res}%s{txt} is nested within cluster variable {res}%s{txt}, so it doesn't reduce DoF\n", absvar, clustervar)
			}
		}
	}
	cluster_data = . // save memory
}


`Void' dof_update_mobility_group(`FixedEffects' S,
								 `RowVector' intercept_index,
								 `StringRowVector' dofadjustments,
								 `Varname' groupvar)
{
	`BipartiteGraph'        bg
	`Integer'               bg_verbose                  // verbose level when calling BipartiteGraph()
	`Integer'               pair_count
	`Integer'               i, j
	`Integer'               g, h
	`Boolean'				save_subgraph
	`String'				grouplabel
	`Integer'               m                           // Mobility groups between a specific pair of FEs

	bg = BipartiteGraph()
	bg_verbose = max(( S.verbose - 1 , 0 ))
	pair_count = 0

	for (i=1; i<=S.G_extended; i++) {
		g = intercept_index[i]
		if (!g) continue

		for (j=i+1; j<=S.G_extended; j++) {
			h = intercept_index[j]
			if (!h) continue

			bg.init(&(S.factors[g]), &(S.factors[h]), bg_verbose)
			++pair_count
			save_subgraph = (pair_count == 1) & (groupvar != "")
			m = bg.init_zigzag(save_subgraph)
			if (S.verbose > 0) printf("{txt}   - mobility groups between FE intercepts #%f and #%f: {res}%f\n", g, h, m)
			if (save_subgraph) {
				if (S.verbose > 0) printf("{txt}   - Saving identifier for the first mobility group: {res}%s\n", groupvar)
				st_store(S.sample, st_addvar("long", groupvar), bg.subgraph_id)
				grouplabel = sprintf("Mobility group between %s and %s", invtokens(S.factors[g].ivars, "#"), invtokens(S.factors[h].ivars, "#"))
				st_varlabel(groupvar, grouplabel)
			}
			S.doflist_M[j] = max(( S.doflist_M[j] , m ))
			if (pair_count==1) S.doflist_M_is_exact[j] = 1
			if (pair_count & anyof(dofadjustments, "firstpair")) break
		}
		if (pair_count & anyof(dofadjustments, "firstpair")) break
	}
}


`Void' dof_update_cvars(`FixedEffects' S)
{
	`Integer'               i, j, g, k
	`Boolean'               has_int
	`Matrix'                tmp, sorted_x
	`RowVector'             zeros, results

	for (i=g=1; g<=S.G; g++) {
		// If model has intercept, redundant cvars are those that are CONSTANT
		// Without intercept, a cvar has to be zero within a FE for it to be redundant
		// Since S.fes[g].x are already demeaned IF they have intercept, we don't have to worry about the two cases
		has_int = S.factors[g].has_intercept
		if (has_int) i++
		k = S.factors[g].num_slopes
		if (!k) continue
		if (S.factors[g].is_individual_fe) continue // currently not adjusting DoFs with slopes of individual FEs
		assert(k == cols(S.factors[g].unsorted_x))
		results = J(1, k, 0)
		// float(1.1) - 1 == 2.384e-08 , so let's pick something bigger, 1e-6
		zeros = J(1, k, 1e-6)

		sorted_x = S.factors[g].sort(S.factors[g].unsorted_x)

		// BUGBUG: This part is probably VERY SLOW !!!
		// This can be speed up by moving the -if- outside the -for-
		// Maybe we can just call a function with fcollapse ???
		for (j = 1; j <= S.factors[g].num_levels; j++) {
			tmp = colminmax(panelsubmatrix(sorted_x, j, S.factors[g].info)) // TODO: this line might be SLOW!!!
			if (has_int) {
				results = results + ((tmp[2, .] - tmp[1, .]) :<= zeros)
			}
			else {
				results = results + (colsum(abs(tmp)) :<= zeros)
			}
		}

		if (sum(results)) {
			if (has_int  & (S.verbose > 0)) printf("{txt}   - the slopes in the FE #%f are constant for {res}%f{txt} levels, which don't reduce DoF\n", g, sum(results))
			if (!has_int & (S.verbose > 0)) printf("{txt}   - the slopes in the FE #%f are zero for {res}%f{txt} levels, which don't reduce DoF\n", g, sum(results))
			S.doflist_M[i..i+S.factors[g].num_slopes-1] = results
		}
		i = i + S.factors[g].num_slopes
	}
}

`Void' dof_update_individual_fe(`FixedEffects' S,
								`Integer' i,
								`Integer' g)
{
	assert(S.doflist_M[i] == 0)
	// TODO: add algorithm(s) here
	// For instance, two authors that worked at the SAME papers are redundant
	// There's probably a nice algorithm or data structure that allows us to see which indivs have the same groups
	// More generally, if indiv A always worked with B or C, and B,C only worked with A, then A is redundant
}

end
