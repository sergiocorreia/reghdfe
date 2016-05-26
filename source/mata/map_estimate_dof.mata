mata:
mata set matastrict on
void map_estimate_dof(`Problem' S, string rowvector adjustments, 
		| `Varname' groupvar, `String' cond) {
	`Boolean' adj_two, adj_pairwise, adj_clusters, adj_continuous, belongs, already_first_constant
	string rowvector all_adjustments
	`String' adj, label, basestring
	`Integer' i, g, SuperG, h, M_due_to_nested, j, m, sum_levels, m2
	`Vector' M, M_is_exact, M_is_nested, is_slope, solved, i_to_g, SubGs

	// Parse list of adjustments/tricks to do
	if (S.verbose>1) printf("\n")
	if (S.verbose>0) printf("{txt}{bf:mata: map_estimate_dof()}\n")
	if (S.verbose>1) printf("{txt} - Estimating degrees-of-freedom used by the fixed effects\n")
	all_adjustments = "two", "three", "pairwise", "clusters", "continuous"
	adjustments = tokens(adjustments)
	for (i=1; i<=length(adjustments);i++) {
		assert_msg(anyof(all_adjustments, adjustments[i]), 
			sprintf("map_estimate_dof error: adjustment %s invalid", adjustments[i]))
	}
	if (S.verbose>1) printf("{txt} - Adjustments:\n")
	for (i=1;i<=length(all_adjustments);i++) {
		adj = all_adjustments[i]
		belongs = anyof(adjustments, adj)
		if (S.verbose>1) printf("{txt}    - %s:  {col 20}{res} %s\n", adj, belongs ? "yes" : "no")
		if (adj=="two") adj_two = belongs
		if (adj=="three") adj_three = belongs
		if (adj=="pairwise") adj_pairwise = belongs
		if (adj=="clusters") adj_clusters = belongs
		if (adj=="continuous") adj_continuous = belongs
	}

	// Assert that the clustervars exist
	for (i=1;i<=S.C;i++) {
		stata(sprintf("confirm numeric variable %s, exact", S.clustervars[i]))
	}

	// Can only save connected group if two or pairwise are active
	if (args()<3) groupvar = ""
	if (groupvar!="") {
		assert_msg(adj_two | adj_pairwise, "map_estimate_dof error: group option requires 'pairwise' or 'two' adjustments")
	}

	// Count all fixed intercepts and slopes
	SubGs = J(S.G, 1, 0) // Intercept + # of slopes in an absvar
	for (g=1;g<=S.G;g++) {
		SubGs[g] = S.fes[g].has_intercept + S.fes[g].num_slopes
	}
	SuperG = sum(SubGs)
	if (S.verbose>1) printf("{txt} - There are %f fixed intercepts and slopes in the %f absvars\n", SuperG, S.G)

	// Initialize result vectors and scalars
	M = J(SuperG, 1, 1)
	M_is_exact = J(SuperG, 1, 0)
	M_is_nested = J(SuperG, 1, 0)
	is_slope = J(SuperG, 1, .)
	solved = J(SuperG, 1, 0)

	// Initial Fill
	h = 0
	already_first_constant = 0
	for (g=1;g<=S.G;g++) {
		for (i=1;i<=SubGs[g];i++) {
			h++
			if (is_slope[h] = i>S.fes[g].has_intercept) {
				M[h] = 0
			}
			else if (!already_first_constant) {
				already_first_constant = 1
				M[h] = 0
			}

		}
	}

	// (Intercept-Only) Look for absvars that are clustervars or are nested within a clustervar
	h = 1
	M_due_to_nested = 0
	if (adj_clusters) {
		for (g=1;g<=S.G;g++) {
			if (S.fes[g].has_intercept & (S.fes[g].is_clustervar | S.fes[g].in_clustervar)) {
				M[h] = S.fes[g].levels
				M_is_exact[h] = M_is_nested[h] = 1
				M_due_to_nested = M_due_to_nested + M[h]
				solved[h] = 1
				if (S.verbose>1 & S.fes[g].is_clustervar) printf("   {txt}(categorical variable {res}%s{txt} is also a cluster variable, so it doesn't count towards DoF)\n", invtokens(S.fes[g].ivars,"#"))
				if (S.verbose>1 & S.fes[g].in_clustervar) printf("   {txt}(categorical variable {res}%s{txt} is nested within cluster {res}%s{txt}, so it doesn't count towards DoF)\n", invtokens(S.fes[g].ivars,"#"), S.clustervars_original[S.fes[g].nesting_clustervar])
			}
			h = h + SubGs[g]
		}
	}

	// (Intercept-only) Excluding those already solved, the first absvar is exact, and the second can be with pairwise/two

	// Note: I dont't include the FEs that are clusters or are nested within cluster when computing redundant coefs
	// On principle, that would be nice to have. EG: reghdfe .... abs(zipcode state##c.time) vce(zipcode)
	// I know state is collinear with zipcode so I would also want to consider state to be redundant

	// However, the number of states should be much smaller than the number of zipcodes, which in turn is smaller
	// Than the number of observations; so I don't worry much about that case (also, there may be possible 
	// complications with that)

	i = 0 // Different from -g- as it iterates over effective FEs (i<=g)
	h = 1 // Different from -g- as it counts multiple FEs in one g (h>=g)
	m2 = 0 // needed for group3hdfe
	i_to_g = J(S.G, 1, 0)
	for (g=1;g<=S.G;g++) {

		assert( (i <= g) & (g <= h) )

		if (!solved[h] & S.fes[g].has_intercept) {
			i++
			if (i==1) {
				M_is_exact[h] = 1
			}
			else if (i==2 & (adj_pairwise | adj_two | adj_three)) {
				M_is_exact[h] = 1
				m = map_connected_groups(S, i_to_g[1], g, groupvar)
				if (S.verbose>2) printf("{txt}    - Mobility groups between fixed intercept #%f and #%f: {res}%f\n", i_to_g[1], g, m)
				M[h] = m
				m2 = m
			}
			else if (i>2 & adj_pairwise) {
				// Call connected in a LOOP (but I need to save the list of those that I needed to check)
				for (j=1; j<i; j++) {
					m = map_connected_groups(S, i_to_g[j], g)
					if (S.verbose>2) printf("{txt}    - Mobility groups between fixed intercept #%f and #%f: {res}%f\n", i_to_g[j], g, m)
					M[h] = max((M[h], m))
				}
				if (S.verbose>2) printf("{txt}    - Maximum of mobility groups wrt fixed intercept #%f: {res}%f\n", g, M[h])

			}
			else if (i==3 & adj_three) {
				M_is_exact[h] = 1
				cmd = "qui group3hdfe __ID%s__ __ID%s__ __ID%s__"
				cmd = sprintf(cmd, strofreal(i_to_g[i-2]), strofreal(i_to_g[i-1]), strofreal(g))
				stata(cmd)
				m = st_numscalar("r(rest)") - m2
				if (S.verbose>2) printf("{txt}    - Redundant coefs in #%f from previous FEs: {res}%f\n", g, m)
				M[h] = m
			}
			i_to_g[i] = g
		}
		h = h + SubGs[g]
	}

	// See if cvars are zero (w/out intercept) or just constant (w/intercept)
	if (adj_continuous) {
		h = 0
		for (g=1;g<=S.G;g++) {
			for (i=1;i<=SubGs[g];i++) {
				h++
				// If model has intercept, redundant cvars are those that are CONSTANT
				// Without intercept, a cvar has to be zero within a FE for it to be redundant
				// Since S.fes[g].x are already demeaned IF they have intercept, we don't have to worry about the two cases
				if (is_slope[h]) {
					M[h] = count_redundant_cvars(S, g, i)
				}
			}
		}
	}

	// Store results
	S.dof_SubGs = SubGs
	
	S.doflist_M = M
	S.doflist_M_is_exact = M_is_exact
	S.doflist_M_is_nested = M_is_nested

	sum_levels = 0
	for (g=1;g<=S.G;g++) sum_levels = sum_levels + S.fes[g].levels * (S.fes[g].has_intercept + S.fes[g].num_slopes)
	S.dof_M = sum(M)
	S.dof_KminusM = sum_levels - S.dof_M
	S.dof_M_due_to_nested = M_due_to_nested
	S.dof_N_hdfe_extended = SuperG

	st_numscalar("e(df_a)", S.dof_KminusM) // We need this in the regression stage!

	// Report results
	if (S.verbose>=2) {
		printf("{txt} - Degrees-of-freedom used by each fixed effect (K=total levels; M=redundant levels)\n")
		h = 0
		for (g=1;g<=S.G;g++) {
			for (i=1;i<=SubGs[g];i++) {
				h++
				label = invtokens(S.fes[g].ivars, "#")
				if (i>S.fes[g].has_intercept) label = label + "#c." + S.fes[g].cvars[i-S.fes[g].has_intercept]
				basestring = "{txt}   - FE%f ({res}%s{txt}): {col 40}K=%f {col 50}M=%f {col 60}is_exact=%f\n"
				printf(basestring, g, label, S.fes[g].levels, M[h], M_is_exact[h])
			}
		}
	}
	if (S.verbose>0) printf(" - Results: N=%f ; K=%f ; M=%f ; (K-M)==df_a=%f\n", S.N, sum_levels, sum(M), sum_levels-sum(M))
}
// -------------------------------------------------------------------------------------------------

void function map_ereturn_dof(`Problem' S) {
	`Integer' h, g, i

	st_numscalar("e(N_hdfe)", S.G)
	st_numscalar("e(N_hdfe_extended)", S.dof_N_hdfe_extended)
	st_numscalar("e(redundant)", S.dof_M)
	st_numscalar("e(M_due_to_nested)", S.dof_M_due_to_nested)
	st_numscalar("e(df_a)", S.dof_KminusM)

	h = 0
	for (g=1;g<=S.G;g++) {
		for (i=1;i<=S.dof_SubGs[g];i++) {
			h++
			st_numscalar( sprintf("e(M%f)",h) , S.doflist_M[h] )
			st_numscalar( sprintf("e(K%f)",h) , S.fes[g].levels )

			st_numscalar( sprintf("e(M%f_exact)",h) , S.doflist_M_is_exact[h])
			st_numscalar( sprintf("e(M%f_nested)",h) , S.doflist_M_is_nested[h])
			st_numscalar( sprintf("e(G%f)",h) , g) // unused?
		}
	}
}

// -------------------------------------------------------------------------------------------------

`Integer' function count_redundant_cvars(`Problem' S, `Integer' g, `Integer' i) {
	`Integer' j, i_lower, i_upper, ans, L, ii
	real rowvector min_max
	`Series' x

	ii = i-S.fes[g].has_intercept
	ans = 0
	L = S.fes[g].levels
	x = S.fes[g].x[., ii]
	
	i_lower = 1
	for (j=1;j<=L; j++) {
		i_upper = S.fes[g].offsets[j]
		min_max = minmax(x[| i_lower \ i_upper |])
		if (sqrt(epsilon(1))>abs(min_max[1]) & sqrt(epsilon(1))>abs(min_max[2])) ans++
		i_lower = i_upper + 1
	}
	if (S.verbose>=2) printf("{txt}    - Fixed slope {res}%s#c.%s {txt}has {res}%f/%f{txt} redundant coefs.\n", invtokens(S.fes[g].ivars,"#"), S.fes[g].cvars[ii], ans, L)
	return(ans)
}

end
