mata:
mata set matastrict on

void function map_precompute(`Problem' S) {
	`Integer' i, g, h, value
	`Varlist' keepvars, cl_ivars
	transmorphic counter, loc
	`Varname' key
	`String' all_clustervars
	if (S.verbose>0) printf("\n{txt}{bf:mata: map_precompute()}\n")

	// Count how many times each var is used, so we can drop them when the counter reaches zero
	counter = asarray_create()
	asarray_notfound(counter, 0)
	// Variables passed through map_init_keepvars()
	keepvars = S.keepvars
	// ivars and cvars of the FEs
	for (g=1; g<=S.G; g++) {
		keepvars = keepvars, S.fes[g].ivars , S.fes[g].cvars
	}
	// Weight var
	if (S.weightvar!="") keepvars = keepvars, S.weightvar
	// Cluster vars
	for (h=1; h<=S.C; h++) {
		cl_ivars = tokens(S.clustervars[h], "#")
		cl_ivars = select(cl_ivars, cl_ivars:!="#")
		keepvars = keepvars, cl_ivars
	}
	// Time and panel vars
	if (S.vce_is_hac & S.timevar!="") keepvars = keepvars, S.timevar
	if (S.vce_is_hac & S.panelvar!="") keepvars = keepvars, S.panelvar

	// Fill
	for (i=1; i<=length(keepvars); i++) {
		asarray(counter, keepvars[i], asarray(counter, keepvars[i])+1)
	}

	// Report
	if (S.verbose>3) printf("{txt}{bf: 0. Usage count of each variable (dropped if it reaches zero)}\n")
	for (i=1; i<=asarray_elements(counter); i++) {
		if (S.verbose>3) printf("{txt}    - key=%s {col 30}count=%f\n", keepvars[i], asarray(counter,keepvars[i]))
	}

	// 1. Store permutation vectors and their invorder, generate ID variables, drop singletons
	if (S.verbose>0) printf("{txt}{bf: 1. Storing permutation vectors, generating ids, dropping singletons}\n")
	if (S.timeit) timer_on(21)
	map_precompute_part1(S, counter)
	if (S.timeit) {
		timer_off(21)
		printf("{res}{col 20}%6.3f{txt}{col 30}precompute 1 (sorts, ids, drop singletons)\n", timer_value(21)[1])
		timer_clear(21)
	} 

	// 2. Store group offsets, group counters; demeaned(x), inv(xx) if num_slopes>0; weightvars
	if (S.verbose>0) printf("{txt}{bf: 2. Storing counters and offsets; processing cvars}\n")
	if (S.timeit) timer_on(22)
	map_precompute_part2(S, counter)
	if (S.timeit) {
		timer_off(22)
		printf("{res}{col 20}%6.3f{txt}{col 30}precompute 2 (counters, offsets, cvars)\n", timer_value(22)[1])
		timer_clear(22)
	} 

	// 3. Create cluster IDs, report whether is/in/nested wrt cluster; store precomputed inv(p)
	if (S.verbose>0) printf("{txt}{bf: 3. Storing reverse permutation vectors, creating cluster IDs}\n")
	if (S.timeit) timer_on(23)
	map_precompute_part3(S, counter)
	if (S.timeit) {
		timer_off(23)
		printf("{res}{col 20}%6.3f{txt}{col 30}precompute 3 (reverse permutations, cluster ids)\n", timer_value(23)[1])
		timer_clear(23)
	} 

	// 4. Keep only the essential variables
	keepvars  = J(1,0,"")
	for (loc=asarray_first(counter); loc!=NULL; loc=asarray_next(counter, loc)) {
		key = asarray_key(counter, loc)
		value = asarray_contents(counter, loc)
		if (value>0) keepvars = keepvars, key
	}
	if (S.verbose>3) printf("{txt}{bf: 4. Keeping the following variables}")
	st_keepvar(keepvars)
	// if (S.verbose>3) printf("\n{txt}    %s\n", invtokens(keepvars))
	if (S.verbose>3) stata("describe _all, numbers") // This resets r(..)!
	if (S.verbose>3) printf("\n")

	// Store N (todo: add other details) to ensure the dataset doesn't change from now on
	S.N = st_nobs()

	// Store updated clustervars (run after describe!)
	for (h=1; h<=S.C;h++) {
		all_clustervars = all_clustervars, S.clustervars[h]
	}
	st_rclear()
	st_global("r(updated_clustervars)", invtokens(all_clustervars))
}
end
