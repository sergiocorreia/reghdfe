mata:
mata set matastrict on
// -------------------------------------------------------------------------------------------------
// Structure of a Fixed Effect
// -------------------------------------------------------------------------------------------------

	struct FixedEffect {
		`Integer'	order 			// "g", the position in varlist
		`Varname'	varlabel		// Original label of this absvar
		`Integer'	num_slopes
		`Integer'	has_intercept
		`Integer'	levels			// Number of categories spanned by the ivars
		`Varlist'	ivars			// number of i.var elements
		`Varlist'	cvars			// number of c.var elements or slopes
		`Boolean'	is_sortedby		// 1 if the dataset is sorted by the ivars of this FE
		`Varname'	idvarname		// (optional) Name of variable with the absvar categories
		`Varlist'	target			// Name of the variable that will hold the estimates for the FE
		
		`Series'	p 				// Permutation vector
		`Series'	inv_p 			// Precompute invorder(p) because we'll use it a lot

		`Vector'	offsets			// Pos. of last obs. with that category (when sorted)
		`Vector'	counts			// number of obs. (weighted) with that value
		`Group'		x				// Vector/Matrix of (optionally demeaned) cvars
		`Matrix'	inv_xx			// Blocks of the inv(x'x) matrix; size KL*K (k=num_slopes, L=levels)
		`Matrix'	xmeans

		`Boolean'	is_clustervar, in_clustervar
		`Integer'	nesting_clustervar // Clustervar that nests this FE, if any

		// Temporary matrices for the stored FEs
		`Matrix'	alphas
		`Matrix'	tmp_alphas
	}
end
