// --------------------------------------------------------------------------
// Regression options
// --------------------------------------------------------------------------
mata:

class FE_Options
{
	`String'		original_varlist	// y x1 x2 (x3 x4 = z1 z2 z3)
	`String'		varlist				// y x1 x2 x3 x4 z1 z2 z3
	
	`String'		original_depvar
	`String'		original_indepvars
	`String'		original_endogvars
	`String'		original_instruments

	`String'		depvar				// y
	`String'		indepvars			// x1 x2
	`String'		endogvars			// x3 x4
	`String'		instruments			// z1 z2 z3
	
	`Boolean'		drop_singletons
	`String'		weight_var			// Weighting variable
	`String'		weight_type			// Weight type (pw, fw, etc)
	`String'		absorb				// contents of absorb()
	`String'		select_if			// If condition
	`String'		select_in			// In condition

	// Parsed input
	`String'		model				// ols, iv
	`String'		summarize_stats
	`Boolean'		summarize_quietly

	`String'		vcetype
	`Integer'		num_clusters
	`Varlist'		clustervars
	`Varlist'		base_clustervars
	`String'		vceextra

	`StringRowVector' dofadjustments // firstpair pairwise cluster continuous
	`Varname'		groupvar
	`String'		residuals
	`String'		diopts

}
end
