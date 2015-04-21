mata:
mata set matastrict on
// -------------------------------------------------------------------------------------------------
// Structure of the MAP Problem
// -------------------------------------------------------------------------------------------------

struct MapProblem {
	struct FixedEffect vector fes	// The G*1 vector of FE structures
	`Integer'		G 				// Number of FEs when bunching slopes
	`Varname'		weightvar 		// Name variable contaning the fw/pw/aw
	`String'		weighttype 		// fw/pw/aw
	`String'		weights 		// "[weighttype=weightvar]"
	`Series'		w 				// Contents of variable contaning the fw/pw/aw
	`Integer'		verbose			// Number of debug messages to show (0=None, 1=A little, 4=A *lot*)			
	`Integer'		N				// Number of obs; after map_precompute() the dataset CANNOT CHANGE!
	`Varlist'		keepvars		// By default we drop cvars and ivars ASAP; this prevents it (useful for clustervars and for timevar+panelvar under HAC errors)

	`Boolean'		will_save_fe	// True if at least one FE will be saved
	`Boolean'		keepsingletons	// If set to 1, it will not drop singletons (do not touch this!)

	`Integer'		C				// Number of cluster variables
	`Varlist'		clustervars
	`Varlist'		clustervars_original // Note: need to apply tokens()
	`Varname'		panelvar
	`Varname'		timevar
	`Boolean'		vce_is_hac

	`Varname' 		by 				// In case we are using reghdfe .. by()
	
	// Optimization parameters	
	`Integer'		groupsize 		// Group variables when demeaning (more is faster but uses more memory)
	`Real'			tolerance
	`Integer'		maxiterations
	`String'		transform		// Kaczmarz Cimmino Symmetric_kaczmarz (k c s)
	`String'		acceleration	// Acceleration method. None/No/Empty is none\
	`Integer'		accel_start		// Iteration where we start to accelerate // set it at 6? 2?3?
	
	// Specific to Aitken's acceleration
	`Integer'		accel_freq		
	`Integer'		stuck_threshold	// Call the improvement "slow" when it's less than e.g. 1%
	`Integer'		bad_loop_threshold	// If acceleration seems stuck X times in a row, pause it
	`Integer'		pause_length	// This is in terms of candidate accelerations, not iterations (i.e. x3)?

	// Temporary
	`Boolean'		storing_betas
	`Varname'		groupvar		// Name of the variable that will hold the mobility group
	`Varname'		grouptype		// Long, double, etc.
	`Varname'		grouplabel
	`Series'		groupseries		// The actual data of the mobility group variable
	`Series'		uid
	`Series'		resid
	`Varname'		residname

	// Temporary storage for DoFs
	`Integer'		dof_M
	`Integer'		dof_M_due_to_nested
	`Integer'		dof_KminusM
	`Integer'		dof_N_hdfe_extended
	`Vector'		doflist_M
	`Vector'		doflist_M_is_exact
	`Vector'		doflist_M_is_nested
	`Vector'		dof_SubGs
}
end
