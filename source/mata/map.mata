// -------------------------------------------------------------------------------------------------
// Mata Code: Method of Alternating Projections with Acceleration
// -------------------------------------------------------------------------------------------------
	// To debug the mata code, uncomment this three lines, and then -do- the file
	//discard
	//pr drop _all
	//clear all

* Type Aliases
	local Boolean 		real scalar
	local Integer 		real scalar
	local Real 			real scalar
	local Vector		real colvector
	local Matrix		real matrix
	local Series		real colvector // Should be N*1
	local Group			real matrix // Should be N*K
	local String 		string scalar
	local Varname 		string scalar
	local Varlist 		string rowvector // rowvector so they match tokens()
	local Problem		struct MapProblem scalar
	local FE			struct FixedEffect scalar
	local FunctionPointer pointer(`Group' function) scalar // Used for the Accelerate & Transform fns

// -------------------------------------------------------------------------------------------------
	include assert_msg.mata
	include FixedEffect.mata
	include MapProblem.mata
	include map_common.mata
	include map_init.mata
	include map_precompute.mata
	include map_precompute_part1.mata
	include map_precompute_part2.mata
	include map_precompute_part3.mata
	include map_projection.mata
	include map_solve.mata
	include map_solve_accelerations.mata
	include map_solve_transformations.mata
	include map_estimate_dof.mata
	include map_connected_groups.mata

	// This is not part of the MAP code but for simplicity we'll put it here
	include fix_psd.mata
// -------------------------------------------------------------------------------------------------
