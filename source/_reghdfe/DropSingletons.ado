// -------------------------------------------------------------
// Iteratively drop singletons for each absvar
// -------------------------------------------------------------
* This could be done iteratively, dropping singletons for each absvar until no progress is made.
* However, that would be extremely slow for a modest gain

cap pr drop DropSingletons
program define DropSingletons, sortpreserve
syntax, num_absvars(integer)

	forv g=1/`num_absvars' {
		reghdfe_absorb, fe2local(`g')
		* ivars cvars target varname varlabel is_interaction is_cont_interaction is_bivariate is_mock levels
		
		* It's either redudant (the second part of i.a##c.b) or tricky (a simple i.a#c.b) to discard singletons with cont. interactions
		local is_slope =  ("`cvars'"!="") & (!`is_bivariate' | `is_mock')
		if (`is_slope') continue

		local N_old = c(N)
		qui bys `ivars': drop if _N==1
		local N_new = c(N)
		local N_dropped = (`N_old' - `N_new')
		if (`N_dropped'>0) Debug, level(0) msg("(dropped `N_dropped' singleton observations for absvar " as result "`varlabel'" as text ")")
	}
end
