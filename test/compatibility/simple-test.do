* ===========================================================================
* Test that the programs that depend on reghdfe *run* (not tested: that they run correctly)
* ===========================================================================

	clear all
	cls
	ms_get_version reghdfe, min_version(6.0.0)

// --------------------------------------------------------------------------
// ivreghdfe
// --------------------------------------------------------------------------

	which ivreghdfe

	sysuse auto, clear
	ivreghdfe price weight (length=gear headroom), a(turn trunk) 
