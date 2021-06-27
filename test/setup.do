* ===========================================================================
* Install reghdfe and ftools from local Git folders
* ===========================================================================
	

// --------------------------------------------------------------------------
// Setup
// --------------------------------------------------------------------------

	set trace off
	set varabbrev on
	log close _all
	clear all
	discard
	cls
	pr drop _all


// --------------------------------------------------------------------------
// Uninstall everything; reinstall dependencies
// --------------------------------------------------------------------------

	cap ado uninstall ftools
	cap ado uninstall reghdfe
	mata: st_local("path", pathresolve(pwd(), "../../ftools/src"))
	di as text `"Installing ftools from "`path'"'
	net install ftools, from("`path'")
	*ftools, compile // need to "clear all" after this


// --------------------------------------------------------------------------
// Ensure that the .mata files don't have trivial errors at compile-time
// --------------------------------------------------------------------------

	* TODO: always inspect by hand for the word "unused"
	cd "../current-code"
	do reghdfe.mata
	cd "../test"


// --------------------------------------------------------------------------
// Build installation files
// --------------------------------------------------------------------------
* This does two things:
* 1) Combine all mata files into a single file reghdfe.mata
* 2) Adds two older versions of reghdfe (v3 and v5)

	cd "../current-code"
	cap rm "../src/reghdfe.ado" // otherwise we won't know if the python script fails
	shell "build.py"
	cd "../test"


// --------------------------------------------------------------------------
// Install reghdfe
// --------------------------------------------------------------------------

	* Reinstall
	* Note: "net install" requires hardcoded paths so we do a workaround
	mata: st_local("path", pathresolve(pwd(), "../src"))
	mata: assert(direxists("`path'"))
	net install reghdfe , from("`path'")


// --------------------------------------------------------------------------
// Quick test to see if it works
// --------------------------------------------------------------------------

	set trace off
	sysuse auto, clear
	which reghdfe
	reghdfe price weight
	reghdfe price weight, a(turn)

	sysuse auto, clear
	clonevar www=weight
	reghdfe price i.rep78 weight www, a(turn)

exit
