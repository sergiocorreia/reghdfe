// --------------------------------------------------------------------------
// Mata Code: FE Estimator (REGHDFE)
// --------------------------------------------------------------------------
// - Project URL:   https://github.com/sergiocorreia/reghdfe
// - Dependency:    https://github.com/sergiocorreia/ftools


// Miscellanea --------------------------------------------------------------
    findfile "ftools_type_aliases.mata"
    include "`r(fn)'"
    local FixedEffects		class FixedEffects scalar
    local Factors           class Factor rowvector
    local Options           class FE_Options scalar
    local Output			class FE_Output scalar
    loc BipartiteGraph      class BipartiteGraph scalar
    loc FactorPointer       pointer(`Factor') scalar


    mata: mata clear
    mata: mata set matastrict on
    mata: mata set mataoptimize on
    mata: mata set matadebug off
    mata: mata set matalnum off

// Versioning ---------------------------------------------------------------
	ms_get_version reghdfe // from parsetools package
	assert("`package_version'" != "")
    mata: string scalar reghdfe_version() return("`package_version'")
    mata: string scalar reghdfe_stata_version() return("`c(stata_version)'")
    mata: string scalar reghdfe_joint_version() return("`package_version'|`c(stata_version)'")


// Include ------------------------------------------------------------------
    local includes options output bipartite class constructor ///
                   common projections transforms accelerations
    foreach include of local includes {
        findfile "reghdfe_`include'.mata"
        include "`r(fn)'"
    }
