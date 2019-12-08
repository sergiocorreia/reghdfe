// --------------------------------------------------------------------------
// Mata Code: FE Estimator (REGHDFE)
// --------------------------------------------------------------------------
// - Project URL:   https://github.com/sergiocorreia/reghdfe
// - Dependency:    https://github.com/sergiocorreia/ftools

    *mata: mata clear
    *mata: mata set matastrict on
    mata: mata set mataoptimize on
    *mata: mata set matadebug off
    *mata: mata set matalnum off

// Include ftools -----------------------------------------------------------
    cap findfile "ftools.mata"
    if (_rc) {
        di as error "reghdfe requires the {bf:ftools} package, which is not installed"
        di as error `"    - install from {stata ssc install ftools:SSC}"'
        di as error `"    - install from {stata `"net install ftools, from("https://github.com/sergiocorreia/ftools/raw/master/src/")"':Github}"'
        exit 9
    }
    include "`r(fn)'"


// Custom types -------------------------------------------------------------
    loc FixedEffects        class FixedEffects scalar
    loc Factors             class Factor rowvector
    loc BipartiteGraph      class BipartiteGraph scalar
    loc FactorPointer       pointer(`Factor') scalar


// Versioning ---------------------------------------------------------------
	ms_get_version reghdfe // from parsetools package
	assert("`package_version'" != "")
    mata: string scalar reghdfe_version() return("`package_version'")
    mata: string scalar reghdfe_stata_version() return("`c(stata_version)'")
    mata: string scalar reghdfe_joint_version() return("`package_version'|`c(stata_version)'")


// Includes -----------------------------------------------------------------
    include "reghdfe_bipartite.mata", adopath
    include "reghdfe_class.mata", adopath
    include "reghdfe_constructor.mata", adopath
    include "reghdfe_common.mata", adopath
    include "reghdfe_projections.mata", adopath
    include "reghdfe_transforms.mata", adopath
    include "reghdfe_accelerations.mata", adopath
    include "reghdfe_lsmr.mata", adopath
