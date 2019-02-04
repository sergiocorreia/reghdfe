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
    findfile "reghdfe_bipartite.mata"
    include "`r(fn)'"

    findfile "reghdfe_class.mata"
    include "`r(fn)'"

    findfile "reghdfe_constructor.mata"
    include "`r(fn)'"

    findfile "reghdfe_common.mata"
    include "`r(fn)'"

    findfile "reghdfe_projections.mata"
    include "`r(fn)'"

    findfile "reghdfe_transforms.mata"
    include "`r(fn)'"

    findfile "reghdfe_accelerations.mata"
    include "`r(fn)'"

    findfile "reghdfe_lsmr.mata"
    include "`r(fn)'"
