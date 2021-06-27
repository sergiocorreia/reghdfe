// --------------------------------------------------------------------------
// Mata Code: regressions with individual FEs and team-level outcomes
// --------------------------------------------------------------------------
// Project URL: ...


// Miscellanea --------------------------------------------------------------
	mata: mata set matastrict on
	mata: mata set mataoptimize on


// Include ftools -----------------------------------------------------------
    cap findfile "ftools.mata"
    if (_rc) {
        di as error "reghdfe requires the {bf:ftools} package, which is not installed"
        di as error `"    - install from {stata ssc install ftools:SSC}"'
        di as error `"    - install from {stata `"net install ftools, from("https://github.com/sergiocorreia/ftools/raw/master/src/")"':Github}"'
        exit 9
    }
    include "`r(fn)'"

// Type aliases -------------------------------------------------------------
	loc Factor				class Factor scalar
	loc OptionalFactor      class Factor
	loc Factors				class Factor rowvector
	loc FactorPointer   	pointer(`Factor') scalar
	loc FactorPointers   	pointer(`Factor') rowvector
	
	loc FE_Factor			class FE_Factor scalar
	loc Optional_FE_Factor	class FE_Factor
	loc FE_Factors			class FE_Factor rowvector
	loc FE_FactorPointer   	pointer(`FE_Factor') scalar
	loc FE_FactorPointers   pointer(`FE_Factor') rowvector
	
	loc Individual_Factor	class Individual_Factor scalar
	//loc Optional_FE_Factor	class FE_Factor
	//loc FE_Factors			class FE_Factor rowvector
	//loc FE_FactorPointer   	pointer(`FE_Factor') scalar
	//loc FE_FactorPointers   pointer(`FE_Factor') rowvector
	
	loc FixedEffects 		class FixedEffects scalar
	loc BipartiteGraph  	class BipartiteGraph scalar
	loc Solution			class Solution scalar

	loc DEBUG				"" // Set to empty string to run slow debugging code, set to "if (0)" for faster runtime

// Includes -----------------------------------------------------------------
	* Notes on include order:
	* - Include "common.mata" first
	* - Include class definitions ("FE.mata") before constructors ("FE_constructor.mata")
	* - Include class definitions before functions that use it ("FE.mata" before "Regression.mata", "DoF.mata", "LSMR.mata", etc)
	//loc includes Mock_Matrix Common Solution FE FE_Constructor Regression Bipartite DoF LSMR LSQR
	

	include "Factor_FE.mata", adopath	// extends Factor() from ftools
	include "Bipartite.mata", adopath		// depends on ftools.mata and ExtendedFactor
	include "Factor_Indiv.mata", adopath	// extends Factor() from ftools; depends on Bipartite
	include "Common.mata", adopath
	include "Solution.mata", adopath
	include "FE.mata", adopath
	include "Regression.mata", adopath
	include "DoF.mata", adopath
	
	include "LSMR.mata", adopath
	include "LSQR.mata", adopath

	include "MAP.mata", adopath
	include "MAP_Accelerations.mata", adopath
	include "MAP_Transformations.mata", adopath
	
	include "Parallel.mata", adopath

exit
