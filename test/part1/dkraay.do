noi cscript "reghdfe: Driscoll-Kraay standard errors" adofile reghdfe

* Check that ivreghdfe is installed
cap which ivreghdfe
if _rc {
	di as error "ivreghdfe not installed; skipping dkraay tests"
	exit
}

* Dataset
	sysuse auto, clear
	bys turn: gen t = _n
	tsset turn t
	drop if missing(rep78)

	* Results to compare (note: ivreghdfe reports different R2s, so we focus on b, V, N, df_r)
	local included_e ///
		scalar: N df_r ///
		matrix: b V

* ==============================================================================
* [TEST 1] Driscoll-Kraay with bandwidth 2 (1 lag)
* ==============================================================================

	local bw 2
	
	local lhs price
	local rhs weight length
	local absvars turn
	local K : list sizeof rhs

	local included ///
		scalar: N rss /// rmse tss F df_r ll ll_0 r2_a df_m 
		matrix: trim_b trim_V ///
		macros: wexp wtype
	
	* 1. Run benchmark: ivreghdfe
	ivreghdfe `lhs' `rhs', absorb(`absvars') dkraay(`bw')
	trim_cons `K'
	mat li e(trim_V)
	storedresults save bench_dk1 e()
	
	* 2. Run reghdfe
	reghdfe `lhs' `rhs', absorb(`absvars') vce(dkraay `bw')
	trim_cons `K' // Need to trim because reghdfe shows a constant
	mat li e(trim_V)
	storedresults compare bench_dk1 e(), tol(1e-12) include(`included') 


* ==============================================================================
* [TEST 2] Driscoll-Kraay with bandwidth 3 (2 lags)
* ==============================================================================

	local bw 3
	
	local lhs price
	local rhs weight length
	local absvars turn
	local K : list sizeof rhs

	local included ///
		scalar: N rss /// rmse tss F df_r ll ll_0 r2_a df_m 
		matrix: trim_b trim_V ///
		macros: wexp wtype
	
	* 1. Run benchmark: ivreghdfe
	ivreghdfe `lhs' `rhs', absorb(`absvars') dkraay(`bw')
	trim_cons `K'
	mat li e(trim_V)
	storedresults save bench_dk2 e()
	
	* 2. Run reghdfe
	reghdfe `lhs' `rhs', absorb(`absvars') vce(dkraay `bw')
	trim_cons `K' // Need to trim because reghdfe shows a constant
	mat li e(trim_V)
	storedresults compare bench_dk2 e(), tol(1e-12) include(`included') 


* ==============================================================================
* [TEST 3] Driscoll-Kraay with bandwidth 4 (3 lags)
* ==============================================================================

	local bw 4
	
	local lhs price
	local rhs weight length
	local absvars turn
	local K : list sizeof rhs

	local included ///
		scalar: N rss /// rmse tss F df_r ll ll_0 r2_a df_m 
		matrix: trim_b trim_V ///
		macros: wexp wtype
	
	* 1. Run benchmark: ivreghdfe
	ivreghdfe `lhs' `rhs', absorb(`absvars') dkraay(`bw')
	trim_cons `K'
	mat li e(trim_V)
	storedresults save bench_dk3 e()
	
	* 2. Run reghdfe
	reghdfe `lhs' `rhs', absorb(`absvars') vce(dkraay `bw')
	trim_cons `K' // Need to trim because reghdfe shows a constant
	mat li e(trim_V)
	storedresults compare bench_dk3 e(), tol(1e-12) include(`included') 


* ==============================================================================
* [TEST 4a] Driscoll-Kraay with bandwidth 1 (0 lags) vs cluster(time)
* ==============================================================================


	local bw 1
	
	local lhs price
	local rhs weight length
	local absvars turn
	local K : list sizeof rhs

	local included ///
		scalar: N rss /// rmse tss F df_r ll ll_0 r2_a df_m 
		matrix: trim_b trim_V ///
		macros: wexp wtype
	
	* 1. Run benchmark: ivreghdfe
	ivreghdfe `lhs' `rhs', absorb(`absvars') dkraay(`bw')
	trim_cons `K'
	mat li e(trim_V)
	storedresults save bench_dk0a e()
	
	* 2. Run reghdfe
	reghdfe `lhs' `rhs', absorb(`absvars') vce(dkraay `bw')
	trim_cons `K' // Need to trim because reghdfe shows a constant
	mat li e(trim_V)
	storedresults compare bench_dk0a e(), tol(1e-12) include(`included') 


* ==============================================================================
* [TEST 4b] Driscoll-Kraay with bandwidth 1 (0 lags) vs cluster(time)
* ==============================================================================

	local bw 1
	
	local lhs price
	local rhs weight length
	local absvars turn
	local K : list sizeof rhs

	local included ///
		scalar: N rss /// rmse tss F df_r ll ll_0 r2_a df_m 
		matrix: trim_b trim_V ///
		macros: wexp wtype
	
	* 1. Run benchmark: ivreghdfe
	ivreghdfe `lhs' `rhs', absorb(`absvars') dkraay(`bw')
	trim_cons `K'
	mat li e(trim_V)
	storedresults save bench_dkmw e()
	
	* 2. Run reghdfe
	reghdfe `lhs' `rhs', absorb(`absvars') vce(dkraay `bw')
	trim_cons `K' // Need to trim because reghdfe shows a constant
	mat li e(trim_V)
	storedresults compare bench_dkmw e(), tol(1e-12) include(`included') 



* ==============================================================================
* [TEST 5] Two-way fixed effects with Driscoll-Kraay
* ==============================================================================

	local absvars turn t
	local bw 2
	
	local lhs price
	local rhs weight length
	local K : list sizeof rhs

	local included ///
		scalar: N rss /// rmse tss F df_r ll ll_0 r2_a df_m 
		matrix: trim_b trim_V ///
		macros: wexp wtype
	
	* 1. Run benchmark: ivreghdfe
	ivreghdfe `lhs' `rhs', absorb(`absvars') dkraay(`bw')
	trim_cons `K'
	mat li e(trim_V)
	storedresults save bench_dk0a e()
	
	* 2. Run reghdfe
	reghdfe `lhs' `rhs', absorb(`absvars') vce(dkraay `bw')
	trim_cons `K' // Need to trim because reghdfe shows a constant
	mat li e(trim_V)
	storedresults compare bench_dk0a e(), tol(1e-12) include(`included') 



* ==============================================================================
* [TEST 6] Different dataset structure - more time periods
* ==============================================================================

	* Create a dataset with more time variation
	clear all
	sysuse auto, clear
	expand 3
	bys turn (foreign price): gen t = _n
	tsset turn t
	drop if missing(rep78)
	
	* Add some noise to make it a proper panel
	replace price = price + rnormal() * 100
	replace weight = weight + rnormal() * 10
	
	local lhs price
	local rhs weight length
	local absvars turn
	local bw 3
	local K : list sizeof rhs

	* 1. Run benchmark: ivreghdfe
	ivreghdfe `lhs' `rhs', absorb(`absvars') dkraay(`bw')
	trim_cons `K'
	mat li e(trim_V)
	storedresults save bench_dk10 e()
	
	* 2. Run reghdfe
	reghdfe `lhs' `rhs', absorb(`absvars') vce(dkraay `bw')
	trim_cons `K' // Need to trim because reghdfe shows a constant
	mat li e(trim_V)
	storedresults compare bench_dk10 e(), tol(1e-12) include(`included') 

exit
