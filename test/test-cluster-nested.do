cscript "reghdfe fixed effect nested in cluster" adofile reghdfe

* Setup
	discard
	clear all
	set more off
	* cls

* Convenience: "Trim <size>" will trim e(b) and e(V)
	capture program drop TrimMatrix
	program define TrimMatrix, eclass
	args size adjdof
		assert `size'>0
		matrix trim_b = e(b)
		matrix trim_V = e(V)
		matrix trim_b = trim_b[1, 1..`size']
		matrix trim_V = trim_V[1..`size',1..`size']
		if ("`adjdof'"!="") {
			matrix trim_V = trim_V * `adjdof'
			ereturn scalar F = e(F) / `adjdof'
		}
		ereturn matrix trim_b = trim_b
		ereturn matrix trim_V = trim_V
	end
	
* Create fake dataset
	sysuse auto
	bys turn: gen t = _n

* [TEST]
	local lhs price
	local rhs weight length
	local absvars turn // make it the same as panelvar
	local clustervar `absvars'
	fvunab tmp : `rhs'
	local K : list sizeof tmp

	tsset `absvars' t

	* Custom adjustment: in this simple case we can compare _cons
	drop if missing(rep)

// -------------------------------------------------------------------------------------------------

	* Check that _xtreg_chk_cl2 is working

	reghdfe `lhs' `rhs', abs(`absvars') vce(cluster `clustervar') dof(clusters) keepsingletons
	assert e(df_a)==0
	reghdfe `lhs' `rhs', abs(`absvars') vce(cluster `clustervar') keepsingletons
	assert e(df_a)==0

// -------------------------------------------------------------------------------------------------

	* Match against xtreg

	* 1. Run benchmark
	xtreg `lhs' `rhs', cluster(`clustervar') fe
	TrimMatrix `K'
	di e(df_a)
	di e(df_m)
	local bench_df_a = e(df_a)
	local bench_within = e(r2_w)
	storedresults save benchmark e()
	
	* 2. Run reghdfe
	reghdfe `lhs' `rhs', absorb(`absvars') vce(cluster `clustervar') keepsingletons
	di e(df_a)
	di e(df_m)
	* NOTE: See statalist post for discussion on df_m discrepancy
	* http://www.stata.com/statalist/archive/2010-03/msg00941.html	
	* "So I think some explanation is necessary. I see no reason, conceptually, why xtreg,fe with small-sample statistics should not be exactly equivalent to areg"
	local adjdof = e(unclustered_df_r) / (e(unclustered_df_r)-1)
	TrimMatrix `K' `adjdof'

	* 3. Compare
	storedresults compare benchmark e(), tol(1e-12) include( ///
		scalar: N tss rss F df_r ///  
		matrix: trim_b trim_V ///
		macros: wexp wtype )
	
	* WILL NOT HOLD B/C xtreg still include the nested ones in e(df_a)
	* assert `bench_df_a'==e(df_a)-1
	
	assert abs(`bench_within'-e(r2_within))<1e-6
	storedresults drop benchmark

// -------------------------------------------------------------------------------------------------

	* Check that _xtreg_chk_cl2 is working

	reghdfe `lhs' `rhs', abs(`absvars') vce(cluster `clustervar') dof(clusters) keepsingletons
	assert e(df_a)==0
	reghdfe `lhs' `rhs', abs(`absvars') vce(cluster `clustervar') keepsingletons
	assert e(df_a)==0

// -------------------------------------------------------------------------------------------------
// Now repeat the above but with a different clustervar name
// -------------------------------------------------------------------------------------------------
	gen samecluster = `absvars'
	local clustervar samecluster

	* Match against xtreg

	* 1. Run benchmark
	xtreg `lhs' `rhs', cluster(`clustervar') fe
	TrimMatrix `K'
	di e(df_a)
	di e(df_m)
	local bench_df_a = e(df_a)
	local bench_within = e(r2_w)
	storedresults save benchmark e()
	
	* 2. Run reghdfe
	reghdfe `lhs' `rhs', absorb(`absvars') vce(cluster `clustervar') keepsingletons
	di e(df_a)
	di e(df_m)
	* NOTE: See statalist post for discussion on df_m discrepancy
	* http://www.stata.com/statalist/archive/2010-03/msg00941.html	
	* "So I think some explanation is necessary. I see no reason, conceptually, why xtreg,fe with small-sample statistics should not be exactly equivalent to areg"
	local adjdof = e(unclustered_df_r) / (e(unclustered_df_r)-1)
	TrimMatrix `K' `adjdof'

	* 3. Compare
	storedresults compare benchmark e(), tol(1e-12) include( ///
		scalar: N tss rss F df_r ///  
		matrix: trim_b trim_V ///
		macros: wexp wtype )
	
	* WILL NOT HOLD B/C xtreg still include the nested ones in e(df_a)
	* assert `bench_df_a'==e(df_a)-1
	
	assert abs(`bench_within'-e(r2_within))<1e-6
	storedresults drop benchmark


// -------------------------------------------------------------------------------------------------
cd "D:/Github/reghdfe/test"
exit
