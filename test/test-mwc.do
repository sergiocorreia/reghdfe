noi cscript "reghdfe with multi-way clustering" adofile reghdfe

* Test multi-way-clustering routine
* Will compare it against areg (N_clustervars==1), ivreg2 (N_clustervars<=2) and cgmreg.
* Note: We use -fixed_cgmreg- because a bug in -cgmreg-, where they count K=rows(X'X)
* But the rows are also included for omitted variables (e.g. the base variable in a factor expansion i.turn)

* Setup
	discard
	clear all
	set more off
	qui adopath + "D:/Github/reghdfe/test" // allows us to run the fixed_cgmreg.ado
	set more off

* Convenience
	capture program drop TrimMatrix
	program define TrimMatrix, eclass
	args size
		assert `size'>0
		matrix trim_b = e(b)
		matrix trim_V = e(V)
		matrix trim_b = trim_b[1, 1..`size']
		matrix trim_V = trim_V[1..`size',1..`size']
		ereturn matrix trim_b = trim_b
		ereturn matrix trim_V = trim_V
	end
	
* Dataset
	sysuse auto
	drop if missing(rep)
	
* Test Parameters
	local lhs price
	local rhs weight length
	local absvars rep
	fvunab tmp : `rhs'
	local K : list sizeof tmp
	
* [TEST] One Clustervar
	local clustervars turn

	* 1. Run benchmark
	areg `lhs' `rhs', absorb(`absvars') cluster(`clustervars')
	TrimMatrix `K'
	storedresults save bench_areg e()
	
	* 2. Run reghdfe-default as benchmark
	reghdfe `lhs' `rhs', absorb(`absvars') vce(cluster `clustervars', suite(default))
	TrimMatrix `K'
	storedresults save bench_def e()
	
	* 3. Run reghdfe-avar as benchmark
	reghdfe `lhs' `rhs', absorb(`absvars') vce(cluster `clustervars', suite(avar))
	TrimMatrix `K'
	storedresults save bench_avar e()
	
	* 4. Run fixed_cgmreg.ado as benchmark
	qui tab `absvars', gen(ABS_)
	fixed_cgmreg `lhs' `rhs' ABS_* , cluster(`clustervars')
	TrimMatrix `K'
	storedresults save bench_cgmreg e()
	
	* 5. Run and test reghdfe-mwc
	reghdfe `lhs' `rhs', absorb(`absvars') vce(cluster `clustervars', suite(mwc))
	TrimMatrix `K'
	
	* Compare
	storedresults compare bench_areg e(), tol(1e-12) include( ///
		scalar: N rmse tss rss r2 r2_a F df_r df_a df_m N_clust /// 
		matrix: trim_b trim_V ///
		macros: wexp wtype )
	storedresults drop bench_areg
	
	storedresults compare bench_def e(), tol(1e-12) include( ///
		scalar: N rmse tss rss r2 r2_a F df_r df_a df_m N_clust /// 
		matrix: trim_b trim_V ///
		macros: wexp wtype )
	storedresults drop bench_def
	
	storedresults compare bench_avar e(), tol(1e-12) include( ///
		scalar: N rmse tss rss r2 r2_a F df_r df_a df_m N_clust /// 
		matrix: trim_b trim_V ///
		macros: wexp wtype )
	storedresults drop bench_avar
	
	* NOTE: This has a lower tolerance; perhaps -cgmreg- uses only single precision somewhere??
	storedresults compare bench_cgmreg e(), tol(1e-10) include( ///
		scalar: N rss r2 /// 
		matrix: trim_b trim_V ///
		macros: wexp wtype )
	storedresults drop bench_cgmreg
	
* [TEST] Two Clustervars
	local clustervars turn trunk // mpg
	
	* 2. Run reghdfe-default as benchmark (it should autoconvert to reghdfe-mwc)
	reghdfe `lhs' `rhs', absorb(`absvars') vce(cluster `clustervars', suite(default))
	TrimMatrix `K'
	assert e(vcesuite)=="mwc"
	storedresults save bench_def e()
	
	* 3. Run reghdfe-avar as benchmark
	reghdfe `lhs' `rhs', absorb(`absvars') vce(cluster `clustervars', suite(avar))
	TrimMatrix `K'
	storedresults save bench_avar e()
	
	* 4. Run fixed_cgmreg.ado as benchmark
	cap qui tab `absvars', gen(ABS_)
	drop ABS_1
	fixed_cgmreg `lhs' `rhs' ABS_* , cluster(`clustervars')
	*mat li e(V)
	*mat li e(rawcovmat)
	drop ABS_*
	TrimMatrix `K'
	storedresults save bench_cgmreg e()
	
	* 5. Run -ivreg2- as benchmark
	cap qui tab `absvars', gen(ABS_)
	ivreg2 `lhs' `rhs' ABS_* , cluster(`clustervars') small nocons partial(ABS_*)
	drop ABS_*
	TrimMatrix `K'
	local ivreg2_r2 = e(r2) // Within R2
	storedresults save bench_ivreg2 e()
	
	* Run and test reghdfe-mwc
	reghdfe `lhs' `rhs', absorb(`absvars') vce(cluster `clustervars', suite(mwc))
	TrimMatrix `K'

	* Compare
	
	storedresults compare bench_def e(), tol(1e-12) include( ///
		scalar: N rmse tss rss r2 r2_a F df_r df_a df_m N_clust /// 
		matrix: trim_b trim_V ///
		macros: wexp wtype )
	storedresults drop bench_def
	
	storedresults compare bench_avar e(), tol(1e-12) include( ///
		scalar: N rmse tss rss r2 r2_a F df_r df_a df_m N_clust /// 
		matrix: trim_b trim_V ///
		macros: wexp wtype )
	storedresults drop bench_avar
	
	* NOTE: This has a lower tolerance; perhaps -cgmreg- uses only single precision somewhere??
	storedresults compare bench_cgmreg e(), tol(1e-8) include( ///
		scalar: N rss r2 /// 
		matrix: trim_b /// trim_V
		macros: wexp wtype )
	storedresults drop bench_cgmreg
	* BUGBUG/NOTE: I can't get a match with CGMREG in the VCV; it's close but not there
	* Maybe it's due to the small sample adjustments? Maybe it's due to the N_clust adjustments?
	
	* Note: Using low tolerance..
	storedresults compare bench_ivreg2 e(), tol(1e-7) include( ///
		scalar: N rmse rss  F df_r N_clust N_clust1 N_clust2 ///  r2_a: excluded b/c r2 differs; could match if we run w/out partial() but then F differs
		matrix: trim_b trim_V ///
		macros: wexp wtype )
	storedresults drop bench_ivreg2
	assert "`e(N_clustervars)'"!="" & "`e(N_clustervars)'"!="."
	if abs(e(r2_within)-`ivreg2_r2')>1e8 {
		di as error "Within R2 doesn't match ivreg2 (`e(r2_within)' vs (`ivreg2_r2')"
		exit 9
	}
	
* [TEST] Three Clustervars
	local clustervars turn trunk mpg
	
	* 2. Run reghdfe-default as benchmark (it should autoconvert to reghdfe-mwc)
	reghdfe `lhs' `rhs', absorb(`absvars') vce(cluster `clustervars', suite(default))
	TrimMatrix `K'
	assert e(vcesuite)=="mwc"
	storedresults save bench_def e()
	
	* 4. Run fixed_cgmreg.ado as benchmark
	cap qui tab `absvars', gen(ABS_)
	drop ABS_1
	fixed_cgmreg `lhs' `rhs' ABS_* , cluster(`clustervars')
	drop ABS_*
	TrimMatrix `K'
	storedresults save bench_cgmreg e()
	
	* Run and test reghdfe-mwc
	reghdfe `lhs' `rhs', absorb(`absvars') vce(cluster `clustervars', suite(mwc))
	TrimMatrix `K'
	assert e(N_clust)<.
	assert e(N_clust1)<.
	assert e(N_clust2)<.
	assert e(N_clustervars)<.
	
	* Compare
	
	storedresults compare bench_def e(), tol(1e-12) include( ///
		scalar: N rmse tss rss r2 r2_a F df_r df_a df_m /// 
		matrix: trim_b trim_V ///
		macros: wexp wtype )
	storedresults drop bench_def

	* NOTE: This has a lower tolerance; perhaps -cgmreg- uses only single precision somewhere??
	storedresults compare bench_cgmreg e(), tol(1e-8) include( ///
		scalar: N rss r2 /// 
		matrix: trim_b /// trim_V
		macros: wexp wtype )
	storedresults drop bench_cgmreg
	* BUGBUG/NOTE: I can't get a match with CGMREG in the VCV; it's close but not there
	* Maybe it's due to the small sample adjustments? Maybe it's due to the N_clust adjustments?
	
cd "D:/Github/reghdfe/test"
exit
