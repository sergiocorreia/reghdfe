noi cscript "reghdfe: ols with multi-way clustering" adofile reghdfe

* Test multi-way-clustering routine
* Will compare it against ivreg2 and cgmreg.
* Note: We use -fixed_cgmreg- because a bug in -cgmreg-, where they count K=rows(X'X)
* But the rows are also included for omitted variables (e.g. the base variable in a factor expansion i.turn)

* Setup
	qui adopath + "C:/git/reghdfe/test" // allows us to run the fixed_cgmreg.ado

* Dataset
	sysuse auto
	bys turn: gen t = _n
	tsset turn t
	drop if missing(rep)
	
	local included_e1 ///
		scalar: N rss r2 /// rmse tss F df_r ll ll_0 r2_a df_m 
		matrix: trim_b  /// trim_V
		macros: wexp wtype

	local included_e2 ///
		scalar: N rmse rss F df_r ll /// tss r2 r2_a df_m ll_0
		matrix: trim_b trim_V ///
		macros: wexp wtype


* [TEST] Two-way cluster

	local lhs price
	local rhs weight length
	local absvars rep
	local clustervars turn trunk
	fvunab tmp : `rhs'
	local K : list sizeof tmp

	* 1. Run fixed_cgmreg.ado as benchmark
	cap qui tab `absvars', gen(ABS_)
	drop ABS_1
	fixed_cgmreg `lhs' `rhs' ABS_* , cluster(`clustervars') // noeigenfix
	mat li e(V)
	//mat li e(rawcovmat)
	drop ABS_*
	trim_cons `K'
	mat li e(trim_V)
	storedresults save bench_cgmreg e()

	* 2. Run -ivreg2- as benchmark
	cap qui tab `absvars', gen(ABS_)
	ivreg2 `lhs' `rhs' ABS_* , cluster(`clustervars') small nocons partial(ABS_*)
	matrix list e(V)
	drop ABS_*
	trim_cons `K'
	local ivreg2_r2 = e(r2) // Within R2
	storedresults save bench_ivreg2 e()
	
	* 3. Run and test reghdfe-mwc
	reghdfe `lhs' `rhs', absorb(`absvars') vce(cluster `clustervars') keepsing
	matrix list e(V)
	notrim

	* NOTE: This has a lower tolerance; perhaps -cgmreg- uses only single precision somewhere??

	storedresults compare bench_cgmreg e(), include(`included_e1') tol(1e-12)
	storedresults drop bench_cgmreg

	* BUGBUG/NOTE: I can't get a match with CGMREG in the VCV; it's close but not there
	* Maybe it's due to the small sample adjustments? Maybe it's due to the N_clust adjustments?
	* UPDATE: almost surely due to the eigenfix


	storedresults compare bench_ivreg2 e(), include(`included_e2') tol(1e-12)
	storedresults drop bench_ivreg2

	assert "`e(N_clustervars)'"!="" & "`e(N_clustervars)'"!="."
	if abs(e(r2_within)-`ivreg2_r2')>1e-8 {
		di as error "Within R2 doesn't match ivreg2 (`e(r2_within)' vs (`ivreg2_r2')"
		exit 9
	}

	
* [TEST] Three Clustervars
	local clustervars turn trunk mpg
	
	* 1. Run fixed_cgmreg.ado as benchmark
	cap qui tab `absvars', gen(ABS_)
	drop ABS_1
	fixed_cgmreg `lhs' `rhs' ABS_* , cluster(`clustervars')
	drop ABS_*
	trim_cons `K'
	storedresults save bench_cgmreg e()
	
	* 2. Run and test reghdfe
	reghdfe `lhs' `rhs', absorb(`absvars') vce(cluster `clustervars') keepsing
	notrim `K'
	assert e(N_clust)<.
	assert e(N_clust1)<.
	assert e(N_clust2)<.
	assert e(N_clust3)<.
	assert e(N_clustervars)<.
	
	* Compare
	* BUGBUG: SAME PROBLEM AS WITH two-way cluster!
	storedresults compare bench_cgmreg e(), include(`included_e1') tol(1e-12)
	storedresults drop bench_cgmreg

exit
