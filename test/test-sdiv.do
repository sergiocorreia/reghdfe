*cscript "reghdfe with (experimental) static first stage, dynamic second stage IV" adofile reghdfe
cls

* Setup
	discard
	clear all
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

* Variables
	local depvar y
	local endogvars z1 z2
	local exogvars x1 x2
	local instruments w1 w2 w3

* Fake dataset
	clear
	local T = 10 // 100
	local N = 100 // 1e5
	local NT = `N' * `T'
	set obs `NT'
	gen long id = 1 + int((_n-1)/`N')
	su id
	bys id: gen int t = _n
	compress
	tsset id t

	* Fixed Effects
	bys id: gen fe1 = rnormal()
	bys id: replace fe1 = fe1[1]

	bys t: gen fe2 = rnormal()
	bys t: replace fe2 = fe2[1]

	* True model:
	* 1) One first stage
	* 2) Dynamic second stage

	gen double u = rnormal()
	gen double v = rnormal()
	gen double latent = rnormal()

	gen double iv = rnormal() + latent
	gen double z  = rnormal() + 0.4 * F.iv + 1.2 * iv + 0.8 * L.iv ///
		+ u + 0.2*L.u + 0.1*L2.u + 0.05*L3.u + 0.05*L4.u

	gen double x = rnormal() + 0.5 * latent
	gen double y = 10 + x ///
		+ z + 0.9*L.z + 0.8*L2.z + 0.7*L3.z + 0.6*L4.z + 0.5*L5.z + 0.4*L6.z + 0.3*L7.z + 0.2*L8.z + 0.1*L9.z ///
		+ u + 0.3*L.u + 0.2*L2.u + 0.1*L3.u




	foreach var in `exogvars' `instruments' error1 error2 {
		gen double `var' = rnormal()
	}
	gen `endogvar' = 1*gear + 1*disp + 1*foreign + 1*rep + 1*length + 15*error1 + 10*error2 + 20 * fe1
	gen `depvar' = 100*`endogvar' + 100*weight + 100*length + 1000*error2 + 20 * fe2

	compress

asdasd
	fvunab tmp : `indepvars' `endogvars'
	local K : list sizeof tmp

	local include  ///
		scalar: N df_r /// rmse rss r2
		matrix: trim_b trim_V ///
		macros: wexp wtype

	qui tab rep, gen(ABS_)

	gen byte cluster2 =  int((1+_n)/10)
	*qui tab cluster2, gen(ABS2_)

* [TEST] UNADJUSTED

	ivreg2 `depvar' `indepvars' ABS_* foreign (`endogvars' = `instruments') , small gmm2s
	TrimMatrix `K'
	storedresults save benchmark e()

	reghdfe `depvar' `indepvars' (`endogvars'=`instruments'), absorb(`absvars') vce(unadjusted) ivsuite(ivreg2) tol(1e-12) est(gmm2s)
	TrimMatrix `K'
	storedresults compare benchmark e(), tol(1e-10) include(`include')

	* Just for reference
	ivregress gmm `depvar' `indepvars' ABS_* foreign (`endogvars'=`instruments'), vce(unadjusted) wmatrix(unadjusted) small

	reghdfe `depvar' `indepvars' (`endogvars'=`instruments'), absorb(`absvars') vce(unadjusted) ivsuite(ivregress) tol(1e-12) est(gmm2s) v(3)
	TrimMatrix `K'
	storedresults compare benchmark e(), tol(1e-8) include(`include')

	storedresults drop benchmark

* [TEST] ROBUST

	ivreg2 `depvar' `indepvars' ABS_* foreign (`endogvars' = `instruments') , small gmm2s robust nocons
	TrimMatrix `K'
	storedresults save benchmark e()

	reghdfe `depvar' `indepvars' (`endogvars'=`instruments'), absorb(`absvars') vce(robust) ivsuite(ivreg2) tol(1e-12) estimator(gmm2s)
	TrimMatrix `K'
	storedresults compare benchmark e(), tol(1e-10) include(`include')

	ivregress gmm `depvar' `indepvars' ABS_* foreign (`endogvars' = `instruments') , small wmatrix(robust) vce(unadjusted)
	TrimMatrix `K'
	storedresults compare benchmark e(), tol(1e-9) include(`include')

	reghdfe `depvar' `indepvars' (`endogvars'=`instruments'), absorb(`absvars') vce(robust) ivsuite(ivregress) tol(1e-12) est(gmm2s)  vceunadjusted
	TrimMatrix `K'
	storedresults compare benchmark e(), tol(1e-9) include(`include')

	storedresults drop benchmark

* [TEST] TWICE ROBUST - DOES NOT GIVE EXACT RESULTS FOR VCE!!!!

	ivregress gmm `depvar' `indepvars' ABS_* foreign (`endogvars' = `instruments') , small wmatrix(robust) vce(robust)
	TrimMatrix `K'
	storedresults save benchmark e()

	ivregress gmm `depvar' `indepvars' ABS_* foreign (`endogvars' = `instruments') , small // wmatrix(robust)
	TrimMatrix `K'
	storedresults compare benchmark e(), tol(1e-10) include(`include')

	* How it looks without it
	reghdfe `depvar' `indepvars' (`endogvars'=`instruments'), absorb(`absvars') vce(robust) ivsuite(ivregress) estimator(gmm2s)

	reghdfe `depvar' `indepvars' (`endogvars'=`instruments'), absorb(`absvars') vce(robust, twicerobust) ivsuite(ivregress) estimator(gmm2s)
	TrimMatrix `K'
	storedresults compare benchmark e(), tol(1e-8) include(scalar: N df_r /// rmse rss r2
		matrix: trim_b ///
		macros: wexp wtype)
	*storedresults compare benchmark e(), tol(1e-2) include(matrix: trim_V) // BUGBUG!!!!!!!!!!!!!!!!!!!!!!!!!

	storedresults drop benchmark
	
	* BUGBUG: I expected this to work but fails; so I'm setting vceunadjusted to TRUE always
	**ivregress gmm `depvar' `indepvars' ABS_* foreign (`endogvars' = `instruments') , small wmatrix(robust) vce(robust) hascons // This one is different b/c it uses vce same as wmatrix!
	**TrimMatrix `K'
	**storedresults save benchmark e()
	**
	**reghdfe `depvar' `indepvars' (`endogvars'=`instruments'), absorb(`absvars') vce(robust) ivsuite(ivregress) tol(1e-12) gmm2s
	**TrimMatrix `K'
	**storedresults compare benchmark e(), tol(1e-6) include(`include')
	**
	**storedresults drop benchmark

* [TEST] CLUSTER

* BUGBUG/WARNING: results do not have good enough precision
* Is that due to the ivreg2.partial bug???

	ivreg2 `depvar' `indepvars' ABS_* foreign (`endogvars' = `instruments') , small gmm2s cluster(`cluster') nocons
	TrimMatrix `K'
	storedresults save benchmark e()

	reghdfe `depvar' `indepvars' (`endogvars'=`instruments'), absorb(`absvars') vce(cluster `cluster') ivsuite(ivreg2) tol(1e-12) est(gmm2s)
	TrimMatrix `K'
	storedresults compare benchmark e(), tol(1e-6) include(`include')

	ivregress gmm `depvar' `indepvars' ABS_* foreign (`endogvars' = `instruments') , small wmatrix(cluster `cluster') // This one is different b/c it uses vce same as wmatrix!
	ivregress gmm `depvar' `indepvars' ABS_* foreign (`endogvars' = `instruments') , small wmatrix(cluster `cluster') vce(unadjusted)
	*TrimMatrix `K'
	*storedresults save benchmark e()

	reghdfe `depvar' `indepvars' (`endogvars'=`instruments'), absorb(`absvars') vce(cluster `cluster') ivsuite(ivregress) tol(1e-12) est(gmm2s) vceunadjusted
	TrimMatrix `K'
	storedresults compare benchmark e(), tol(1e-6) include(`include')

	storedresults drop benchmark

* [TEST] MWC

	local cluster turn cluster2

	ivreg2 `depvar' `indepvars' ABS_* foreign (`endogvars' = `instruments') , small gmm2s cluster(`cluster') nocons partial(ABS_* foreign)
	TrimMatrix `K'
	storedresults save benchmark e()

	reghdfe `depvar' `indepvars' (`endogvars'=`instruments'), absorb(`absvars') vce(cluster `cluster') ivsuite(ivreg2) tol(1e-12) est(gmm2s)
	TrimMatrix `K'
	storedresults compare benchmark e(), tol(1e-6) include(`include')

	storedresults drop benchmark

cd "C:/Git/reghdfe/test"
exit
