cd "D:/Github/reghdfe" // /source
cscript "reghdfe with ivreg2/ivregress and two-step gmm" adofile reghdfe

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
	
* Create fake dataset
	sysuse auto

	local depvar price
	local indepvars weight
	local absvars rep foreign
	local endogvars trunk mpg
	local instruments disp gear head
	local cluster turn

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

	ivreg2 `depvar' `indepvars' ABS_* foreign (`endogvars' = `instruments') , small cue
	TrimMatrix `K'
	storedresults save benchmark e()

	*reghdfe `depvar' `indepvars' (`endogvars'=`instruments'), absorb(`absvars') vce(unadjusted) ivsuite(ivreg2) tol(1e-12) nocon verbose(3) // dropsingletons
	*TrimMatrix `K'
	*storedresults compare benchmark e(), tol(1e-6) include(`include')

	reghdfe `depvar' `indepvars' (`endogvars'=`instruments'), absorb(`absvars') vce(unadjusted) ivsuite(ivreg2) tol(1e-12) nocon cue
	TrimMatrix `K'
	storedresults compare benchmark e(), tol(1e-6) include(`include')

	storedresults drop benchmark

* [TEST] ROBUST

	ivreg2 `depvar' `indepvars' ABS_* foreign (`endogvars' = `instruments') , small cue robust nocons
	TrimMatrix `K'
	storedresults save benchmark e()

	reghdfe `depvar' `indepvars' (`endogvars'=`instruments'), absorb(`absvars') vce(robust) ivsuite(ivreg2) tol(1e-12) cue
	TrimMatrix `K'
	storedresults compare benchmark e(), tol(1e-6) include(`include')

	storedresults drop benchmark

* [TEST] CLUSTER

	ivreg2 `depvar' `indepvars' ABS_* foreign (`endogvars' = `instruments') , small cue cluster(`cluster') nocons
	TrimMatrix `K'
	storedresults save benchmark e()

	reghdfe `depvar' `indepvars' (`endogvars'=`instruments'), absorb(`absvars') vce(cluster `cluster') ivsuite(ivreg2) tol(1e-12) cue
	TrimMatrix `K'
	storedresults compare benchmark e(), tol(1e-6) include(`include')

	storedresults drop benchmark

* [TEST] MWC

	local cluster turn cluster2

	ivreg2 `depvar' `indepvars' ABS_* foreign (`endogvars' = `instruments') , small cue cluster(`cluster') nocons partial(ABS_* foreign)
	TrimMatrix `K'
	storedresults save benchmark e()

	reghdfe `depvar' `indepvars' (`endogvars'=`instruments'), absorb(`absvars') vce(cluster `cluster') ivsuite(ivreg2) tol(1e-12) cue
	TrimMatrix `K'
	storedresults compare benchmark e(), tol(1e-6) include(`include')

	storedresults drop benchmark

cd "D:/Github/reghdfe/test"
exit
