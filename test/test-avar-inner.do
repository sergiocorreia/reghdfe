* CALL THIS FROM test-avar
args i
di as text "{title:i=`i'}"
assert inlist(`i',1,2)

* [TEST] Will compare reghdfe against ivreg2
	
	local lhs price
	local rhs weight length gear
	local absvars trunk
	qui tab `absvars', gen(ABS_)
	fvunab tmp : `rhs'
	local K : list sizeof tmp

	local included ///
		scalar: N rmse rss F df_r /// r2 r2_a mss df_a df_m 
		matrix: trim_b trim_V ///
		macros: wexp wtype clustvar

	local included_subset /// Used when not full rank due to clusters
		scalar: N rmse rss df_r /// F
		matrix: trim_b /// trim_V
		macros: wexp wtype clustvar

	* Unadjusted
	ivreg2 `lhs' `rhs' ABS_*, small partial(ABS_*)
	TrimMatrix `K'
	storedresults save benchmark e()

	reghdfe `lhs' `rhs', absorb(`absvars') nocons vce(unadjusted)
	TrimMatrix `K'
	storedresults compare benchmark e(), tol(1e-6) include(`included')

	reghdfe `lhs' `rhs', absorb(`absvars') nocons vce(unadjusted, suite(avar))
	TrimMatrix `K'
	storedresults compare benchmark e(), tol(1e-6) include(`included')

	storedresults drop benchmark

	* Robust
	ivreg2 `lhs' `rhs' ABS_*, small partial(ABS_*) robust
	TrimMatrix `K'
	storedresults save benchmark e()

	reghdfe `lhs' `rhs', absorb(`absvars') nocons vce(robust)
	TrimMatrix `K'
	storedresults compare benchmark e(), tol(1e-6) include(`included')

	reghdfe `lhs' `rhs', absorb(`absvars') nocons vce(robust, suite(avar))
	TrimMatrix `K'
	storedresults compare benchmark e(), tol(1e-6) include(`included')
	
	storedresults drop benchmark

	* Autoco
	ivreg2 `lhs' `rhs' ABS_*, small partial(ABS_*) bw(2)
	TrimMatrix `K'
	storedresults save benchmark e()

	reghdfe `lhs' `rhs', absorb(`absvars') nocons vce(, bw(2))
	TrimMatrix `K'
	storedresults compare benchmark e(), tol(1e-6) include(`included')

	reghdfe `lhs' `rhs', absorb(`absvars') nocons vce(, bw(2) suite(avar))
	TrimMatrix `K'
	storedresults compare benchmark e(), tol(1e-6) include(`included')

	storedresults drop benchmark

	* Autoco + Kernel
	ivreg2 `lhs' `rhs' ABS_*, small partial(ABS_*) bw(2) kernel(par)
	TrimMatrix `K'
	storedresults save benchmark e()
	reghdfe `lhs' `rhs', absorb(`absvars') nocons vce(, bw(2) kernel(par))
	TrimMatrix `K'
	storedresults compare benchmark e(), tol(1e-6) include(`included')
	storedresults drop benchmark

	* dkraay
	ivreg2 `lhs' `rhs' ABS_*, small partial(ABS_*) dkraay(3)
	TrimMatrix `K'
	storedresults save benchmark e()
	reghdfe `lhs' `rhs', absorb(`absvars') nocons vce(, dkraay(3))
	TrimMatrix `K'
	storedresults compare benchmark e(), tol(1e-6) include(`included')
	storedresults drop benchmark

	* kiefer
	ivreg2 `lhs' `rhs' ABS_*, small partial(ABS_*) kiefer
	TrimMatrix `K'
	storedresults save benchmark e()
	reghdfe `lhs' `rhs', absorb(`absvars') nocons vce(, kiefer)
	TrimMatrix `K'
	storedresults compare benchmark e(), tol(1e-6) include(`included')
	storedresults drop benchmark

	* Cluster
	ivreg2 `lhs' `rhs' ABS_*, small partial(ABS_*) cluster(turn t`i')
	TrimMatrix `K'
	local fullrank = (e(rankS)==e(rankzz))
	storedresults save benchmark e()

	reghdfe `lhs' `rhs', absorb(`absvars') nocons vce(cluster turn t`i', suite(default)) dof(none) tol(1e-12)
	TrimMatrix `K'
	if (`fullrank') storedresults compare benchmark e(), tol(1e-6) include(`included')
	if (!`fullrank') storedresults compare benchmark e(), tol(1e-6) include(`included_subset')

	reghdfe `lhs' `rhs', absorb(`absvars') nocons vce(cluster turn t`i', suite(mwc)) dof(none) tol(1e-12)
	TrimMatrix `K'
	if (`fullrank') storedresults compare benchmark e(), tol(1e-6) include(`included')
	if (!`fullrank') storedresults compare benchmark e(), tol(1e-6) include(`included_subset')

	reghdfe `lhs' `rhs', absorb(`absvars') nocons vce(cluster turn t`i', suite(avar)) dof(none) tol(1e-12)
	TrimMatrix `K'
	if (`fullrank') storedresults compare benchmark e(), tol(1e-6) include(`included')
	if (!`fullrank') storedresults compare benchmark e(), tol(1e-6) include(`included_subset')

	storedresults drop benchmark

	* Cluster + Autoco
	ivreg2 `lhs' `rhs' ABS_*, small partial(ABS_*) cluster(turn t`i') bw(5)
	TrimMatrix `K'
	local fullrank = (e(rankS)==e(rankzz))
	storedresults save benchmark e()

	reghdfe `lhs' `rhs', absorb(`absvars') nocons vce(cluster turn t`i', bw(5)) dof(none)
	TrimMatrix `K'
	if (`fullrank') storedresults compare benchmark e(), tol(1e-6) include(`included')
	if (!`fullrank') storedresults compare benchmark e(), tol(1e-6) include(`included_subset')
	storedresults drop benchmark

	* Cluster + Autoco + kernel
	ivreg2 `lhs' `rhs' ABS_*, small partial(ABS_*) cluster(turn t`i') bw(5) kernel(thann)
	TrimMatrix `K'
	storedresults save benchmark e()
	reghdfe `lhs' `rhs', absorb(`absvars') nocons vce(cluster turn t`i', bw(5) kernel(thann)) dof(none)
	TrimMatrix `K'
	if (`fullrank') storedresults compare benchmark e(), tol(1e-6) include(`included')
	if (!`fullrank') storedresults compare benchmark e(), tol(1e-6) include(`included_subset')
	storedresults drop benchmark
