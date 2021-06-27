noi cscript "reghdfe: cloning the data shouldn't help" adofile reghdfe
* For consistency, running vce(cluster _n) should be the same as running vce(robust)

* Dataset
	sysuse auto
	bys turn: gen t = _n
	tsset turn t
	drop if missing(rep)
	gen i = _n
	
	local included_e ///
		macros: wexp wtype ///
		matrix: trim_b trim_V ///
		scalar: N rmse tss rss r2 r2_a F df_m ll ll_0 // mss is not reported by areg if clustervar!=absvar

* [TEST] Cluster is absvar (reghdfe robust vs reghdfe cluster)
	local lhs price
	local rhs weight length
	local absvars turn
	local clustervar i
	fvunab tmp : `rhs'
	local K : list sizeof tmp
	
	* 1. Run benchmark
	reghdfe `lhs' `rhs', absorb(`absvars') keepsingletons verbose(-1) vce(robust)
	storedresults save benchmark e()

	* Expand dataset without gaining info
	expand 50
	
	* 2. Run reghdfe
	reghdfe `lhs' `rhs', absorb(`absvars') keepsingletons verbose(-1) vce(cluster `clustervar')
	loc exclude macros: vce cmdline
	storedresults compare benchmark e(), tol(1e-12) exclude(`exclude')
	loc exclude `exclude' clustvar1 clustvar title3 scalar: N_clustervars N_clust1
	storedresults compare benchmark e(), tol(1e-12) reverse exclude(`exclude')

	storedresults drop benchmark

exit
