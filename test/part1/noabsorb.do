noi cscript "reghdfe: noabsorb option" adofile reghdfe

* Dataset
	sysuse auto
	bys turn: gen t = _n
	tsset turn t
	drop if missing(rep)
	gen byte c = 1
	
	local excluded ///
		macros: cmdline extended_absvars absvars

* [TEST]

	local lhs price
	local rhs weight length
	local absvars c
	fvunab tmp : `rhs'
	local K : list sizeof tmp
	
	* 1. Run benchmark
	reghdfe `lhs' `rhs', absorb(`absvars') keepsingletons verbose(-1)
	storedresults save benchmark e()
	
	* 2. Run reghdfe
	reghdfe `lhs' `rhs', noabsorb keepsingletons verbose(-1)
	storedresults compare benchmark e(), tol(1e-12) exclude(`excluded')


storedresults drop benchmark
exit
