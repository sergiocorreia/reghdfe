noi cscript "reghdfe: bug in st_data" adofile reghdfe

* Dataset
	sysuse auto
	recast double length weight
	bys turn: gen int t = _n
	tsset turn t
	drop if missing(rep)

	local excluded ///
		macros: cmdline title marginsok marginsprop cmd predict model estat_cmd ///
		scalar: rank


* [TEST] Prevent regression of Stata 11-14 issue
* Here we stay below the limit of four parenthesis

	local lhs price
	local rhs 1.rep78 2.rep78 3.rep78#1.foreign

	* 1. Run benchmark
	reg `lhs' `rhs'
	storedresults save benchmark e()

	* 2. Run reghdfe
	reghdfe `lhs' `rhs', noa v(1)
	storedresults compare benchmark e(), tol(1e-9) exclude(`excluded')

	* 2. Run reghdfe with pool(1)
	reghdfe `lhs' `rhs', noa pool(1) v(1)
	storedresults compare benchmark e(), tol(1e-9) exclude(`excluded')

	* 3. Run reghdfe with pool(5)
	reghdfe `lhs' `rhs', noa pool(5) v(1)
	storedresults compare benchmark e(), tol(1e-9) exclude(`excluded')

	* Done!
	storedresults drop benchmark


* [TEST] Prevent regression of Stata 11-14 issue
	* Here we go above the limit of four parenthesis

	local lhs price
	local rhs 31.turn#c.weight 43.turn#c.weight 4.rep#c.weight 1.rep78 2.rep78 3.rep78#1.foreign 

	* 1. Run benchmark
	reg `lhs' `rhs'
	storedresults save benchmark e()

	* 2. Run reghdfe
	reghdfe `lhs' `rhs', noa v(1)
	storedresults compare benchmark e(), tol(1e-8) exclude(`excluded')

	* 2. Run reghdfe with pool(1)
	reghdfe `lhs' `rhs', noa pool(1) v(1)
	storedresults compare benchmark e(), tol(1e-8) exclude(`excluded')

	* 3. Run reghdfe with pool(5)
	reghdfe `lhs' `rhs', noa pool(5) v(1)
	storedresults compare benchmark e(), tol(1e-8) exclude(`excluded')

	* Done!
	storedresults drop benchmark


exit
