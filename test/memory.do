noi cscript "reghdfe: memory-related options" adofile reghdfe

* TEST 1: One set of FEs

	* Dataset
	sysuse auto, clear
	bys turn: gen int t = _n
	xtset turn t
	gen byte c = 1
	gen int i = _n

	* Options	
	loc cmd 	reghdfe price length i.rep##c.weight gear [fw=c], a(FE=trunk) vce(cluster trunk#turn i) tol(1e-10) resid summ
	loc predict predict resid, resid
	loc exclude macro: cmdline // scalar: ic

	* Benchmark
	`cmd' nocompact pool(0)
	storedresults save benchmark e()
	`predict'
	rename (FE resid _reghdfe_resid) B_=

	* Low-memory alternative
	`cmd' compact pool(1)
	storedresults compare benchmark e(), tol(1e-12) exclude(`exclude')
	`predict'

	* Assertions
	_vassert FE B_FE
	_vassert resid B_resid
	_vassert _reghdfe_resid B__reghdfe_resid
	storedresults drop benchmark


* TEST 2: Two sets of FEs
* With 2+ FEs, each col might iterate a different number of times, so results might vary slightly with pool(#)
* In particular, e(V) and e(F) might vary slightly

	* Dataset
	sysuse auto, clear
	bys turn: gen int t = _n
	xtset turn t
	gen byte c = 1
	gen int i = _n

	* Options	
	loc cmd 	reghdfe price length i.rep##c.weight gear [fw=c], a(trunk turn) vce(cluster trunk#turn i) tol(1e-12) // v(1)
	loc exclude macro: cmdline // scalar: ic

	* Benchmark
	`cmd' nocompact pool(0)
	storedresults save benchmark e()

	* Low-memory alternative
	`cmd' compact pool(1)
	storedresults compare benchmark e(), tol(1e-12) exclude(`exclude')
	storedresults drop benchmark


* TEST 2: Two sets of FEs - LSMR method

	* Dataset
	sysuse auto, clear
	bys turn: gen int t = _n
	xtset turn t
	gen byte c = 1
	gen int i = _n

	* Options	
	loc cmd 	reghdfe price length i.rep##c.weight gear [fw=c], a(trunk turn) vce(cluster trunk#turn i) tol(1e-12) accel(lsmr) // v(1)
	loc exclude macro: cmdline // scalar: ic

	* Benchmark
	`cmd' nocompact pool(0)
	storedresults save benchmark e()

	* Low-memory alternative
	`cmd' compact pool(1)
	storedresults compare benchmark e(), tol(1e-12) exclude(`exclude')
	storedresults drop benchmark


exit
