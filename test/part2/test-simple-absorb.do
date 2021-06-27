* ===========================================================================
* Simple tests of a() to compare against reghdfe
* ===========================================================================

	clear all
	sysuse auto
	set trace off
	set tracedepth 2

	loc iter 0 // just to keep track of how many specs we tested in total

	loc absorbs `""noa" "a(turn)" `"a(turn trunk)"' `a(turn#foreign trunk)'"'
	loc techniques map lsmr lsqr
	loc transforms cimmino symmetric_kaczmarz

	* LOOP
	foreach absorb of local absorbs {

	* Benchmark
	reghdfe price weight length, `absorb' tol(1e-12) accel(cg) transf(sym)
	trim_coefs 2
	storedresults save benchmark e()

	* LOOP
	foreach technique of local techniques {

	* LOOP
	foreach transform of local transforms {

	if ("`technique'" != "map" & "`transform'" != "cimmino") continue

		loc ++iter
		di as text _n "{bf}{hline 64}"
		di as text `"{bf}(`iter') absorb="`absorb'" technique=`technique' transform=`transform'"'
		di as text "{bf}{hline 64}" _n

		loc cmd `"reghdfe price weight length, `absorb' tol(1e-12) tech(`technique') transform(`transform') accel(cg)"'
		di as text `"CMD: {inp}`cmd'"'
		`cmd'
		trim_coefs 2
		storedresults compare benchmark e(), tol(1e-10) include(scalar: r2 rss matrix: trim_b trim_V)

	}
	}
		storedresults drop benchmark
	}

	di as text "DONE! Tested `iter' specifications"
exit
