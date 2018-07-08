*do setup
noi cscript "reghdfe: test -prune- option for correctness" adofile reghdfe

loc cmd reghdfe y x, a(id1 id2) tol(1e-10)
loc exclude macro: cmdline scalar: ic

* [TESTS]
forv i=1/5 {
	use prune`i', clear
	di as input "[CMD] `cmd'"
	
	`cmd' noprune
	loc ic`i'a = e(ic)
	storedresults save benchmark e()
	
	`cmd' prune
	loc ic`i'b = e(ic)
	storedresults compare benchmark e(), tol(1e-12) exclude(`exclude')
	storedresults drop benchmark
	
	assert `ic`i'a' >= `ic`i'b' // pruning reduces the number of iterations
}

di as text "IC:"
forv i=1/5 {
	di as text "`ic`i'a' | `ic`i'b'"
}

exit
