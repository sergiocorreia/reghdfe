* ===========================================================================
* Test how individual fixed effects get aggregated: agg(sum) & agg(mean)
* ===========================================================================

	clear all

	* do create_toy_patent_dataset


// --------------------------------------------------------------------------
// Tiny dataset (just to test it runs)
// --------------------------------------------------------------------------

	clear
	input int(y x) byte(patent inventor)
	1582 1960 1 1
	1582 1960 1 2
	1582 1960 1 3
	 867  961 2 4
	5106 6341 3 1
	5106 6341 3 2
	5106 6341 3 3
	7135 5468 4 5
	7135 5468 4 4
	end
	list
	sort patent inventor, stable

	reghdfe y x, a(inventor) group(patent) indiv(inventor) agg(sum)
	reghdfe y x, a(inventor) group(patent) indiv(inventor) agg(mean)


// --------------------------------------------------------------------------
// Toy dataset
// --------------------------------------------------------------------------
// Test for:
// a) With absorb(constant IndivFE) vs absorb(IndivFE year)
// b) With function(sum) and function(mean)
// c) With preconditioner(none) and preconditioner(diagonal)
// d) With keepsingletons on and off 

	loc iter 0 // just to keep track of how many specs we tested in total

	* LOOP
	forval i = 1/2 {
	if (`i'==1) loc absorb_long "c inventor_id"
	if (`i'==2) loc absorb_long "inventor_id year"
	if (`i'==1) loc absorb_wide ""
	if (`i'==2) loc absorb_wide "a(year)"

	* LOOP
	foreach preconditioner in none diagonal {

	* LOOP
	foreach ks in "" "keepsingletons" {

	* LOOP
	foreach fun in sum average {
	if ("`fun'" == "sum") loc wide_fn "toy-patents-wide"
	if ("`fun'" == "average") loc wide_fn "toy-patents-wide-average"
		
		loc ++iter
		di as text _n "{bf}{hline 64}"
		di as text "{bf}(`iter') absorb=`i' fun=`fun' precond=`preconditioner' keepsing=<`keepsingletons'>"
		di as text "{bf}{hline 64}" _n

		use "./toy-patents-long", clear
		reghdfe citations funding lab_size, a(`absorb_long') group(patent_id) individual(inventor_id) aggregation(`fun') precond(`preconditioner') `ks'
		trim_coefs 2
		storedresults save long1 e()

		use "./`wide_fn'", clear
		reghdfe citations funding lab_size dummy*, `absorb_wide' `ks'
		trim_coefs 2
		storedresults compare long1 e(), tol(1e-10) include(scalar: r2 rss matrix: trim_b)
		storedresults drop long1

	}
	}
	}
	}

	di as text "DONE! Tested `iter' specifications"
exit
