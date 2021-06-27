* ===========================================================================
* Test a bug related to sorting order
* ===========================================================================

	clear all
	* do create_toy_patent_dataset


// --------------------------------------------------------------------------
// Tiny dataset (just to test it runs)
// --------------------------------------------------------------------------

	clear
	input int(y x) byte(patent inventor id)
	1582 1960 1 1 1
	1582 1960 1 2 1
	1582 1960 1 3 1
	 867  961 2 1 1
	 867  961 2 2 1
	5106 6341 3 1 2
	5106 6341 3 2 2
	5106 6341 3 3 2
	7135 5468 4 2 2
	7135 5468 4 3 2
	3735  568 5 1 2
	3735  568 5 2 2
	  35 1244 6 1 2
	  35 1244 6 3 2
	2421   44 7 3 2
	 317  144 8 1 2
	 917  644 9 2 2
	end

	cls
	list

	sort patent inventor, stable
	reghdfe y x, a(id inventor) group(patent) indiv(inventor) agg(sum) precond(no)
	storedresults save sort1 e()

	set seed 1234
	gen double e = runiform()
	sort e
	li

	reghdfe y x, a(id inventor) group(patent) indiv(inventor) agg(sum) precond(no)
	storedresults compare sort1 e(), tol(1e-10) // include(scalar: r2 rss matrix: trim_b)
	storedresults drop sort1

/*
cls
reghdfe y x, a(id inventor) group(patent) indiv(inventor) agg(sum) keepsing miniter(1) // <whatever>
sort patent inventor
reghdfe y x, a(id inventor) group(patent) indiv(inventor) agg(sum) keepsing miniter(1) // .5008358
egen tag = tag(patent)
tab inventor, gen(DUM)
fcollapse (sum) DUM*, by(patent) fast merge
reg y x sum_DUM* ibn.id if tag // .5008358

*/



// --------------------------------------------------------------------------
// Toy dataset
// --------------------------------------------------------------------------

	use "./toy-patents-long", clear

	sort patent_id inventor_id
	reghdfe citations funding lab_size, a(year inventor_id) group(patent_id) individual(inventor_id) precond(none)
	storedresults save sort1 e()

	sort inventor_id patent_id
	reghdfe citations funding lab_size, a(year inventor_id) group(patent_id) individual(inventor_id) precond(none)
	storedresults save sort2 e()

	set seed 1234
	gen double e = runiform()
	sort e
	reghdfe citations funding lab_size, a(year inventor_id) group(patent_id) individual(inventor_id) precond(none)
	storedresults save unsorted1 e()

	storedresults compare sort1 e(), tol(1e-10) exclude(scalar: ic) // include(scalar: r2 rss matrix: trim_b)
	storedresults compare sort2 e(), tol(1e-10) exclude(scalar: ic) // include(scalar: r2 rss matrix: trim_b)
	storedresults drop sort1 sort2


exit
