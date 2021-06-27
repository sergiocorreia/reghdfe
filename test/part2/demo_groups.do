* ===========================================================================
* Test options specific to reghdfe with individual FEs
* ===========================================================================

	clear all

	* do create_toy_patent_dataset


// --------------------------------------------------------------------------
// Test singletons dropping
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

	reghdfe y x, a(inventor) group(patent) indiv(inventor) precond(none) keepsing


// --------------------------------------------------------------------------
// Test group() option
// --------------------------------------------------------------------------

	clear

	use "./toy-patents-long"
	reghdfe citations funding lab_size, a(year) group(patent_id) aggreg(sum)
	use "./toy-patents-wide"
	reghdfe citations funding lab_size, a(year) aggreg(sum)

	use "./toy-patents-long"
	reghdfe citations funding lab_size [fw=c], a(year) group(patent_id) aggreg(sum)
	use "./toy-patents-wide"
	reghdfe citations funding lab_size [fw=c], a(year) aggreg(sum)
	
	use "./toy-patents-long"
	reghdfe citations funding lab_size [aw=c], a(year) group(patent_id) aggreg(sum)
	use "./toy-patents-wide"
	reghdfe citations funding lab_size [aw=c], a(year) aggreg(sum)
	
	use "./toy-patents-long"
	reghdfe citations funding lab_size [pw=w1], a(year) group(patent_id) aggreg(sum)
	use "./toy-patents-wide"
	reghdfe citations funding lab_size [pw=w1], a(year) aggreg(sum)
	
	use "./toy-patents-long"
	reghdfe citations funding lab_size [fw=w2], a(year) group(patent_id) aggreg(sum)
	use "./toy-patents-wide"
	reghdfe citations funding lab_size [fw=w2], a(year) aggreg(sum)


// --------------------------------------------------------------------------
// Test individual() option
// --------------------------------------------------------------------------

	use "./toy-patents-long", clear

	reghdfe citations, a(c inventor_id) group(patent_id) individual(inventor_id) dof(none) keepsing precond(none) aggreg(sum)
	storedresults save long1 e()

	reghdfe citations c, a(inventor_id) group(patent_id) individual(inventor_id) dof(none) keepsing precond(none) aggreg(sum)
	storedresults save long2 e()

	use "./toy-patents-wide", clear
	reghdfe citations dummy*, noa keepsing precond(no)
	di e(rss) // 385.577
	reghdfe citations dummy*, a(c) keepsing
	di e(rss) // 385.577

	storedresults compare long1 e(), tol(1e-10) include(scalar: r2 rss)
	storedresults compare long2 e(), tol(1e-10) include(scalar: r2 rss)
	storedresults drop long1 long2

// --------------------------------------------------------------------------

	use "./toy-patents-long", clear
	sort patent inventor
	reghdfe citations funding lab_size, a(c inventor_id) group(patent_id) individual(inventor_id) dof(none) keepsing precond(none) aggreg(sum)
	trim_coefs 2
	storedresults save long1 e()

	use "./toy-patents-wide", clear
	reghdfe citations funding lab_size dummy*, noa keepsing
	trim_coefs 2
	storedresults save wide1 e()

	reghdfe citations funding lab_size dummy*, noa keepsing
	trim_coefs 2
	storedresults compare long1 e(), tol(1e-10) include(scalar: r2 rss N matrix: trim_b)
	storedresults compare wide1 e(), tol(1e-10) include(scalar: r2 rss N matrix: trim_b)
	storedresults drop long1 wide1


// --------------------------------------------------------------------------

	use "./toy-patents-long", clear
	sort patent inventor
	reghdfe citations funding lab_size, a(inventor_id year) group(patent_id) individual(inventor_id) dof(none) keepsing precond(none) aggreg(sum)
	trim_coefs 2
	storedresults save long1 e()

	use "./toy-patents-wide", clear
	reghdfe citations funding lab_size dummy*, a(year) keepsing
	trim_coefs 2
	storedresults save wide1 e()

	reghdfe citations funding lab_size dummy*, a(year) keepsing
	trim_coefs 2
	storedresults compare long1 e(), tol(1e-10) include(scalar: r2 rss N matrix: trim_b)
	storedresults compare wide1 e(), tol(1e-10) include(scalar: r2 rss N matrix: trim_b)
	storedresults drop long1 wide1


// --------------------------------------------------------------------------

	use "./toy-patents-long", clear
	sort patent inventor
	// IT SEEMS I NEED A MUCH MORE STRICTER TOL FOR LSQR!??!?
	reghdfe citations funding lab_size, a(inventor_id year) tech(lsqr) tol(1e-14) group(patent_id) individual(inventor_id) dof(none) keepsing precond(none) aggreg(sum)
	trim_coefs 2
	storedresults save long1 e()

	use "./toy-patents-wide", clear
	reghdfe citations funding lab_size dummy*, a(year) keepsing
	trim_coefs 2
	storedresults save wide1 e()

	reghdfe citations funding lab_size dummy*, a(year) keepsing
	trim_coefs 2
	storedresults compare long1 e(), tol(1e-10) include(scalar: r2 rss N matrix: trim_b)
	storedresults compare wide1 e(), tol(1e-10) include(scalar: r2 rss N matrix: trim_b)
	storedresults drop long1 wide1

// --------------------------------------------------------------------------

	* use "./toy-patents-long", clear
	* sort patent inventor
	* reghdfe citations funding lab_size, a(inventor_id##c.age year) group(patent_id) individual(inventor_id) dof(none) keepsing precond(none)
	* * unsure how to run in wide

exit
