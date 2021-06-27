* ===========================================================================
* Test reghdfe preconditioner for individual FEs
* ===========================================================================

	clear all

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

	reghdfe y x, a(inventor) group(patent) indiv(inventor) precond(no) keepsing
	reghdfe y x, a(inventor) group(patent) indiv(inventor) precond(diag) keepsing


// --------------------------------------------------------------------------
// Test group() option
// --------------------------------------------------------------------------

	use "./toy-patents-long", clear
	reghdfe citations funding lab_size, a(year) group(patent_id) precond(no)
	storedresults save pre_no e()

	reghdfe citations funding lab_size, a(year) group(patent_id) precond(diag)
	storedresults compare pre_no e(), tol(1e-10) include(scalar: r2 rss)
	storedresults drop pre_no

// --------------------------------------------------------------------------
// Test individual() option
// --------------------------------------------------------------------------

	use "./toy-patents-long", clear
	reghdfe citations funding lab_size, a(c inventor_id) group(patent_id) individual(inventor_id) precond(no) dof(none)
	trim_coefs 2
	storedresults save pre_no e()
	
	reghdfe citations funding lab_size, a(c inventor_id) group(patent_id) individual(inventor_id) precond(diag) dof(none)
	trim_coefs 2
	storedresults compare pre_no e(), tol(1e-10) include(scalar: r2 rss matrix: trim_b)
	storedresults drop pre_no


exit
