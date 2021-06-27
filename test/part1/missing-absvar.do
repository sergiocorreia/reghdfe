noi cscript "reghdfe: ols with robust VCE" adofile reghdfe

* Dataset
	sysuse auto
	bys turn: gen t = _n
	tsset turn t
	drop if missing(rep)
	replace foreign = . if foreign
	
* [TEST] Ensure we drop obs. where absvar is missing
	
	cou if !missing(foreign)
	loc N = r(N)
	reghdfe price weight gear, absorb(foreign) keepsing
	assert e(N)==`N'

exit
