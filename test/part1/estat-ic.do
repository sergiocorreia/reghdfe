noi cscript "reghdfe postestimation: estat ic" adofile reghdfe

	sysuse auto
	areg price weight length foreign, a(turn)
	estat ic
	return list
	storedresults save benchmark r()

	reghdfe price weight length foreign, a(turn) keepsing
	estat ic
	return list
	storedresults compare benchmark r() // , tol(1e-12) include(`included_e')

storedresults drop benchmark
exit
