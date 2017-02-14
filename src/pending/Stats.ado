// -----------------------------------------------------------------------------
// Matrix of summary statistics
// -----------------------------------------------------------------------------
cap pr drop Stats
pr Stats
	syntax varlist(numeric), [weightexp(string)] stats(string) statsmatrix(string)

	* Optional weights
	mata: st_local("weight", sprintf("[%s=%s]", HDFE.weighttype, HDFE.weightvar))
	assert "`weight'" != ""
	if ("`weight'" == "[=]") loc weight
	loc weight : subinstr local weight "[pweight" "[aweight"

	mata: st_local("summarize_stats", HDFE.options.summarize_stats)
	mata: st_local("varlist", HDFE.options.varlist)
	mata: st_local("cvars", HDFE.cvars)
	loc full_varlist `varlist' `cvars'

	qui tabstat `full_varlist' `weight' , stat(`stats') col(stat) save
	matrix reghdfe_statsmatrix = r(StatTotal)


matrix list reghdfe_statsmatrix

asdasdsad
g
	* Fix names (__L__.price -> L.price)
	local colnames : colnames `statsmatrix'
	FixVarnames `colnames'
	local colnames "`r(newnames)'"
	matrix colnames `statsmatrix' = `colnames'
end
