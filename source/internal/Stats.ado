// -----------------------------------------------------------------------------
// Matrix of summary statistics
// -----------------------------------------------------------------------------
capture program drop Stats
program define Stats
	syntax varlist(numeric), [weightexp(string)] stats(string) statsmatrix(string) [USEcache]

	if ("`usecache'"=="") {
		local tabstat_weight : subinstr local weightexp "[pweight" "[aweight"
		qui tabstat `varlist' `tabstat_weight' , stat(`stats') col(stat) save
		matrix `statsmatrix' = r(StatTotal)

		* Fix names (__L__.price -> L.price)
		local colnames : colnames `statsmatrix'
		FixVarnames `colnames'
		local colnames "`r(newnames)'"
		matrix colnames `statsmatrix' = `colnames'
	}
	else {
		cap conf matrix reghdfe_statsmatrix
		
		* Fix names
		FixVarnames `varlist'
		local sample_names "`r(newnames)'"

		* Trim matrix
		local all_names : colnames reghdfe_statsmatrix
		local first 1 // 1 if `statsmatrix' is still empty
		foreach name of local all_names {
			local is_match : list name in sample_names
			if (`is_match' & `first') {
				local first 0
				matrix `statsmatrix' = reghdfe_statsmatrix[1..., "`name'"]
			}
			else if (`is_match') {
				matrix `statsmatrix' = `statsmatrix' , reghdfe_statsmatrix[1..., "`name'"]	
			}
		}
	}
end
