*! version 1.1.0 10jul2014
* predict after reghdfe
* TODO: Not tested for -avge- variables!

program define reghdfe_p, // sortpreserve properties(default_xb)
	local version `clip(`c(version)', 11.2, 13.1)' // 11.2 minimum, 13+ preferred
	qui version `version'
	
	if "`e(cmd)'" != "reghdfe" {
		error 301
	}

	syntax newvarname [if] [in] , [XB XBD D Residuals SCores]

	local option `xb' `xbd' `d' `residuals' `scores'
	if ("`option'"=="") local option xb // The default, as in -areg-
	local numoptions : word count `option'
	if (`numoptions'!=1) {
		di as error "(predict reghdfe) syntax error; specify one and only one option"
		exit 112
	}
	if ("`option'"=="scores") local option residuals

	local fixed_effects = e(absvars)

	* We need to have saved FEs and AvgEs for every option except -xb-
	if ("`option'"!="xb") {
		
		* Construct -d- if needed (sum of FEs)
		tempvar d
		qui gen double `d' = 0 `if' `in'

		forv g=1/`e(N_hdfe)' {
			local ok 0
			cap conf e `e(hdfe_target`g')'
			if (_rc==0) {
				cap conf numeric var `e(hdfe_target`g')'
			}
			if (_rc!=0) {
				di as error "In order to predict, all the FEs need to be saved with the absorb option (#`g' was not)" _n "For instance, instead of {it:absorb(i.year i.firm)}, set absorb(FE_YEAR=i.year FE_FIRM=i.firm)"
				exit 112
			}

			if missing("`e(hdfe_cvar`g')'") {
				qui replace `d' = `d' + `e(hdfe_target`g')'
			}
			else {
				qui replace `d' = `d' + `e(hdfe_target`g')' * `e(hdfe_cvar`g')'
			}
		}
		local K = cond( e(N_avge)==. , 0 , e(N_avge) )
		forv g=1/`K' {
			local ok 0
			cap conf e `e(avge_target`g')'
			if (_rc==0) {
				cap conf numeric var `e(avge_target`g')'
			}
			if (_rc!=0) {
				di as error "(predict reghdfe) you need to save all AvgEs in reghdfe, AvgE`g' not saved"
				exit 112
			}
			qui replace `d' = `d' + `e(avge_target`g')'
		}
	} // Finished creating `d' if needed
	
	
	* Construct -xb- if needed
	if ("`option'"!="d") {
		tempvar xb
		_predict double `xb' `if' `in', xb
	}
	
	* Adjusting for -noconstant- option
	local adj_cons = cond(e(_cons)<., e(_cons), 0)
	cap replace `xb' = `xb' + `adj_cons'
	cap replace `d' = `d' - `adj_cons'
	
	if ("`option'"=="xb") {
		rename `xb' `varlist'
	}
	else if ("`option'"=="d") {
		rename `d' `varlist'
	}
	else {
		qui replace `xb' = `xb' + `d' `if' `in'
		if ("`option'"=="xbd") {
			rename `xb' `varlist'
		}
		else if ("`option'"=="residuals") {
			qui replace `xb' = `e(depvar)' - `xb' `if' `in'
			rename `xb' `varlist'
		}
		else {
			error 112
		}
	}

	fvrevar `e(depvar)', list
	local format : format `r(varlist)'
	format `format' `varlist'
	* TODO: Allow [type] and recast to that, as -predict- usually does
end
