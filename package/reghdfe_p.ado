*! version 1.1.0 10jul2014
* predict after reghdfe
* TODO: Not tested for -avge- variables!

program define reghdfe_p
	local version `clip(`c(version)', 11.2, 13.1)' // 11.2 minimum, 13+ preferred
	qui version `version'
	
	*if "`e(cmd)'" != "reghdfe" {
	*	error 301
	*}
	syntax anything [if] [in] , [XB XBD D Residuals SCores]
	if (`"`scores'"' != "") {
		_score_spec	`anything'
		local varlist `s(varlist)'
	}
	else {
		local 0 `anything'
		syntax newvarname  // [if] [in] , [XB XBD D Residuals SCores]
	}

	local weight "[`e(wtype)'`e(wexp)']" // After -syntax-!!!
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
		
		* Only estimate using e(sample) except when computing xb (when we don't need -d- and can predict out-of-sample)
		if (`"`if'"'!="") {
			local if `if' & e(sample)==1
		}
		else {
			local if "if e(sample)==1"
		}

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
				qui replace `d' = `d' + `e(hdfe_target`g')' `if' `in'
			}
			else {
				qui replace `d' = `d' + `e(hdfe_target`g')' * `e(hdfe_cvar`g')' `if' `in'
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
			qui replace `d' = `d' + `e(avge_target`g')' `if' `in'
		}
	} // Finished creating `d' if needed

	tempvar xb // XB will eventually contain XBD and RESID if that's the output
	_predict double `xb' `if' `in', xb

	if ("`option'"=="xb") {
		rename `xb' `varlist'
	}
	else {
		* Make residual have mean zero (and add that to -d-)
		su `e(depvar)' `if' `in' `weight', mean
		local mean = r(mean)
		su `xb' `if' `in' `weight', mean
		local mean = `mean' - r(mean)
		su `d' `if' `in' `weight', mean
		local mean = `mean' - r(mean)
		qui replace `d' = `d' + `mean' `if' `in'

		if ("`option'"=="d") {
			rename `d' `varlist'
			la var `varlist' "d[`fixed_effects']"
		}
		else if ("`option'"=="xbd") {
			qui replace `xb' = `xb' + `d' `if' `in'
			rename `xb' `varlist'
			la var `varlist' "Xb + d[`fixed_effects']"
		}
		else if ("`option'"=="residuals") {
			qui replace `xb' = `e(depvar)' - `xb' - `d' `if' `in'
			rename `xb' `varlist'
			la var `varlist' "Residuals"
		}
		else {
			error 112
		}
	}

	fvrevar `e(depvar)', list
	local format : format `r(varlist)'
	format `format' `varlist'
end
