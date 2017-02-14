program define reghdfe_p
	* Note: we IGNORE typlist and generate the newvar as double
	* Note: e(resid) is missing outside of e(sample), so we don't need to condition on e(sample)

	loc opt_list XB STDP Residuals D XBD DResiduals
	syntax newvarname [if] [in] [, XB STDP Residuals D XBD DResiduals]

	* Ensure there is only one option
	opts_exclusive "`xb' `stdp' `residuals' `d' `xbd' `dresiduals'"

	* Default option is xb
	cap opts_exclusive "`xb' `stdp' `residuals' `d' `xbd' `dresiduals' foobar"
	if (!c(rc)) {
		di as text "(option xb assumed; fitted values)"
		loc xb "xb"
	}

	* Except for xb and stdp, we need the previously computed residuals
	if ("`xb'" == "" & "`stdp'" == "") {
		_assert ("`e(resid)'" != ""), msg("you must add the {bf:resid} option to reghdfe before running this prediction")
		conf numeric var `e(resid)', exact
	}

	if ("`xb'" != "" | "`stdp'" != "") {
		* xb: normal treatment
		_predict double `varlist' `if' `in', `xb' `stdp'
	}
	else if ("`residuals'" != "") {
		* resid: just return the preexisting variable
		gen double `varlist' = `e(resid)' `if' `in'
	}
	else if ("`d'" != "") {
		* d: y - xb - resid
		tempvar xb
		_predict double `xb' `if' `in', xb
		gen double `varlist' = `e(depvar)' -  `xb' - `e(resid)' `if' `in'
	}
	else if ("`xbd'" != "") {
		* xbd: y - resid
		gen double `varlist' = `e(depvar)' - `e(resid)' `if' `in'
	}
	else if ("`dresiduals'" != "") {
		* dresid:	y - xb
		tempvar xb
		_predict double `xb' `if' `in', xb
		gen double `varlist' = `e(depvar)' -  `xb' `if' `in'
	}
	else {
		error 100
	}
end
