program define reghdfe_p, rclass
	* Note: we IGNORE typlist and generate the newvar as double
	* Note: e(resid) is missing outside of e(sample), so we don't need to condition on e(sample)

	* HACK: Intersect -score- and replace with -residuals-
	cap syntax anything [if] [in], SCore
	loc was_score = !c(rc)
	if (`was_score') {
		* Call _score_spec to get newvarname; discard type
		* - This resolves wildcards that -margins- sends to predict (e.g. var* -> var1)
		* - Do we really need to pass it `if' and `in' ?
		_score_spec `anything', score
		loc 0 `s(varlist)' `if' `in' , residuals
	}

	syntax newvarname [if] [in] [, XB STDP Residuals D XBD DResiduals]

	* Ensure there is only one option
	opts_exclusive "`xb' `stdp' `residuals' `d' `xbd' `dresiduals'"

	* Default option is xb
	cap opts_exclusive "`xb' `stdp' `residuals' `d' `xbd' `dresiduals' placeholder"
	if (!c(rc)) {
		di as text "(option xb assumed; fitted values)"
		loc xb "xb"
	}

	local fixed_effects "`e(absvars)'"

	* Except for xb and stdp, we need the previously computed residuals
	if ("`xb'" == "" & "`stdp'" == "") {
		_assert ("`e(resid)'" != ""), msg("you must add the {bf:resid} option to reghdfe before running this prediction")
		conf numeric var `e(resid)', exact
	}

	if ("`xb'" != "" | "`stdp'" != "") {
		* xb: normal treatment
		PredictXB `varlist' `if' `in', `xb' `stdp'
	}
	else if ("`residuals'" != "") {
		* resid: just return the preexisting variable
		gen double `varlist' = `e(resid)' `if' `in'
		la var `varlist' "Residuals"
		if (`was_score') return local scorevars `varlist'
	}
	else if ("`d'" != "") {
		* d: y - xb - resid
		tempvar xb
		PredictXB `xb' `if' `in', xb
		gen double `varlist' = `e(depvar)' -  `xb' - `e(resid)' `if' `in'
		la var `varlist' "d[`fixed_effects']"
	}
	else if ("`xbd'" != "") {
		* xbd: y - resid
		gen double `varlist' = `e(depvar)' - `e(resid)' `if' `in'
		la var `varlist' "Xb + d[`fixed_effects']"
	}
	else if ("`dresiduals'" != "") {
		* dresid:	y - xb
		tempvar xb
		PredictXB `xb' `if' `in', xb
		gen double `varlist' = `e(depvar)' -  `xb' `if' `in'
	}
	else {
		error 100
	}
end

program PredictXB
	syntax newvarname [if] [in], [*]
	cap matrix list e(b) // if there are no regressors, _predict fails
	if (c(rc)) {
		gen double `varlist' = 0 `if' `in'
	}
	else {
		_predict double `varlist' `if' `in', `options'
	}
end
