// -------------------------------------------------------------
// Display Regression Table
// -------------------------------------------------------------
cap pr drop Replay
 program define Replay, eclass
	syntax , [stored] [*]
	Assert e(cmd)=="reghdfe"
	local subcmd = e(subcmd)
	Assert "`subcmd'"!="" , msg("e(subcmd) is empty")
	if (`c(version)'>=12) local hidden hidden

	if ("`stored'"!="" & "`e(stored_estimates)'"!="" & "`e(stage)'"=="iv") {
		local est_list "`e(stored_estimates)'"
		tempname hold
		estimates store `hold'
		foreach est of local est_list {
			cap estimates restore `est'
			if (!c(rc)) Replay
		}
		ereturn clear // Need this because -estimates restore- behaves oddly
		qui estimates restore `hold'
		assert e(cmd)=="reghdfe"
		estimates drop `hold'
	}

	if ("`e(stage)'"=="first") local first_depvar " - `e(depvar)'"
	if ("`e(stage)'"!="") di as text _n "{inp}{title:Stage: `e(stage)'`first_depvar'}"

	local diopts = "`e(diopts)'"
	if ("`options'"!="") { // Override
		_get_diopts diopts /* options */, `options'
	}

	if ("`subcmd'"=="ivregress") {
		* Don't want to display anova table or footnote
		_coef_table_header
		_coef_table, `diopts'
	}
	else if (strpos("`subcmd'","ivreg2")==1) {
		cap conf matrix e(first)
		if (c(rc)==0) local ffirst ffirst
		ereturn local cmd = "`subcmd'"
		`subcmd' , `diopts' `ffirst'
		ereturn local cmd = "reghdfe"
	}
	else {

		* Regress-specific code, because it doesn't play nice with ereturn
		sreturn clear 

		if "`e(prefix)'" != "" {
			_prefix_display, `diopts'
			exit
		}
		
		Header // _coef_table_header

		di
		local plus = cond("`e(model)'"=="ols" & inlist("`e(vce)'", "unadjusted", "ols"), "plus", "")
		_coef_table, `plus' `diopts'
	}
	reghdfe_footnote
end
