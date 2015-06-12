// -------------------------------------------------------------
// Display Regression Table
// -------------------------------------------------------------
cap pr drop Replay
 program define Replay, eclass
	syntax , [*]
	Assert e(cmd)=="reghdfe"
	local subcmd = e(subcmd)
	Assert "`subcmd'"!="" , msg("e(subcmd) is empty")
	if (`c(version)'>=12) local hidden hidden

	if ("`e(stored_estimates)'"!="" & "`stage'"=="iv") {
		local est_list = e(stored_estimates)
		tempname hold
		estimates store `hold'
		foreach est of local est_list {
			qui estimates restore `est'
			Replay			
		}
		ereturn clear // Need this because -estimates restore- behaves oddly
		qui estimates restore `hold'
		assert e(cmd)=="reghdfe"
		estimates drop `hold'
	}


	* conf matrix e(first)


	local diopts = "`e(diopts)'"
	if ("`options'"!="") { // Override
		_get_diopts diopts /* options */, `options'
	}

	if ("`subcmd'"=="ivregress") {
		* Don't want to display anova table or footnote
		_coef_table_header
		_coef_table, `diopts' bmatrix(`b') vmatrix(e(V)) // plus 
	}
	else if ("`subcmd'"=="ivreg2") {
		* Backup before showing both first and second stage
		tempname hold


		// BUGBUG: Update this part
		estimates store `hold'
		// ereturn repost b=`b', rename
		ereturn local cmd = "`subcmd'"
		`subcmd' , `diopts'
		ereturn clear // Need this because -estimates restore- behaves oddly
		qui estimates restore `hold'
		assert e(cmd)=="reghdfe"
		estimates drop `hold'

		*ereturn local cmd = "reghdfe"
		*matrix `b' = e(b)
		*matrix colnames `b' = `backup_colnames'
		*ereturn repost b=`b', rename
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
		local plus = cond(e(model)=="ols" & inlist("`e(vce)'", "unadjusted", "ols"), "plus", "")
		_coef_table, `plus' `diopts' bmatrix(`b') vmatrix(e(V))
	}
	mata: reghdfe_width = max(strlen(st_matrixcolstripe_split("r(table)", 32, 0)))
	mata: st_local("width" , strofreal(reghdfe_width))
	mata: mata drop reghdfe_width
	if (`width'<12) local width 12
	ereturn `hidden' scalar width = `width'
	reghdfe_footnote
end
