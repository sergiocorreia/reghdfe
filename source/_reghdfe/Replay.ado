// -------------------------------------------------------------
// Display Regression Table
// -------------------------------------------------------------
cap pr drop Replay
 program define Replay, eclass
	syntax , [*]
	Assert e(cmd)=="reghdfe"
	local subcmd = e(subcmd)
	Assert "`subcmd'"!="" , msg("e(subcmd) is empty")

	* Add pretty names for AvgE variables
	tempname b
	matrix `b' = e(b)
	local backup_colnames : colnames `b'
	matrix colnames `b' = `e(prettynames)'
	local savefirst = e(savefirst)
	local suboptions = e(suboptions)

	di as error "We need to add something like ivreg2 for the clusters.."
	di as error "We can use the LHS of header: Number of clusters (turn) = as result 18"
	
	di as error "Same for explaining unusual SEs with -avar- and -mwc-"
	di as error "like <Statistics robust to heteroskedasticity and clustering on turn and t>"
	di as error "like <and kernel-robust to common correlated disturbances (Driscoll-Kraay)>"
	di as error "<  kernel=Bartlett; bandwidth=2>"
	di as error "<  time variable (t):  t>"
	di as error "<  group variable (i): turn>"

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
		
		if ("`e(ivreg2_firsteqs)'"!="") {
			estimates store `hold'

			local i 0
			foreach fs_eqn in `e(ivreg2_firsteqs)' {
				local instrument  : word `++i' of `e(instd)'
				di as input _n "{title:First stage for `instrument'}"
				estimates replay `fs_eqn' , nohead `diopts'
				if (!`savefirst') estimates drop `fs_eqn'
			}

			ereturn clear
			qui estimates restore `hold'
			di as input _n "{title:Second stage}"
		}

		estimates store `hold'
		ereturn repost b=`b', rename
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
		_coef_table_header
		di
		local plus = cond(e(model)=="ols" & inlist("`e(vce)'", "unadjusted", "ols"), "plus", "")
		_coef_table, `plus' `diopts' bmatrix(`b') vmatrix(e(V))
	}

	reghdfe_footnote
	* Revert AvgE else -predict- and other commands will choke


end
