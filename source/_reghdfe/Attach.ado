capture program drop Attach
program define Attach, eclass
	syntax, [NOTES(string)] [statsmatrix(string)] summarize_quietly(integer)
	
	* Summarize
	* This needs to happen after all the missing obs. have been dropped and the only obs. are those that *WILL* be in the regression
	if ("`statsmatrix'"!="") {
		* Update beta vector
		* ...

		if (!`summarize_quietly' & "`statsmatrix'"!="") {
			di as text _n "{ul:Regression Summary Statistics}" _c
			tempname stats_transpose
			matrix `stats_transpose' = (`statsmatrix')'
			matrix list `stats_transpose', noheader
		}
		ereturn matrix summarize = `statsmatrix'
	}

	* Parse key=value options and append to ereturn as hidden
	mata: st_local("notes", strtrim(`"`notes'"')) // trim (supports large strings)
	local keys
	while (`"`notes'"'!="") {
		gettoken key notes : notes, parse(" =")
		Assert !inlist("`key'","sample","time"), msg("Key cannot be -sample- or -time-") // Else -estimates- will fail
		gettoken _ notes : notes, parse("=")
		gettoken value notes : notes, quotes
		local keys `keys' `key'
		ereturn hidden local `key' `value'
	}
	if ("`keys'"!="") ereturn hidden local keys `keys'

end
