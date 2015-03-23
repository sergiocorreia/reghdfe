cap pr drop Save
program define Save, rclass
	* Run this after -Demean .. , save_fe(1)-
	* For each FE, if it has a -target-, add label, chars, and demean or divide
	CheckCorrectOrder save
	syntax , original_depvar(string)

	mata: st_local("G", strofreal(G))
	mata: st_local("weightexp", weightexp)
	forv g=1/`G' {

		// ivars cvars target varname varlabel is_interaction is_cont_interaction is_bivariate is_mock levels
		mata: fe2local(`g')
		if ("`target'"=="") continue

		* Rename, add label and chars
		rename __Z`g'__ `target'
		local label `varlabel'
		la var `target' "Fixed effects of `label' on `original_depvar'"
		char `target'[label] `label'
		char `target'[levels] `levels'

		* Substract mean, or divide by cvar (fixing division by zero errors)
		if ("`cvars'"!="" & !(`is_bivariate' & !`is_mock')) {
			char `target'[cvars] `cvars'
			qui replace `target' = cond(abs(`cvars')<epsfloat(), 0,  `target'/`cvars')
			// BUGBUG BUGBUG float(`target'/`cvars')) -> this makes them have the same FE but loses precision!
		}
		else {
			qui su `target' `weightexp', mean
			qui replace `target' = `target' - r(mean)
			// BUGBUG BUGBUG -> WHAT WAS THE PROBLEM WITH THIS?
		}

		local keepvars `keepvars' `target'
	}

	cap drop __Z*__
	return local keepvars " `keepvars'" // the space prevents MVs
end

