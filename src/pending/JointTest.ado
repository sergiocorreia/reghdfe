* Compute model F-test; called by regress/mwc/avar wrappers
cap pr drop JointTest
pr JointTest, eclass
	args K
	if (`K'>0) {
		RemoveOmitted
		qui test `r(indepvars)' // Wald test
		if (r(drop)==1) {
			Debug, level(0) msg("WARNING: Missing F statistic (dropped variables due to collinearity or too few clusters).")
			ereturn scalar F = .
		}
		else {
			ereturn scalar F = r(F)
			if missing(e(F)) di as error "WARNING! Missing FStat"
		}
		ereturn scalar df_m = r(df)
		ereturn scalar rank = r(df) // Not adding constant anymore
	}
	else {
		ereturn scalar F = 0
		ereturn scalar df_m = 0
		ereturn scalar rank = 0 // Not adding constant anymore
	}
end

* Remove omitted variables from a beta matrix, and return remaining indepvars
cap pr drop RemoveOmitted
pr RemoveOmitted, rclass
	tempname b
	matrix `b' = e(b)
	local names : colnames `b'
	foreach name of local names {
		_ms_parse_parts `name'
		assert inlist(r(omit),0,1)
		if !r(omit) {
			local indepvars `indepvars' `name'
		}
	}
	return local indepvars `indepvars'
end
