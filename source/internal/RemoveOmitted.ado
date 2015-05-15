* Remove omitted variables from a beta matrix, and return remaining indepvars
capture program drop RemoveOmitted
program define RemoveOmitted, rclass
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
