* Convenience function:
* "Trim <size>" will trim e(b) and e(V)

cap pr drop trim_coefs
pr trim_coefs, eclass
args size adjdof

	matrix trim_b = e(b)
	matrix trim_V = e(V)

	if ("`size'"=="") {
		loc size = colsof(trim_b) - 1
	}
	assert `size'>0

	matrix trim_b = trim_b[1, 1..`size']
	matrix trim_V = trim_V[1..`size', 1..`size']

	if ("`adjdof'"!="") {
		matrix trim_V = trim_V * `adjdof'
		ereturn scalar F = e(F) / `adjdof'
	}

	ereturn matrix trim_b = trim_b
	ereturn matrix trim_V = trim_V
end
