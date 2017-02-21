* Convenience function used b/c triM_cons renames e(b) into e(b_trim)

cap pr drop notrim
pr notrim, eclass
args size adjdof
	matrix trim_b = e(b)
	matrix trim_V = e(V)
	if ("`adjdof'"!="") {
		matrix trim_V = trim_V * `adjdof'
		// ereturn scalar F = e(F) / `adjdof'
	}
	ereturn matrix trim_b = trim_b
	ereturn matrix trim_V = trim_V
end
