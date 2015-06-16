capture program drop GenUID
program define GenUID
	args uid
	local uid_type = cond(c(N)>c(maxlong), "double", "long")
	gen `uid_type' `uid' = _n // Useful for later merges
	la var `uid' "[UID]"
end
