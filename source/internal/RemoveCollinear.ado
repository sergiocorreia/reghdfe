* Tag Collinear Variables with an o. and compute correct e(df_m)
	* Obtain K so we can obtain DoF = N - K - kk
	* This is already done by regress EXCEPT when clustering
	* (but we still need the unclustered version for r2_a, etc.)

capture program drop RemoveCollinear
program define RemoveCollinear, rclass
	syntax, depvar(varname numeric) [indepvars(varlist numeric) weightexp(string)]

	qui _rmcoll `indepvars' `weightexp', forcedrop
	local okvars = r(varlist)
	if ("`okvars'"==".") local okvars
	local df_m : list sizeof okvars

	foreach var of local indepvars {
		local ok : list var in okvars
		local prefix = cond(`ok', "", "o.")
		local label : char `var'[name]
		if (!`ok') di as text "note: `label' omitted because of collinearity"
		local varlist `varlist' `prefix'`var'
	}

	mata: st_local("vars", strtrim(stritrim( "`depvar' `varlist'" )) ) // Just for aesthetic purposes
	return local vars "`vars'"
	return scalar df_m = `df_m'

end
