program reghdfe_estat, rclass
	version `=cond(c(version)<14, c(version), 13)'
	if !inlist("`e(cmd)'", "reghdfe", "ppmlhdfe") {
		error 301
	}
	
	gettoken key 0 : 0, parse(", ")
	local lkey = length(`"`key'"')

	if `"`key'"' == substr("summarize",1,max(2,`lkey')) {

		local 0 `rest'
		syntax [anything] , [*] [noheader] // -noheader- gets silently ignored b/c it will always be -on-

		**if ("`anything'"=="") {
		**	* By default include the instruments
		**	local anything // `e(depvar)' `e(indepvars)' `e(endogvars)' `e(instruments)'
		**}

		* Need to use -noheader- as a workaround to the bug in -estat_summ-
		estat_summ `anything' , `options' noheader

	}
	else if `"`key'"' == "vce" {
		vce `0'
	}
	else if `"`key'"' == "ic" {
		syntax, [*]
		estat_default ic, df(`=e(df_m)+1') `options'
	}
	else {
		di as error `"invalid subcommand `key'"'
		exit 321
	}
	return add // ?
end
