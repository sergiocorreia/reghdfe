*! reghdfe 3.0.3 12may2015
*! Sergio Correia (sergio.correia@duke.edu)

*! version 1.1.0 10jul2014
* TODO: Not tested for -avge- variables!
program reghdfe_estat, rclass
	if (c(version)<12) di as error "(warning: -reghdfe- has only been throughly tested with Stata 12 and above)"
	version `=cond(c(version)<14, c(version), 13)'
	if "`e(cmd)'" != "reghdfe" {
		error 301
	}
	
	gettoken key 0 : 0, parse(", ")
	local lkey = length(`"`key'"')

	if `"`key'"' == substr("summarize",1,max(2,`lkey')) {

		local 0 `rest'
		syntax [anything] , [*] [noheader] // -noheader- gets silently ignored b/c it will always be -on-

		if ("`anything'"=="") {
			* By default include the instruments
			local anything `e(depvar)' `e(indepvars)' `e(endogvars)' `e(instruments)'
		}

		* Need to use -noheader- as a workaround to the bug in -estat_summ-
		estat_summ `anything' , `options' noheader

	}
	else if `"`key'"' == "vce" {
		vce `0'
	}
	else {
		di as error `"invalid subcommand `key'"'
		exit 321
	}
	return add // ?
end
