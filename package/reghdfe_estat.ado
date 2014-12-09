*! version 1.1.0 10jul2014
* TODO: Not tested for -avge- variables!
program reghdfe_estat, rclass
	if (c(version)<12) di as error "(warning: -reghdfe- has only been throughly tested with Stata 12 and above)"
	version `=cond(c(version)<14, c(version), 13)'
	if "`e(cmd)'" != "reghdfe" {
		error 301
	}
	
	gettoken key rest : 0, parse(", ")
	local lkey = length(`"`key'"')

	if `"`key'"' == substr("summarize",1,max(2,`lkey')) {
		* estat_summ `rest'
		* There is a bug in estat_summ , so we'll do it by hand
		if ("`e(wtype)'"!="") local weight [`e(wtype)'`e(wexp)']
		su `e(depvar)' `e(indepvars)' `e()' `e(endogvars)' `e(instruments)' if e(sample) `weight'
	}
	else if `"`key'"' == "vce" {
		vce `rest'
	}
	return add // ?
end
