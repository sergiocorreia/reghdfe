cap pr drop ParseOneAbsvar
program define ParseOneAbsvar, rclass
	syntax, ABSVAR(string)

	Assert !strpos("`absvar'","###"), msg("error parsing <`absvar'> : ### is invalid")
	Assert regexm("`absvar'", "^[a-zA-Z0-9_=.#]+$"), msg("error parsing <`absvar'> : illegal characters ")
	Assert !regexm("`absvar'", "##([^c]|(c[^.]))"), msg("error parsing <`absvar'> : expected c. after ##")
	local original_absvar `absvar'

* Split at equal sign
	local equalsign = strpos("`absvar'","=")
	local target = substr("`absvar'",1,`equalsign'-1)
	local absvar = substr("`absvar'",`equalsign'+1, .)
	if ("`target'"!="") conf new var `target'

	local is_interaction = strpos("`absvar'", "#")>0
	local is_bivariate = strpos("`absvar'", "##")>0

* Split interactions
	mata: st_local("vars", subinstr("`absvar'", "#", " ") )
	foreach var of local vars {

		local dot = strpos("`var'", ".")
		local root = substr("`var'", `dot'+1, .)
		unab root : `root' , max(1)
		conf numeric var `root'
		
		local prefix = substr("`var'", 1, `dot'-1)
		local prefix = lower( cond("`prefix'"=="", "i", "`prefix'") ) // -i.- is default prefix

		Assert inlist("`prefix'", "i", "c") , msg("error parsing <`absvar'><`var'> : only i. and c. allowed, not `prefix'.")
		Assert !strpos("`root'", ".") , msg("error parsing <`absvar'><`var'> : no time series operators allowed")
		
		if ("`prefix'"=="i") {
			local ivars `ivars' `root'
		}
		else {
			Assert "`cvars'"=="", msg("error: can't have more than one continuous variable in the interaction")
			local cvars `cvars' `root'
		}
	}
	local tab  "        "
	Debug, level(3) msg(as text "    Parsing " as result "`original_absvar'")
	Debug, level(3) msg(as text "`tab'ivars = " as result "`ivars'")
	if ("`cvars'"!="") Debug, level(3) msg(as text "`tab'cvars = " as result "`cvars'")
	if ("`target'"!="") Debug, level(3) msg(as text "`tab'target = " as result "`target'")
	Debug, level(3) msg(as text "`tab'is_interaction = " as result "`is_interaction'")
	Debug, level(3) msg(as text "`tab'is_bivariate = " as result "`is_bivariate'")
	// Debug, level(3) newline

	return scalar is_interaction = `is_interaction'
	return scalar is_cont_interaction = `is_interaction' & ("`cvars'"!="")
	return scalar is_bivariate = `is_bivariate'
	if ("`target'"!="") return local target "`target'"
	if ("`cvars'"!="") return local cvars "`cvars'"
	return local ivars "`ivars'"
end
