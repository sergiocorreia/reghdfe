//------------------------------------------------------------------------------
// Name tempvars into e.g. L.x i1.y i2.y AvgE:z , etc.
//------------------------------------------------------------------------------
cap pr drop FixVarnames
program define FixVarnames, rclass
local vars `0'

	foreach var of local vars {
		* Note: -var- can be <o.var>
		_ms_parse_parts `var'
		local is_omitted = r(omit)
		local name = r(name)

		local is_temp = substr("`name'",1,2)=="__"
		local newname : char `name'[name]
		*local label : var label `basevar'

		* Stata requires all interaction elements to have an o.
		if (`is_omitted' & `is_temp') {
			while regexm("`newname'", "^(.*[^o])\.(.*)$") {
				local newname = regexs(1) + "o." + regexs(2)
			}
		}
		else if (`is_omitted') {
			local newname "o.`name'" // same as initial `var'!
		}

		Assert ("`newname'"!=""), msg("var=<`var'> --> new=<`newname'>")
		local newnames `newnames' `newname'
	}

	local A : word count `vars'
	local B : word count `newnames'
	Assert `A'==`B', msg("`A' vars but `B' newnames")
	return local newnames "`newnames'"
end
