//------------------------------------------------------------------------------
// Name tempvars into e.g. L.x i1.y i2.y AvgE:z , etc.
//------------------------------------------------------------------------------
cap pr drop FixVarnames
program define FixVarnames, rclass
local vars `0'

	foreach var of local vars {
		local newname `var'

		* -var- can be <o.__W1__>
		if ("`var'"=="_cons") {
			local newname `var'
		}
		else {
			fvrevar `var', list
			local basevar "`r(varlist)'"
			local label : var label `basevar'
			local is_temp = substr("`basevar'",1,2)=="__"
			local is_omitted = strpos("`var'", "o.")
			local prefix = cond(`is_omitted'>0, "o.", "")
			local name : char `basevar'[name]

			 if (`is_temp' & "`name'"!="") {
				local newname `prefix'`name'
				
				* Fix bug when the var is omitted:
				local bugmatch = regexm("`newname'", "^o\.([0-9]+)b?\.(.+)$")
				if (`bugmatch') {
					local newname = regexs(1) + "o." + regexs(2) // EG: 1o.var
				}
			}

		}
		
		Assert ("`newname'"!=""), msg("var=<`var'> --> new=<`newname'>")
		local newnames `newnames' `newname'
	}

	local A : word count `vars'
	local B : word count `newnames'
	Assert `A'==`B', msg("`A' vars but `B' newnames")
	
	***di as error "newnames=`newnames'"
	return local newnames "`newnames'"
end
