//------------------------------------------------------------------------------
// Expand Factor Variables, interactions, and time-series vars
//------------------------------------------------------------------------------
// This basically wraps -fvrevar-, adds labels, and drops omitted/base
cap pr drop ExpandFactorVariables
program define ExpandFactorVariables, rclass
syntax varlist(min=1 numeric fv ts) [if] [,setname(string)] [CACHE]

	local expanded_msg `"" - variable expansion for `setname': " as result "`varlist'" as text " ->""'

	* It's (usually) a waste to add base and omitted categories
	* EG: if we use i.foreign#i.rep78 , several categories will be redundant, seen as e.g. "0b.foreign" in -char list-
	* We'll also exclude base categories that don't have the "bn" option (to have no base)

	* However, if there is a cont. interaction, then we DO want the base categories!

	* Loop for each var and then expand them into i.var -> 1.var.. and loop
	* Why two loops? B/c I want to save each var expansion to allow for a cache

	if ("`cache'"!="") mata: varlist_cache = asarray_create()

	local newvarlist
	* I can't do a simple foreach!
	* Because a factor expression could be io(3 4).rep78
	* and foreach would split the parens in half
	while (1) {
	gettoken fvvar varlist : varlist, bind
	if ("`fvvar'"=="") continue, break

		fvrevar `fvvar' `if' // , stub(__V__) // stub doesn't work in 11.2
		local contents

		foreach var of varlist `r(varlist)' {
			
			* Get readable varname
			local fvchar : char `var'[fvrevar]
			local tschar : char `var'[tsrevar]
			local name `fvchar'`tschar'
			local color input
			if ("`name'"=="") {
				local name `var'
				local color result
			}
			char `var'[name] `name'
			la var `var' "[Tempvar] `name'"

			* See if the factor can be dropped safely
			if (substr("`var'", 1, 2)=="__") {
				local color result
				local parts : subinstr local fvchar "#" " ", all
				local continteraction = strpos("`fvchar'", "c.")>0
				foreach part of local parts {
					*di as error "part=<`part'> cont=`continteraction' all=<`fvchar'>"
					* "^[0-9]+b\." -> "b.*\."
					if (regexm("`part'", "b.*\.") & !`continteraction') | regexm("`part'", "o.*\.") {
						local color error	
					}
				}
				if ("`color'"=="error") {
					local color result
				}


				* Need to rename it, or else it gets dropped since its a tempvar
				if ("`color'"!="error") {
					local newvarbase : subinstr local name "." "__", all // pray that no variable has three _
					local newvarbase : subinstr local newvarbase "#" "_X_", all // idem
					local newvarbase : permname __`newvarbase', length(30)

					* In what cases will just using newvarbase fail???
					local i // Empty
					while (1) {
						local newvar "`newvarbase'`i'"
					
						if ("`i'"=="") {
							local i 1
						}
						else {
							local ++i
						}

						Assert `i'<1000, msg("Couldn't create tempvar for `var' (`name')")
						cap conf new var `newvar', exact
						if _rc==0 {
							continue, break
						}
					}

					rename `var' `newvar'
					local var `newvar'
				}
			}

			* Save contents of the expansion for optional -cache-			
			if ("`color'"!="error") {
				local contents `contents' `var'
			}
			
			* Set debug message
			local expanded_msg `"`expanded_msg' as `color' " `name'" as text " (`var')""'
		}

		if ("`cache'"!="") mata: asarray(varlist_cache, "`fvvar'", "`contents'")
		Assert "`contents'"!="", msg("error: variable -`fvvar'- in varlist -`varlist'- in category -`setname'- is  empty after factor expansion")
		local newvarlist `newvarlist' `contents'
	}

	* Yellow=Already existed, White=Created, Red=NotCreated (omitted or base)
	Debug, level(3) msg(`expanded_msg')
	return local varlist "`newvarlist'"
end
