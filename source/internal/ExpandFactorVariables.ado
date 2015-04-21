// -------------------------------------------------------------------------------------------------
// Expand factor time-series variables
// -------------------------------------------------------------------------------------------------
* Steps:
* 1) Call -fvrevar-
* 2) Label newly generated temporary variables
* 3) Drop i) omitted variables, and ii) base variables (if not part of a #c.var interaction)

capture program drop ExpandFactorVariables
program define ExpandFactorVariables, rclass
syntax varlist(min=1 numeric fv ts) [if] [,setname(string)] [SAVECACHE(integer 0)] verbose(integer)
	
	* If saving the data for later regressions -savecache(..)- we will need to match each expansion to its newvars
	* The mata array is used for that
	* Note: This explains why we need to wrap -fvrevar- in a loop
	if (`savecache') {
		mata: varlist_cache = asarray_create()
		mata: asarray_notfound(varlist_cache, "")
	}

	local expanded_msg `"" - variable expansion for `setname': {res}`varlist'{txt} ->""'
	while (1) {
		gettoken factorvar varlist : varlist, bind
		if ("`factorvar'"=="") continue, break

		fvrevar `factorvar' `if' // , stub(__V__) // stub doesn't work in Stata 11.2
		local contents
		foreach var of varlist `r(varlist)' {
			LabelRenameVariable `var' // Tempvars not renamed will be dropped automatically
			if !r(is_dropped) {
				local contents `contents' `r(varname)'
				// if (`savecache') di as error `"<mata: asarray(varlist_cache, "`factorvar'", "`r(varname)'")>"'
				if (`savecache') mata: asarray(varlist_cache, "`factorvar'", "`r(varname)'")
			}
			* Yellow=Already existed, White=Created, Red=NotCreated (omitted or base)
			local color = cond(r(is_dropped), "error", cond(r(is_newvar), "input", "result"))
			if (`verbose'>3) {
				local expanded_msg `"`expanded_msg' as `color' " `r(name)'" as text " (`r(varname)')""'
			}
		}
		Assert "`contents'"!="", msg("error: variable -`factorvar'- in varlist -`varlist'- in category -`setname'- is  empty after factor/time expansion")
		local newvarlist `newvarlist' `contents'
	}

	Debug, level(4) msg(`expanded_msg')
	return clear
	return local varlist "`newvarlist'"
end

capture program drop LabelRenameVariable
program define LabelRenameVariable, rclass
syntax varname
	local var `varlist'
	local fvchar : char `var'[fvrevar]
	local tschar : char `var'[tsrevar]
	local is_newvar = ("`fvchar'`tschar'"!="") & substr("`var'", 1, 2)=="__"
	local name "`var'"
	local will_drop 0

	if (`is_newvar') {
		local name "`fvchar'`tschar'"
		local parts : subinstr local fvchar "#" " ", all
		local has_cont_interaction = strpos("`fvchar'", "c.")>0
		local is_omitted 0
		local is_base 0
		foreach part of local parts {
			if (regexm("`part'", "b.*\.")) local is_base 1
			if (regexm("`part'", "o.*\.")) local is_omitted 1
		}

		local will_drop = (`is_omitted') | (`is_base' & !`has_cont_interaction')
		if (!`will_drop') {
			char `var'[name] `name'
			la var `var' "[TEMPVAR] `name'"
			local newvar : subinstr local name "." "__", all
			local newvar : subinstr local newvar "#" "_X_", all
			* -permname- selects newname# if newname is taken (# is the first number available)
			local newvar : permname __`newvar', length(30)
			rename `var' `newvar'
			local var `newvar'
		}
	}

	return scalar is_newvar = `is_newvar'
	return scalar is_dropped = `will_drop'
	return local varname "`var'"
	return local name "`name'"
end
