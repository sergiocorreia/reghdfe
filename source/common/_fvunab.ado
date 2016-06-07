cap pr drop _fvunab
pr _fvunab, sclass
	sreturn clear
	syntax anything(name=remainder equalok) [, NOIsily]
	loc is_numlist 0 // Will match inside L(1 2 3) or L(-1/1)
	
	while ("`remainder'" != "") {
		* gettoken won't place spaces in 0;
		* but we can see if a space is coming with `next_char'
		gettoken 0 remainder: remainder, parse(" #.()=")
		loc next_char = substr("`remainder'", 1, 1)
		
		* Match common delimiters
		loc delim1 = inlist("`0'", "#", ".", "(", ")", "=")

		* Match "i" and "L" in "i.turn L(1 2)"
		loc delim2 = inlist("`next_char'", ".", "(")
		
		if !(`delim1' | `delim2' | `is_numlist') {
			syntax varlist(numeric fv ts)
			loc 0 `varlist'
			loc unique `unique' `varlist'
		}
		if ("`0'" == ")") loc is_numlist 0
		if ("`0'" != "." & "`next_char'" == "(") loc is_numlist 1
		
		if ("`next_char'" != " ") loc next_char // add back spaces
		loc answer "`answer'`0'`next_char'"
		if ("`noisily'" != "") di as result "{bf:`answer'}"
	}
	local unique : list uniq unique
	sreturn local basevars `unique' // similar to fvrevar,list
	sreturn local varlist `answer'
end

/* _FVUNAB

Description:
	Variant of -fvunab- that does not expand "x##y" into "x y x#y"
	Also does not expand "x#y" into "i.x#i.y"

Example:
	sysuse auto
	_fvunab F2.pri   tu##c.L.trun#ibn.foreign (pri	= tu##for#c.pri) weigh
	di "`s(varlist)'"
*/
