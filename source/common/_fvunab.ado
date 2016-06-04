cap pr drop _fvunab
pr _fvunab, sclass
	* Variant of -fvunab- that does not expand x##y
	loc remainder `0'
	
	while ("`remainder'" != "") {
		gettoken 0 remainder: remainder, parse(" #.()=")
		loc next_char = substr("`remainder'", 1, 1)
		if !inlist("`0'", "#", ".", "(", ")", "=") & ("`next_char'" != ".") {
			syntax varlist(numeric fv ts)
			loc 0 `varlist'
		}
		if ("`next_char'" != " ") loc next_char
		loc answer "`answer'`0'`next_char'"
	}
	sreturn local varlist `answer'
end

/* USAGE:
sysuse auto
_fvunab F2.pri   tu##c.L.trun#ibn.foreign (pri	= tu##for#c.pri) weigh
di "`s(varlist)'"
*/
