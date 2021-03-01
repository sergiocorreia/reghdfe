// -------------------------------------------------------------
// Display Regression Footnote
// -------------------------------------------------------------
program reghdfe5_footnote
syntax [, width(int 13)]

	if (`"`e(absvars)'"' == "_cons") {
		exit
	}

	tempname table
	matrix `table' = e(dof_table)
	mata: st_local("var_width", strofreal(max(strlen(st_matrixrowstripe("`table'")[., 2]))))
	if (`var_width' > `width') loc width = `var_width'
	loc rows = rowsof("`table'")
	loc cols = rowsof("`table'")
	local vars : rownames `table'

	// Setup table
	di as text _n "Absorbed degrees of freedom:"
	tempname mytab
	.`mytab' = ._tab.new, col(5) lmargin(0)
	.`mytab'.width	 `width'  | 12  	    12    		14 			1 |
	.`mytab'.pad		.		 1     		 1		   	 1			0
	.`mytab'.numfmt		.		%9.0g		%9.0g		%9.0g	  	.
	.`mytab'.numcolor	.		text 		text		result		.
	.`mytab'.sep, top

	local explain_exact 0
	local explain_nested 0
	
	// Header
	.`mytab'.titles "Absorbed FE" "Categories" " - Redundant" "  = Num. Coefs" ""
	.`mytab'.sep, middle

	// Body	
	forval i = 1/`rows' {
		local var : word `i' of `vars'
		loc var = subinstr("`var'", "1.", "", .)
		loc note " "
		if (`=`table'[`i', 4]'==1) {
			loc note "?"
			loc explain_exact 1
		}
		if (`=`table'[`i', 5]'==1) {
			loc note "*"
			loc explain_nested 1
		}

		// noabsorb
		if (`rows'==1 & `=`table'[`i', 1]'==1 & strpos("`var'", "__")==1) loc var "_cons"

		.`mytab'.row "`var'" `=`table'[`i', 1]' `=`table'[`i', 2]' `=`table'[`i', 3]' "`note'"
	}

	// Bottom
	.`mytab'.sep, bottom
	if (`explain_exact') di as text "? = number of redundant parameters may be higher"
	if (`explain_nested') di as text `"* = FE nested within cluster; treated as redundant for DoF computation"'
end
