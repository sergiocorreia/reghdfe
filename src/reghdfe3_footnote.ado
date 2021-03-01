// -------------------------------------------------------------
// Display Regression Footnote
// -------------------------------------------------------------

program reghdfe3_footnote
syntax [, linesize(int 79)]

	local skip1 = max(`s(width_col1)'-1, 12) // works with both _coef_table, ivreg2 and ivregress

if ("`e(model)'"=="ols" & inlist("`e(vce)'", "unadjusted", "ols")) {
	local dfa1  = e(df_a) + 1
	local todisp `"F(`=e(df_a)-1', `e(df_r)') = "'
	local skip3 = max(23-length(`"`todisp'"')-2,0)
	local skip2 = max(14-length(`"`dfa1'"')-2,0)
	local skip0 `skip1'

	foreach fe in `e(extended_absvars)' {
		local skip1 = max(`skip1', length("`fe'"))
	}

	di as text %`skip0's "Absorbed" " {c |}" ///
		_skip(`skip3') `"`todisp'"' ///
		as res %10.3f e(F_absorb) %8.3f fprob(e(df_a),e(df_r),e(F_absorb)) ///
		as text _skip(13) `"(Joint test)"'

	* Col width
	local WX = `skip1' + 1

	* Show by-fe FStats
	* Relevant macros: NUM_FE, FE1, .., FE_TARGET1, .., FE_VARLIST
	local r2 = 1 - e(rss0)/e(tss)
	local r2_report %4.3f `r2'
	forval i = 1/`e(N_hdfe_extended)' {
		local fe : word `i' of `e(extended_absvars)'
		if (e(F_absorb`i')<.) {
			di as text %`skip1's "`fe'" " {c |}" _continue
			
			local df_a_i = e(df_a`i') - (`i'==1)
			local df_r_i = e(df_r`i')
			local todisp `"F(`df_a_i', `df_r_i') = "'
			local skip3 = max(23-length(`"`todisp'"')-2,0)
			di as text _skip(`skip3') `"`todisp'"' _continue
			
			di as res %10.3f e(F_absorb`i') %8.3f fprob(e(df_a`i'),e(df_r`i'),e(F_absorb`i')) _continue
			di as text _skip(12) `"(Nested test)"'

			local r2 = 1 - e(rss`i')/e(tss)
			local r2_report `r2_report' " -> " %4.3f `r2'
			*local cats = e(K`i') - e(M`i')
			*local data = "`e(K`i')' categories, `e(M`i')' collinear, `cats' unique"
			*local skip = 62 - length("`data'")
			*di as text _skip(`skip') `"(`data')"'
		}
	}
	di as text "{hline `=1+`skip0''}{c BT}{hline 64}"
	if (e(rss0)<.) di as text " R-squared as we add HDFEs: " `r2_report'
} // regress-unadjusted specific
else {
	foreach fe in `e(absvars)' {
		local skip1 = max(`skip1', length("`fe'"))
	}
	local WX = `skip1' + 1
}

* Show category data
di as text
di as text "Absorbed degrees of freedom:"
di as text "{hline `WX'}{c TT}{hline 49}{c TRC}"   // {c TT}{hline 14}"
di as text %`skip1's "Absorbed FE" " {c |}" ///
	%13s "Num. Coefs." ///
	%16s "=   Categories" ///
	%15s "-   Redundant" ///
	"     {c |} " _continue

// if ("`e(corr1)'"!="") di as text %13s "Corr. w/xb" _continue
di as text _n "{hline `WX'}{c +}{hline 49}{c RT}"  // {c +}{hline 14}"

	local i 0
	local explain_exact 0
	local explain_nested 0

	forval i = 1/`e(N_hdfe_extended)' {
		local fe : word `i' of `e(extended_absvars)'


		di as text %`skip1's "`fe'" " {c |}" _continue
		local numcoefs = e(K`i') - e(M`i')
		assert `numcoefs'<. & `numcoefs'>=0
		local note = cond(`e(M`i'_exact)'==0, "?", " ")
		if ("`note'"=="?") {
			local explain_exact 1
		}
		else if (`e(M`i'_nested)'==1) {
			local note *
			local explain_nested 1
		}
		
		di as text %13s "`numcoefs'" _continue
		di as text %16s "`e(K`i')'" _continue
		
		di as text %15s "`e(M`i')'" _continue
		di as text %2s "`note'" "   {c |} " _continue
		//if ("`e(corr`i')'"!="") {
		//	di as text %13.4f `e(corr`i')' _continue
		//}
		di
	}
di as text "{hline `WX'}{c BT}{hline 49}{c BRC}" // {c BT}{hline 14}"
if (`explain_exact') di as text "? = number of redundant parameters may be higher"
if (`explain_nested') di as text `"* = fixed effect nested within cluster; treated as redundant for DoF computation"'
// di as text _skip(4) "Fixed effect indicators: " in ye "`e(absvars)'"

end
