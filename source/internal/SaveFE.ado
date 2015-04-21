capture program drop SaveFE
program define SaveFE
	syntax, model(string) depvar(string) untransformed(string) subpredict(string) [weightexp(string)]

	Debug, level(2) msg("(calculating fixed effects)")
	tempvar resid
	local score = cond("`model'"=="ols", "score", "resid")
	Debug, level(3) msg(" - predicting resid (equation: y=xb+d+cons+resid)")
	if e(df_m)>0 {
		`subpredict' double `resid', `score' // equation: y = xb + d + e, we recovered "e"
	}
	else {
		gen double `resid' = `depvar'
	}
	mata: store_resid(HDFE_S, "`resid'")

	Debug, level(3) msg(" - reloading untransformed dataset")
	qui use "`untransformed'", clear
	erase "`untransformed'"
	mata: resid2dta(HDFE_S)

	Debug, level(3) msg(" - predicting resid+d+cons (equation: y=xb+d+cons+resid)")
	tempvar resid_d
	if e(df_m)>0 {
		`subpredict' double `resid_d', `score' // This is "d+e" (including constant)
	}
	else {
		gen double `resid_d' = `depvar'
	}

	Debug, level(3) msg(" - computing d = resid_d - mean(resid_d) - resid")
	tempvar d
	local tmpweightexp = subinstr("`weightexp'", "[pweight=", "[aweight=", 1)
	// if ("`weightvar'"!="") assert "`tmpweightexp'"!=""
	su `resid_d' `tmpweightexp', mean
	gen double `d' = `resid_d' - r(mean) - `resid'
	drop `resid' `resid_d'
	//clonevar dd = `d'
	Debug, level(3) msg(" - disaggregating d = z1 + z2 + ...")
	mata: map_solve(HDFE_S, "`d'", "", "", 1) // Store FEs in Mata (will fail if partial is set)
	//regress dd __hdfe*, nocons
	drop `d'
end
