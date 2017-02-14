cap pr drop SaveFE
pr SaveFE
	syntax, model(string) depvar(string) untransformed(string) subpredict(string) ///
		has_intercept(integer) ///
		[weightexp(string)] [drop_resid_vector(integer 1)]

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
	mata: resid2dta(HDFE_S, 0, `drop_resid_vector')

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

	* Computing mean(resid_d), the constant term (only if there is an intercept in absorb)
	if (`has_intercept') {
		local tmpweightexp = subinstr("`weightexp'", "[pweight=", "[aweight=", 1)
		su `resid_d' `tmpweightexp', mean
		gen double `d' = `resid_d' - r(mean) - `resid'
	}
	else {
		gen double `d' = `resid_d' - `resid'
	}	
	
	drop `resid' `resid_d'
	//clonevar dd = `d'

	Debug, level(3) msg(" - disaggregating d = z1 + z2 + ...")
	mata: map_save_fe(HDFE_S, "`d'")
	//regress dd __hdfe*, nocons
	drop `d'
end
