program define reghdfe_store_alphas, eclass
	mata: st_local("save_any_fe", strofreal(HDFE.save_any_fe))
	assert inlist(`save_any_fe', 0, 1)
	if (`save_any_fe') {
		_assert e(depvar) != "", msg("e(depvar) is empty")
		_assert e(resid) != "", msg("e(resid) is empty")
		// we can't use -confirm var- because it might have TS operators
		fvrevar `e(depvar)', list
		confirm numeric var `e(resid)', exact
		tempvar d
		if (e(rank)) {
			qui _predict double `d' if e(sample), xb
		}
		else if (e(report_constant)) {
			gen double `d' = _b[_cons] if e(sample)
		}
		else {
			gen double `d' = 0 if e(sample)
		}
		qui replace `d' = `e(depvar)' - `d' - `e(resid)' if e(sample)

		mata: HDFE.store_alphas("`d'")
		drop `d'

		// Drop resid if we don't want to save it; and update e(resid)
		cap drop __temp_reghdfe_resid__
		if (!c(rc)) ereturn local resid
	}
end
