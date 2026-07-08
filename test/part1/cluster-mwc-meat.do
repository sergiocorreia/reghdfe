noi cscript "reghdfe: multi-way cluster meat matches legacy algorithm" adofile reghdfe

clear all
set seed 123456
set obs 4317

gen byte touse = 1
gen int id_single = ceil(_n / 10)
bys id_single: gen byte t_single = _n
gen int id_unbal = floor(sqrt(_n * 13)) + mod(_n, 7)
gen int t_unbal = mod(_n + floor(_n / 17), 41)
gen double resid = cond(mod(_n, 2), 1, -1) * (sin(_n / 7) + cos(_n / 13) + mod(_n, 19) / 10)

forval i = 1/80 {
	gen double x`i' = cond(mod(_n, 17)==0, 1e7, 1) * (sin(_n / (`i' + 1)) + cos(_n * (`i' + 1) / 100)) + rnormal() / 100
}
unab xvars : x*

* Force-load the reghdfe Mata source before calling the private helper directly.
findfile reghdfe.mata
qui do "`r(fn)'"

mata:
real matrix reghdfe_test_meat_legacy(pointer(class Factor scalar) scalar FP,
                                     real matrix X,
                                     real colvector resid,
                                     real scalar Msize,
                                     real scalar report_constant)
{
	real scalar				i, N_clust
	real matrix				X_sorted
	real colvector			resid_sorted
	real matrix				X_tmp
	real colvector			resid_tmp
	real rowvector			Xe_tmp
	real matrix				M

	if (cols(X)==0 & !report_constant) return(J(0,0,0))

	N_clust = (*FP).num_levels
	(*FP).panelsetup()
	X_sorted = (*FP).sort(X)
	resid_sorted = (*FP).sort(resid)
	M = J(Msize, Msize, 0)

	if (cols(X)) {
		for (i=1; i<=N_clust; i++) {
			X_tmp = panelsubmatrix(X_sorted, i, (*FP).info)
			resid_tmp = panelsubmatrix(resid_sorted, i, (*FP).info)
			Xe_tmp = quadcross(1, 0, resid_tmp, X_tmp, report_constant)
			M = M + quadcross(Xe_tmp, Xe_tmp)
		}
	}
	else {
		assert(report_constant)
		for (i=1; i<=N_clust; i++) {
			resid_tmp = panelsubmatrix(resid_sorted, i, (*FP).info)
			M = M + quadsum(resid_tmp) ^ 2
		}
	}

	return(M)
}

void function reghdfe_test_meat_cmp(pointer(class Factor scalar) scalar FP,
                                    real matrix X,
                                    real colvector resid,
                                    real scalar report_constant,
                                    string scalar label)
{
	real scalar				delta, Msize
	real matrix				M_old, M_new

	Msize = cols(X) + report_constant
	M_old = reghdfe_test_meat_legacy(FP, X, resid, Msize, report_constant)
	M_new = reghdfe_vce_cluster_meat(FP, X, resid, Msize, report_constant)
	delta = mreldif(M_old, M_new)

	if (missing(delta) | delta >= 1e-10) {
		printf("{err}%s mismatch: mreldif = %21.16e\n", label, delta)
		_error(9)
	}
}

void function reghdfe_test_meat_pair(string scalar idvar,
                                     string scalar tvar,
                                     string scalar label)
{
	real matrix				X, X0
	real colvector			resid
	class Factor scalar		F_id, F_t, F_id_t

	X = st_data(., tokens(st_local("xvars")))
	X0 = J(rows(X), 0, .)
	resid = st_data(., "resid")

	F_id = factor(idvar, "touse", 0, "", 1, 1, ., 0)
	F_t = factor(tvar, "touse", 0, "", 1, 1, ., 0)
	F_id_t = join_factors(F_id, F_t, ., ., 1)

	reghdfe_test_meat_cmp(&F_id, X, resid, 0, label + " id/no constant")
	reghdfe_test_meat_cmp(&F_id, X, resid, 1, label + " id/constant")
	reghdfe_test_meat_cmp(&F_id, X0, resid, 1, label + " id/constant only")

	reghdfe_test_meat_cmp(&F_t, X, resid, 0, label + " t/no constant")
	reghdfe_test_meat_cmp(&F_t, X, resid, 1, label + " t/constant")
	reghdfe_test_meat_cmp(&F_t, X0, resid, 1, label + " t/constant only")

	reghdfe_test_meat_cmp(&F_id_t, X, resid, 0, label + " id#t/no constant")
	reghdfe_test_meat_cmp(&F_id_t, X, resid, 1, label + " id#t/constant")
	reghdfe_test_meat_cmp(&F_id_t, X0, resid, 1, label + " id#t/constant only")
}

void function reghdfe_test_cluster_meat()
{
	reghdfe_test_meat_pair("id_single", "t_single", "singleton")
	reghdfe_test_meat_pair("id_unbal", "t_unbal", "unbalanced")
}

reghdfe_test_cluster_meat()
end

exit
