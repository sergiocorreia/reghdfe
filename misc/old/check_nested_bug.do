clear all
set more off
cls
* ssc install moremata

set obs 100
gen x = rnormal()
gen y = 10 + 2 * x + 5 * rnormal()
gen id = (runiform()>0.05)
replace id = sum(id)
replace id = id+1
bys id: gen byte singleton = _N==1
tab singleton
bys id: gen t = _n
xtset id t

* Unadjusted S.E.
qui xtreg y x, fe
matrix V = e(V)

* (by hand)
predict double resid, e
mata: KK = st_numscalar("e(df_a)") + 1
qui areg x, absorb(id)
predict double x_resid, resid

mata: N = 100
mata: K = 1
mata: X = st_data(., "x_resid")
mata: e = st_data(., "resid")
mata: invXX = invsym(quadcross(X, 0 , X, 0))

mata: e'*e / (N-K-KK) * invXX
matrix list V, format("%12.11g")

* Robust S.E.
* (by hand)
mata: S = X :* e
mata: invXX * quadcross(S,S) * invXX
mata: N, K, KK
mata: N / (N-K-KK) * invXX * quadcross(S,S) * invXX

* Clustered S.E.
qui xtreg y x, fe vce(cluster id)
matrix V = e(V)

* (by hand)
gen double s1 = x_resid * resid
preserve
collapse (sum) s1, by(id)
mata: S = st_data(., tokens("s1"))
mata: N_clust = rows(S)
mata: N_clust
mata: N_clust / (N_clust-1) * (N-1) / (N-K-1) * invXX * quadcross(S,S) * invXX
// WHY N-K-1??? What's the point of the -1?
// IF we are already in within transformation, why not do (N)/(N-K)
matrix list V, format("%12.11g")
restore

**** DROP SINGLETONS ****
drop if singleton
drop singleton resid s1 x_resid



* Unadjusted S.E.
qui xtreg y x, fe
matrix V = e(V)



* (by hand)
mata: N = st_numscalar("c(N)")
predict double resid, e
mata: e = st_data(., "resid")
mata: KK = st_numscalar("e(df_a)") + 1

qui areg x, absorb(id)
predict double x_resid, resid
mata: X = st_data(., "x_resid")

mata: K = 1
mata: invXX = invsym(quadcross(X, 0 , X, 0))

mata: e'*e / (N-K-KK) * invXX
matrix list V, format("%12.11g")


* Robust S.E.

* (by hand)
mata: S = X :* e
mata: invXX * quadcross(S,S) * invXX
mata: N / (N-K-KK) * invXX * quadcross(S,S) * invXX

* Clustered S.E.
qui xtreg y x, fe vce(cluster id)
matrix V1 = e(V)
qui areg y x, absorb(id) vce(cluster id)
matrix V2 = e(V)

* (by hand)
gen double s1 = x_resid * resid
preserve
collapse (sum) s1, by(id)
mata: S = st_data(., tokens("s1"))
mata: N_clust = rows(S)
mata: N_clust
mata: N_clust / (N_clust-1) * (N-1) / (N-K-1) * invXX * quadcross(S,S) * invXX
// WHY N-K-1??? What's the point of the -1?
// IF we are already in within transformation, why not do (N)/(N-K)
matrix list V1, format("%12.11g")
matrix list V2, format("%12.11g")
restore



