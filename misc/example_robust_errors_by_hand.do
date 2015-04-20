clear all
set more off
cls
* ssc install moremata

set obs 100
gen x = rnormal()
gen y = 10 + 2 * x + 5 * rnormal()
gen id = (runiform()>0.8)
replace id = sum(id)
replace id = id+1
tab id

* Unadjusted S.E.
qui regress y x
matrix list e(V), format("%12.11g")

* (by hand)
predict double resid, resid
mata: N = 100
mata: K = 2
mata: X = st_data(., "x")
mata: e = st_data(., "resid")
mata: invXX = invsym(quadcross(X, 1 , X, 1))
mata: e'*e / (N-K) * invXX

* Robust S.E.
qui regress y x, vce(robust)
matrix list e(V), format("%12.11g")

* (by hand)
mata: S = X :* e , e
mata: N / (N-K) * invXX * quadcross(S,S) * invXX

* Clustered S.E.
qui regress y x, vce(cluster id)
matrix list e(V), format("%12.11g")

* (by hand)
gen double s1 = x * resid
gen double s2 = resid
collapse (sum) s1 s2, by(id)
mata: S = st_data(., tokens("s1 s2"))
mata: N_clust = rows(S)
mata: N_clust
mata: N_clust / (N_clust-1) * (N-1) / (N-K) * invXX * quadcross(S,S) * invXX

*
*li
