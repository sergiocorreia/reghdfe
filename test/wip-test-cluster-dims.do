discard
pr drop _all
clear
sysuse auto
drop if rep==.
cls

bys turn: gen t = _n
tsset turn t

* Assert that mwc and avar give the same results for 1
* and for 1 add default suite
reghdfe price weight, a(foreign) vce(cluster turn t, bw(2) kernel(par) suite(avar))

exit
ereturn list

di e(vce)
di e(vcetype)
di e(vcesuite)
di e(title2)
di e(dofadjustments)
di e(clustvar)

di e(N_clust)
di e(N_clust1)
di e(N_clust2)
di e(N_clustervars)
di e(N_hdfe)
