clear all
include reghdfe.mata, adopath

sysuse auto, clear
mata: colsum(st_data(., tokens("2.rep78 3.rep78 4.rep78 5.rep78")))
mata: colsum(st_data(., tokens("2bn.rep78 3bn.rep78 4bn.rep78 5bn.rep78")))

asd








drop if mi(rep78)
keep rep78 price weight
gen www = weight
loc vars "price   2.rep78   3.rep78   4.rep78   5.rep78    weight www"
mata: data = st_data_pool(1::st_nobs(), tokens("`vars'"), 10)
mata: cols(data)
mata: value = sum(data[., 2])
mata: check = sum(st_data(., "2.rep78"))
mata: value, ., check
mata: assert(value == check)
