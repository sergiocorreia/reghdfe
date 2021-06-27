* ===========================================================================
* 
* ===========================================================================

clear
input byte(y x year group_id indiv_id)
1 1 1 4 2
0 1 1 6 2
0 1 . 7 3
0 1 1 8 3
end

reghdfe y x if !mi(year), a(year indiv_id) group(group_id) indiv(indiv_id) // works
gen sample = e(sample)
li

reghdfe y x, a(year indiv_id) group(group_id) indiv(indiv_id) // does not work!
