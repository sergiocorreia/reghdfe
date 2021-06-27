clear all
cls
adopath + "C:\Git\groupreg\test"

use "C:\Users\Sergio\Downloads\error_data_set"

rename highly_cited y
rename team_created x
rename tech_class id
rename refined_id indiv_id
rename patnum group_id

* reghdfe y x, a(id indiv_id) group(group_id) indiv(indiv_id)

save "C:\Users\Sergio\Downloads\error04mar2021", replace

chaos_drop_obs, f("C:\Users\Sergio\Downloads\error04mar2021") prob(0.1) maxiter(100) minobs(10) allow(3456): ///
	reghdfe y x, a(id indiv_id) group(group_id) indiv(indiv_id)
