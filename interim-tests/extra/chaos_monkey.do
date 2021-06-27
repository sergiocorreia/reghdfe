exit
* https://www.wikiwand.com/en/Chaos_Monkey
* COPIED FROM PPMLHDFE SOURCE; USE WITH CARE

ppmlhdfe, reload
pr drop _all

*cd "C:\Git\ppml_hdfe_demo"
cls
set more off
*use "C:\Git\ppml_hdfe_demo\slopes3", clear
*gen double index = _n
use "C:\Git\ppml_hdfe_demo\slopes3c", clear
cou

while (1) {
	loc N = c(N)
	di as text "."
	preserve
	loc seed = ceil(runiform()*100000)
	set seed `seed'
	qui drop if runiform() > .9

	if (`N'==c(N)) {
		di as text "ABORTING 1"
		restore
		continue
	}

	cou
	*cap noi ppmlhdfe y x, absorb(id1 id2 id3##c.z, maxiter(400)) maxiter(50) // v(-1) 
	cap noi ppmlhdfe y x, absorb(id1 id3##c.z, maxiter(500)) maxiter(100) // v(-1) 
	loc rc = c(rc)
	assert `rc'!=1
	if (`rc'!=1234) {
		di as text "ABORTING 2"
		restore
		continue
	}
	//if (inlist(c(rc), 0, 2001)) {
	//	di as text "ABORTING"
	//	restore
	//	continue
	//}
	//else {
	//	di c(rc)
	//	assert(inlist(c(rc), 3498, 430))
	//	di as text "FOUND"
	//}
	restore, not

	foreach id of varlist id* {
		*di as text "`id'"
		qui bys `id': drop if _N==1
		gegen max_y = max(y), by(`id')
		qui drop if max_y==0
		drop max_y
	}
	foreach id of varlist id* {
		*di as text "`id'"
		qui bys `id': drop if _N==1
		gegen max_y = max(y), by(`id')
		qui drop if max_y==0
		drop max_y
	}
	foreach id of varlist id* {
		*di as text "`id'"
		qui bys `id': drop if _N==1
		gegen max_y = max(y), by(`id')
		qui drop if max_y==0
		drop max_y
	}	
	sort index
	cou

	*assert(0)
	qui save "C:\Git\ppml_hdfe_demo\slopes3c", replace
}


exit

// --------------------------------------------------------------------------
// Drop obs one by one 
// --------------------------------------------------------------------------

clear all
cls
use "C:\Git\ppml_hdfe_demo\slopes3d", clear
loc n = c(N)
forv i=1/`n' {
	qui use "C:\Git\ppml_hdfe_demo\slopes3d", clear
	qui drop in `i'
	cap ppmlhdfe y x, absorb(id1 id2##c.z)
	di as text "[`i'] `c(rc)'"
	if (c(rc) !=1234) continue
	assert(0)
}
cou	

// --------------------------------------------------------------------------
// Drop IDs one by one
// --------------------------------------------------------------------------

clear all
cls
use "C:\Git\ppml_hdfe_demo\slopes3d", clear
levelsof id2, loc(levels)

foreach level of local levels {
	qui use "C:\Git\ppml_hdfe_demo\slopes3d", clear
	qui drop if id2==`level'
	cap ppmlhdfe y x, absorb(id1 id2##c.z)
	di as text "[`i'] `c(rc)'"
	if (c(rc) !=1234) continue
	assert(0)
}
