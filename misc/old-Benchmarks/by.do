// QUESTION: Is by() worth it?
// ANSWER: It seems NO; slower than doing the regression by parts
// (Note however than savecache+usecache is still useful)

* Cleanup
	clear all
	cls
	set trace off
	set more off

* Parameters
	local N = 2e6
	local K = 3
	local G = 3
	local H 5
	local pool 5
	
* Build dataset	
	set obs `N'
	
	gen byte group = 1 + int((_n-1)/c(N)*`H')
	tab group

	local by_absvars group
	forv i=1/`G' {
		gen u = uniform()
		sort u
		gen long id`i' = int((_n-1)/1000)
		local absvars `absvars' id`i'
		local by_absvars `by_absvars' group#id`i'
		drop u
	}

	local vars y
	forv i=1/`K' {
		local vars `vars' x`i'
	}

	foreach var of local vars {
		gen double `var' = uniform()
	}

* Test only the demeaning part (not the regression which is unrelated

* Options:
* 1) reghdfe+save+by, then reghdfe_use+if or preserve+keep+regress+restore
* 2) Use -if- , reghdfe_save, then reghdfe_use, then use -if- again ...

tempfile data
save "`data'"

* 1)
*use group using "`data'" 
*levelsof group
timer on 1
forv i=1/`H' {
	use if group==`i' using "`data'", clear
	drop group
	local cmd reghdfe `vars' , a(`absvars') groupsize(`pool') savecache dof(none) keepsingletons
	di as error `i' // `"[CMD] `cmd'"'
	qui `cmd'
	if (`i'==2) reghdfe `vars', a(`absvars') usecache
}
timer off 1

* 2)
timer on 2
	use "`data'", clear
	local cmd reghdfe `vars' , a(`by_absvars') groupsize(`pool') savecache dof(none) keepsingletons
	di as error `"[CMD] `cmd'"'
	qui `cmd'
	regress `vars' if __ID1__==2 , nocons
timer off 2

timer list
timer clear
exit

asd
foreach pool of local groupsizes {
	preserve
	local cmd reghdfe `vars' , a(`absvars') timeit groupsize(`pool') savecache dof(none) keepsingletons
	di as error `"[CMD] `cmd'"'
	timer on 1
	`cmd'

	di as error "POOL was `pool'" _n
	restore
}
exit

/* RESULTS

Parameters used:
	local N = 2e6
	local K = 3
	local G = 3
	local H 5
	local pool 5
	
Independent regressions took 21s, by() regressions took 32s

*/
