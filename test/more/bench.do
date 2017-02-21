cls
clear all
discard
set more off
set trace off


/*cap ado uninstall moresyntax
net install moresyntax, from("C:/git/moresyntax/src")

cap ado uninstall ftools
net install ftools, from("C:/git/ftools/src")
ftools, check
*/

cap ado uninstall reghdfe
net install reghdfe , from(C:/git/reghdfe/src)
reghdfe, check
mata: mata desc using lreghdfe
reghdfe, compile

sysuse auto
timer clear


// -------------------------------------------------------------------------------------------------
// Comparison of Easy and Hard Datasets
// -------------------------------------------------------------------------------------------------
* See http://cran.r-project.org/web/packages/lfe/vignettes/speed.pdf
* For the lfe (R) equivalent; example taken from there

	clear all
	set trace off
	set more off
	timer clear

	local N 1000000 // Don't exactly use 1MM due to Stata Bug
	local rep 0
	local G 100000
	local H 10000 // 500 or 300
	local M 5 // 4 or 10

cap pr drop smpl
pr smpl
	syntax newvarname, size(integer)
	gen long `varlist' = 1 + int(uniform() * `size')
	compress `varlist'
end

* ---------------------


	set obs `N'
	set seed 1234
	gen x1 = uniform()
	gen x2 = uniform()

	smpl id1, size(`G')
	smpl id2, size(`H')

	gen y = 10 + x1 + x2 + sin(id1) + log(10+id2) + 100 * rnormal()
	assert !missing(y)

	su id1 id2
	sort id1 id2

	gen byte weights = 1












timer clear


timer on 1
reghdfe y x1 x2 /*[fw=weights]*/, a(id1 id2) fast verbose(2) dof(none) // accel(none) transform(kac) // maxiter(1) tol(.9)
timer off 1

timer on 2
reghdfe y x1 x2 , a(id1 id2) verbose(2)
timer off 2

timer on 3
reghdfe y x1 x2 /*[fw=weights]*/, a(id1 id2) fast verbose(2) dof(none) // accel(none) transform(kac) // maxiter(1) tol(.9)
timer off 3

timer on 4
reghdfe y x1 x2 , a(id1 id2) verbose(2)
timer off 4




timer list
timer clear
exit

timer on 2
	// fixed_effects(ABSVARS, TOUSE, WEIGHTTYPE, WEIGHTVAR, DROPSING, VERBOSE)
	timer on 11
	//mata: HDFE = fixed_effects("id1 id2", "", "fweight", "weights", 1, 1)
	mata: HDFE = fixed_effects("id1 id2", "", "", "", 1, 1)
	timer off 11

	timer on 12
	//mata: HDFE.acceleration = "steepest_descent"
	// why is none faster?!?!
	mata: HDFE.acceleration = "none"
	mata: HDFE.acceleration = "test"
	mata: HDFE.transform = "kaczmarz"
	mata: y = HDFE.partial_out("y x1 x2")
	timer off 12

	timer on 13
	mata: mata: quicksolve(HDFE, y)
	timer off 13

timer off 2

timer list	

timer clear
exit
