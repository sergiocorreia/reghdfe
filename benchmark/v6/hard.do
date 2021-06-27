// -------------------------------------------------------------------------------------------------
// Comparison of Easy and Hard Datasets
// -------------------------------------------------------------------------------------------------
* See http://cran.r-project.org/web/packages/lfe/vignettes/speed.pdf
* For the lfe (R) equivalent; example taken from there

	clear all
	set trace off
	set more off
	set matsize 11000 // for areg
	timer clear

	local N 1000000
	local rep 2
	local G 50000
	local H 2000 // 500 or 300
	local M 20 // 4 or 10

/*
--------------------------------------------------------------------------------
INITIAL RUN circa 2016:
--------------------------------------------------------------------------------

	local N 100000
	local rep 2

	* Hardest
	local G 10000
	local H 500
	local M 5
	.         timer list
	  10:      1.61 /        1 =       1.6090 easy cg
	  11:      1.03 /        1 =       1.0300 easy aitken
	  19:      4.63 /        1 =       4.6350 easy areg
	  20:     16.42 /        1 =      16.4200 hardest cg
	  21:     25.47 /        1 =      25.4670 hardest aitken
	  29:      5.05 /        1 =       5.0540 hardest areg

	* Hard G=10k H=300 M=10
	  10:      1.51 /        1 =       1.5130
	  11:      1.06 /        1 =       1.0650
	  19:      1.45 /        1 =       1.4500
	  20:      5.88 /        1 =       5.8790
	  21:     16.45 /        1 =      16.4540
	  29:      1.56 /        1 =       1.5600

	  
 -------------------------------------------------------------------------------
 UPDATE feb2017:
--------------------------------------------------------------------------------

 N=100,000 rep=3 (the times below are for three reps!!!)

	G=10,000 H=500 M=5
	  10:      1.87 /        1 =       1.8730
	  11:      1.03 /        1 =       1.0330
	  19:     10.96 /        1 =      10.9600
	  20:      9.67 /        1 =       9.6700
	  21:     25.06 /        1 =      25.0640
	  29:     11.49 /        1 =      11.4860

	G=10,000 H=300 M=10
	  10:      1.86 /        1 =       1.8590
	  11:      1.11 /        1 =       1.1100
	  19:      4.92 /        1 =       4.9220
	  20:      3.78 /        1 =       3.7780
	  21:     12.06 /        1 =      12.0630
	  29:      4.72 /        1 =       4.7190

*/


// -------------------------------------------------------------------------------------------------
// Programs
// -------------------------------------------------------------------------------------------------
capture program drop smpl
program define smpl
	syntax newvarname, size(integer)
	gen long `varlist' = 1 + int(uniform() * `size')
	compress `varlist'
end


// -------------------------------------------------------------------------------------------------
// Easy-to-Compute Dataset
// -------------------------------------------------------------------------------------------------
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



// -------------------------------------------------------------------------------------------------
// Hard-to-Compute Dataset
// -------------------------------------------------------------------------------------------------
	clear
	set obs `N'
	set seed 1234
	gen x1 = uniform()
	gen x2 = uniform()

	smpl id1, size(`G')
	smpl id2, size(`M')
	replace id2 = 1 + mod((id1 + id2) , `H')

	gen y = 10 + x1 + x2 + sin(id1) + log(10+id2) + 100 * rnormal()
	assert !missing(y)

	su id1 id2
	sort id1 id2

global tol "1e-8"

// -------------------------------------------------------------------------------------------------

**timer clear
**set rmsg on
**reghdfe y x1 x2, a(id1 id2) nosample dof(none) tol($tol) tech(map) accel(cg) transf(sym) fastreg timeit
**exit
	
// -------------------------------------------------------------------------------------------------

	* True values:
	reghdfe y x1 x2, a(id1 id2) nosample dof(none) tol(1e-14) tech(map) accel(cg) transf(sym)
	global b1 = _b[x1] 
	global b2 = _b[x2] 

	loc opt_parallel `"parallel(3, cores(2) tmp("C:\Git\asd\borrar"))"'
	di

	foreach timer in 10 11 12 13 14 15 16 17 18 20 21 22 23 24 35 36 39 {
		global delta1_`timer' = 0
	}

	forv i=1/`rep' {
		timeit 10 : reghdfe y x1 x2, a(id1 id2) nosample dof(none) tol($tol) tech(map) accel(cg) transf(sym)
		timeit 11 : reghdfe y x1 x2, a(id1 id2) nosample dof(none) tol($tol) tech(map) accel(cg) transf(cim)
		*timeit 12 : reghdfe y x1 x2, a(id1 id2) nosample dof(none) tol($tol) tech(map) accel(a)  transf(kac)
		timeit 13 : reghdfe y x1 x2, a(id1 id2) nosample dof(none) tol($tol) tech(lsmr)
		timeit 18 : reghdfe y x1 x2, a(id1 id2) nosample dof(none) tol(1e-12) tech(lsqr) // need to use way more tolerance

		timeit 20 : reghdfe y x1 x2, a(id1 id2) nosample dof(none) tol($tol) tech(map) accel(cg) transf(sym) `opt_parallel' 
		timeit 21 : reghdfe y x1 x2, a(id1 id2) nosample dof(none) tol($tol) tech(map) accel(cg) transf(cim) `opt_parallel'
		*timeit 22 : reghdfe y x1 x2, a(id1 id2) nosample dof(none) tol($tol) tech(map) accel(a)  transf(kac) `opt_parallel'
		timeit 23 : reghdfe y x1 x2, a(id1 id2) nosample dof(none) tol($tol) tech(lsmr) `opt_parallel'
		timeit 24 : reghdfe y x1 x2, a(id1 id2) nosample dof(none) tol(1e-12) tech(lsqr) `opt_parallel'

		timeit 35 : reghdfe y x1 x2, a(id1 id2) nosample dof(none) tol($tol) accel(cg) transf(sym) version(5) nowarn
		*timeit 36 : reghdfe y x1 x2, a(id1 id2)          dof(none) tol($tol) accel(cg) transf(sym) version(3) nowarn

		*timeit 39 : areg y x* i.id2, a(id1)
		di
	}

	timer list

	foreach timer in 10 11 12 13 14 15 16 17 18 20 21 22 23 24 35 36 39 {
		loc delta "?"
		cap loc delta = string(${delta1_`timer'}*1e8, "%20.10f")
		di as text "delta1 `timer': `delta'"	
	}

exit

// -------------------------------------------------------------------------------------------------
// Hard-to-Compute Dataset
// -------------------------------------------------------------------------------------------------
	clear
	set obs `N'
	set seed 1234
	gen x1 = uniform()
	gen x2 = uniform()

	smpl id1, size(`G')
	smpl id2, size(`M')
	replace id2 = 1 + mod((id1 + id2) , `H')
	tab id2

	gen y = 10 + x1 + x2 + sin(id1) + log(10+id2) + 100 * rnormal()
	assert !missing(y)

	su id1 id2
	sort id1 id2

// -------------------------------------------------------------------------------------------------

	forv i=1/`rep' {
		timer on 20
		di as text "." _c
		qui reghdfe y x1 x2, a(id1 id2) nosample dof(none) tol($tol) accel(cg) transf(sym)
		timer off 20
	}

	forv i=1/`rep' {
		timer on 21
		di as text "." _c
		qui reghdfe y x1 x2, a(id1 id2) nosample dof(none) tol($tol) accel(a) transf(kac)
		timer off 21
	}


	forv i=1/`rep' {
		timer on 29
		qui areg y x* i.id2, a(id1)
		timer off 29
	}

	
	timer list

exit

// -------------------------------------------------------------------------------------------------

	* This shows why CG is better
	set rmsg on
	reghdfe y x1 x2, a(id1 id2) fast dof(none) tol(1e-8) accel(cg) transf(sym) old
	reghdfe y x1 x2, a(id1 id2) nosample dof(none) tol(1e-8) accel(cg) transf(sym) // v(4) time
	reghdfe y x1 x2, a(id1 id2) nosample dof(none) tol(1e-8) accel(cg) transf(cim) // v(4) time
	reghdfe y x1 x2, a(id1 id2) nosample dof(none) tol(1e-8) accel(sd) transf(sym) // v(4) time
	cap noi reghdfe y x1 x2, a(id1 id2) nosample dof(none) tol(1e-8) accel(a) transf(kac) // Stuck
	* slow reghdfe y x1 x2, a(id1 id2) fast dof(none) tol(1e-8) accel(none) transf(kac) // v(4) time // 46
	set rmsg off

exit

