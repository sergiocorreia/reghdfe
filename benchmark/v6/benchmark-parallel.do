// -------------------------------------------------------------------------------------------------
// based on "hard.do"
// -------------------------------------------------------------------------------------------------
* See http://cran.r-project.org/web/packages/lfe/vignettes/speed.pdf
* For the lfe (R) equivalent; example taken from there

	clear all
	set trace off
	set more off
	set matsize 11000 // for areg
	timer clear
	set type double

	local N 50000 // 1000000
	local rep 3 // 3
	local G 50000
	local H 300 // 2000 500 300
	local M 4 // 20 10 4


/* Results 01Mar2021

10:     86.47 /        2 =      43.2350
11:    149.54 /        2 =      74.7685
13:    176.93 /        2 =      88.4625
14:    183.96 /        2 =      91.9815
20:     83.94 /        2 =      41.9690
21:    115.42 /        2 =      57.7075
22:    119.73 /        2 =      59.8650
23:    113.76 /        2 =      56.8800


delta1 10: 0.0034568126
delta1 11: 0.0027094882
delta1 13: 0.0046480708
delta1 14: 0.0532000555
delta1 20: 0.0028292702
delta1 21: 0.0027264191
delta1 22: 0.0046531445
delta1 23: 0.0535316569

Basically: LSQR still not that accurate with 1e-12 tolerance, but not *that* bad

However, it seems to be always dominated by LSMR (at each given accuracy)

Some speed notes:

- MAP+CG+SYM still fastest, by a lot
- Paralle imporvements are small for MAP+CG+SYM, but quite a bit larger for all others

All in all, MAP performs much better (2x LSMR!), and even after parallelizing performs better
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
// Hard-to-Compute Dataset
// -------------------------------------------------------------------------------------------------
	clear
	set obs `N'
	set seed 1234
	forval i = 1/5 {
		gen x`i' = uniform()
	}

	smpl id1, size(`G')
	smpl id2, size(`M')
	replace id2 = 1 + mod((id1 + id2) , `H')

	gen y = 10 + x1 + x2 + x3 + x4 + sin(id1) + log(10+id2) + 50 * rnormal()
	assert !missing(y)

	su id1 id2
	sort id1 id2

global tol "1e-8"

// -------------------------------------------------------------------------------------------------

	loc cmd "reghdfe y x1 x2 x3 x4 x5, a(id1 id2) nosample dof(none)"
	loc opt_parallel `"parallel(3, cores(2) tmp("C:\Git\asd\borrar"))"'

	* True values:
	`cmd' tol(1e-14) tech(map) accel(cg) transf(sym)
	global b1 = _b[x1] 
	global b2 = _b[x2] 

	foreach timer in 10 11 12 13 14 15 16 17 18 20 21 22 23 24 35 36 39 {
		global delta1_`timer' = 0
	}

	di
	forv i=1/`rep' {
		timeit 10 : `cmd' tol($tol)  tech(map) accel(cg) transf(sym)
		*timeit 11 : `cmd' tol($tol)  tech(map) accel(cg) transf(cim)
		*timeit 13 : `cmd' tol($tol)  tech(lsmr)
		//timeit 14 : `cmd' tol(1e-12) tech(lsqr) // need to use way more tolerance
		timeit 15 : `cmd' tol($tol)  tech(map) accel(cg) transf(unrolled)

		//timeit 20 : `cmd'  `opt_parallel' tol($tol)  tech(map)  accel(cg) transf(sym)
		//timeit 21 : `cmd'  `opt_parallel' tol($tol)  tech(map)  accel(cg) transf(cim)
		//timeit 22 : `cmd'  `opt_parallel' tol($tol)  tech(lsmr)
		//timeit 23 : `cmd'  `opt_parallel' tol(1e-12) tech(lsqr)
		di
	}

	timer list

	foreach timer in 10 11 12 13 14 15 16 17 18 20 21 22 23 24 35 36 39 {
		loc delta "?"
		cap loc delta = string(${delta1_`timer'}*1e8, "%20.10f")
		if ("`delta'" == "?") continue
		di as text "delta1 `timer': `delta'"	
	}

exit
