* Cleanup
	clear all
	cls
	set trace off
	set more off

* Parameters
	local N = 2e6 // 6
	local K = 20
	local G = 3
	local groupsizes 5 10 15 20
	
* Build dataset	
	set obs `N'

	forv i=1/`G' {
		gen u = uniform()
		sort u
		gen long id`i' = int((_n-1)/1000)
		local absvars `absvars' id`i'
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

foreach pool of local groupsizes {
	preserve
	local cmd reghdfe `vars' , a(`absvars') timeit groupsize(`pool') savecache dof(none) keepsingletons
	di as error `"[CMD] `cmd'"'
	timer on 1
	`cmd'
	timer off 1
	timer list 1
	timer clear 1
	di as error "POOL was `pool'" _n
	restore
}
exit




/* RESULTS

Parameters used:
	N = 2MM obs
	K = 10
	#FE = 3
	
Time by groupsize:
	gs=1 : 81.1s
	gs=2 : 58.4s
	gs=5 : 41.8s
	gs=10: 39.2s
	
What's the pattern?
Instead of tinking of the groupsize, think of the number of repetitions (K/gs)
Then, we have a linear trend, and increasing number of repetitions by 2 increases
the time by 40% (up to the point where the other costs dominate)

Thus, there may not be much difference at all between a gs of 10 and a gs of 5?

To see this last question, I increased K to 20 and tested again:

gs	reps	time	% fastest
5	4.00	82.00	131%
10	2.00	74.68	120%
15	1.33	66.50	107%
20	1.00	62.43	100%

So there is still an important speedup from a higher gs)

*/
