* Useful to track memory usage (manually through task manager)
/* Results with 5mm obs and 10 regressors:

(in thousand K bytes)
- baseline dataset		 645
- compact pool(2)		1320		74s
- default				3100		50s

So an aggressive memory option (compact+pool) uses at most 2x memory
while the default option uses up to 4.8x memory, running at 2/3ds of the time

*/

clear all
set niceness 8
set more off
cls
sysuse auto
*reghdfe price weight, noa

clear
cls

set obs `=2e5'
gen double y = runiform()
forval i = 1/3 {
	gen double x`i' = runiform()
}
gen double v1 = runiform()
gen double v2 = runiform()
gen double v3 = runiform()

gen id1 = ceil(runiform()*1000)
gen id2 = ceil(runiform()*1000)
sort id1 id2

set rmsg on

qui groupreg y x* in 1/10, 
groupreg y x*, a(id1##c.v1 id2##c.(v2 v3)) precond(none)		// 10.8s
groupreg y x*, a(id1##c.v1 id2##c.(v2 v3)) precond(diag)		//  4.2s
groupreg y x*, a(id1##c.v1 id2##c.(v2 v3)) precond(bloc)		//  2.1s
reghdfe y x*, a(id1##c.v1 id2##c.(v2 v3))						//  6.8s



exit


reghdfe  y x* in 1/10, noa
reghdfe  y x*, a(id1 id2)

exit
