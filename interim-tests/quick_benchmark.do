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

set obs `=3e5'
gen double y = runiform()
forval i = 1/10 {
	gen double x`i' = runiform()
}
gen id1 = ceil(runiform()*10000)
gen id2 = ceil(runiform()*1000)

set rmsg on


qui groupreg y x* in 1/10, 

groupreg y x*, a(id1 id2) precond(none)
groupreg y x*, a(id1 id2) precond(diag)
groupreg y x*, a(id1 id2) precond(bloc)
* reghdfe y x*, a(id1 id2) // still way slower than reghdfe (but not if we use reghdfe..pool(1))



exit


reghdfe  y x* in 1/10, noa
reghdfe  y x*, a(id1 id2)

exit
