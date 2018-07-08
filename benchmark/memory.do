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
reghdfe price weight, noa

clear
cls

set obs 5000000
gen double y = runiform()
forval i = 1/10 {
	gen double x`i' = runiform()
}
gen id1 = ceil(runiform()*10000)
gen id2 = ceil(runiform()*1000)

set rmsg on
sleep 10000
*reghdfe y x*, a(id1 id2) compact
*sleep 4000
reghdfe y x*, a(id1 id2) compact pool(2)
sleep 4000
reghdfe y x*, a(id1 id2)

exit
