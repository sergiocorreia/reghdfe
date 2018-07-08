* Useful to track memory usage

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
sleep 2000
reghdfe y x*, a(id1 id2) compact
sleep 2000
reghdfe y x*, a(id1 id2) compact pool(5)
sleep 2000
reghdfe y x*, a(id1 id2)

exit
