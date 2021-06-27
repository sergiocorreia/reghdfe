clear all
cls
set rmsg off

sysuse auto
reghdfe price weight, noa
clear
cls


set obs 2000000
gen double y = runiform()
forval i = 1/10 {
	gen double x`i' = runiform()
}
gen long id1 = ceil(runiform()*10000)
gen long id2 = ceil(runiform()*1000)

timer clear
forval i = 1/3 {
	timer on 1`i'
	di as text "."
	qui reghdfe y x*, a(id1 id2) dof(none) tol(1e-8) tech(lsmr) nosample
	timer off 1`i'
}

forval i = 1/3 {
	timer on 2`i'
	di as text "."
	qui reghdfe y x*, a(id1 id2) dof(none) tol(1e-8) tech(lsmr) nosample par(4, cores(2) tmp("C:\Git\asd\borrar"))
	timer off 2`i'
}

timer list
timer clear

exit

