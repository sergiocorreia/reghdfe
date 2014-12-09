clear all
discard
set more off
sysuse auto
cls
set trace off
set tracedepth 3
rename head headroo
reghdfe price head, a(B=tu A=rep##c.length) verbose(3) avge(C=foreign) vce(cluster tu#foreign) // 



predict resid, resid
exit


reghdfe, alt
exit

