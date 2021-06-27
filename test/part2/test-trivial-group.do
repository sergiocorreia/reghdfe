sysuse auto, clear

gen long i = _n
reghdfe price weight, a(turn) indiv(turn) group(i) precond(no)
reghdfe price weight, a(turn) version(5)

exit
