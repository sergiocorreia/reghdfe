pr drop _all
discard

clear all
set more off
set trace off
cls

cap ado uninstall ftools
net install ftools, from("C:/git/ftools/src")
ftools, compile

cap ado uninstall reghdfe
net install reghdfe , from("C:/git/reghdfe/src")
reghdfe, compile

sysuse auto

gen WEIGHT=1
gen omega = _n
set trace off
reghdfe 43.turn price  weight gear head if turn>10 in 1/73 [fw=W], a(for##c.disp trunk) verbose(1) transform() cluster(omega)


exit 

reghdfe 43.turn price  (weight=gear head) if turn>10 in 1/73 [fw=W], a(for##c.disp trunk) verbose(1) transform() cluster(omega)

exit

reghdfe price weight, a(turn trunk) verbose(1) transform(rand)



exit


findfile reghdfe.mata
do "`r(fn)'"
