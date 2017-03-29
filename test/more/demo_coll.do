cls
clear all
discard
set more off
set trace off

cap ado uninstall ftools
net install ftools, from("C:/git/ftools/src")
ftools, check
ftools, compile

cap ado uninstall reghdfe
net install reghdfe , from(C:/git/reghdfe/src)
//set trace on
reghdfe, compile

sysuse auto
gen WEIGHT = weight
gen HEAD = head + rnormal()*1e-8
gen foo = (turn==35) + (foreign==1)
gen wh = weight-head
egen long TT = group(turn trunk)

loc vars "price weight head gear  HEAD WEIGHT wh"
loc absvars "foreign turn#trunk"
//loc absvars "foreign##c.gear"
//loc absvars "turn foreign##c.(gear length)"


areg price ibn.foreign weight head gear wh HEAD WEIGHT, a(TT)
reghdfe `vars', a(`absvars')


reghdfe gear_ratio price, a(`absvars')
//areg gear_ratio price ibn.foreign, a(TT)

loc vce robust // cluster turn
loc vce cluster turn trunk
loc vce cluster trunk
reghdfe price weight i.trunk if trunk <= 15, a(turn) keepsing tol(1e-12) vce(`vce')
matrix list e(V)
areg price weight i.trunk if trunk <= 15, a(turn) vce(`vce')
matrix list e(V)

//reghdfe `vars' [fw=displacement], a(`absvars') fast  tol(1e-10) keepsing

exit
