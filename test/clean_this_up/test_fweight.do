/*
clear all
discard
cls
set more off

sysuse auto
keep if rep<.

replace turn = int(turn/14)
collapse (count) n=price (mean) price, by(turn rep foreign)

reghdfe price foreign [aw=n] , a(i.turn##c.rep)
local b = _b[foreign]
reg price foreign i.turn i.turn#c.rep [aw=n]
assert abs(`b'-_b[foreign])<0.01


exit
*/

clear all
discard
cls
set more off
sysuse auto

gen n = 100000 // 0000 // int(uniform()*10+3)
reghdfe price weight [fw=n], a(turn##c.length)
exit
areg price weight turn#c.length [fw=n], a(turn)





* THIS IS THE SIMPLEST ILLUSTRATION OF THE BUG

clear all
discard
cls
set more off

sysuse auto
keep if rep<.

collapse (count) n=price (mean) price, by(foreign turn)

gen c = 1
reghdfe price foreign [fw=n] , a(turn turn#c.turn)
reghdfe price foreign [fw=n] , a(turn##c.turn)
reg price foreign i.turn##c.turn [fw=n]



discard
clear all
set obs 10
gen c = 1
set seed 32453245
gen y = uniform()
gen x = uniform()
gen v = uniform()
gen n = 1 + int(10*uniform())
reghdfe y x [fw=n] , a(i.c i.c#c.v)
reg y x v [fw=n]





