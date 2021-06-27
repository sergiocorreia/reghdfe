* Test issue #194
* Compact option causes error if time/panel var has been missed

pr drop _all
clear all

sysuse auto, clear
bys turn: gen t = _n
xtset turn t

reghdfe price weight gear, noa
reghdfe price weight gear, noa compact

preserve
	drop t
	reghdfe price weight gear, noa
	reghdfe price weight gear, noa compact
restore

preserve
	drop turn
	reghdfe price weight gear, noa
	reghdfe price weight gear, noa compact
restore

preserve
	drop turn t
	reghdfe price weight gear, noa
	reghdfe price weight gear, noa compact
restore

preserve
	reghdfe price weight gear, noa
	reghdfe price weight L.gear, noa compact
	drop t
	cap noi reghdfe price weight L.gear, noa compact
	assert c(rc)
restore

di as text "Done! Test passed"
