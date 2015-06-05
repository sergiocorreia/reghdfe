cls
set more off
clear all
*set processors 1

cap log close _all
log using res2fe.log, replace
version

cap ado uninstall reghdfe
net from https://raw.githubusercontent.com/sergiocorreia/reghdfe/updated_mata/package/
net install reghdfe

set obs 2000000
gen y = uniform()
gen x1 = uniform()
gen x2 = uniform()
gen x3 = uniform()
gen id = 1 + int((_n-1)/1000)
bys id: gen long t = _n
compress
xtset id t

local vce vce(cluster id)

set rmsg on

areg y x* , absorb(id) `vce'
xtreg y x*, fe `vce'
reghdfe y x* , absorb(id) `vce' old // v2
reghdfe y x* , absorb(id) `vce' // v3-slow
reghdfe y x* , absorb(id) `vce' fast // dof(none) keepsingletons // v3-fast

log close _all

set rmsg off
exit
