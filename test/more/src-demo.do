cls
clear all
discard
set more off
set trace on

cap ado uninstall ftools
net install ftools, from("C:/git/ftools/src")
ftools, check
ftools, compile

cap ado uninstall reghdfe
net install reghdfe , from(c:\git\reghdfe/src)

reghdfe, check
mata: mata desc using lreghdfe


sysuse auto, clear

set trace off
reghdfe, compile
/*
mata: HDFE = fixed_effects("turn", "", 0, 1)
mata: HDFE = fixed_effects("turn", "", 1, 1)
mata: HDFE = fixed_effects("turn trunk foreign", "", 1, 1)
mata: HDFE = fixed_effects("turn trunk", "", 1, 1)
mata: HDFE = fixed_effects("turn trunk", "foreign", 1, 1)
mata: HDFE = fixed_effects("turn trunk foreign", "", 0, 1)
*/
reghdfe, check



* ---------------------
findfile "ftools_type_aliases.mata"
include "`r(fn)'"
local FixedEffects      class FixedEffects scalar
mata:
mata set matastrict off
void quicksolve(`FixedEffects' HDFE, `Variables' y)
{
    w = HDFE.has_weights ? HDFE.weight : 1
    X = y[., 2..cols(y)]
    y = y[., 1]
    //invsym(X' * X) * X' * y

    // see http://fmwww.bc.edu/repec/bocode/r/ranktest.ado wf
    xx = quadcross(X, w, X)
    xy = quadcross(X, w, y)
    reghdfe_cholqrsolve(xx, xy)
    rows(X), cols(X)
}
end
* ---------------------




loc vars "price weight head"
loc absvars "foreign turn#trunk"
loc absvars "foreign turn#trunk"
loc absvars "foreign##c.gear"
loc absvars "turn foreign##c.(gear length)"
loc touse ""

// fixed_effects(ABSVARS, TOUSE, WEIGHTTYPE, WEIGHTVAR, DROPSING, VERBOSE)

reghdfe `vars', a(`absvars') fast keepsing
mata: HDFE = fixed_effects("`absvars'", "", "", "", 0, 1)
mata: y = HDFE.partial_out("`vars'")
mata: quicksolve(HDFE, y)

reghdfe `vars' [fw=displacement], a(`absvars') fast  tol(1e-10) keepsing
mata: HDFE = fixed_effects("`absvars'", "", "fweight", "disp", 0, 1)
//mata: HDFE.acceleration = "steepest_descent" // steepest_descent
//mata: HDFE.maxiter=10
//mata: HDFE.verbose=4
//mata: HDFE.transform = "cimmino"
mata: y = HDFE.partial_out("`vars'")
mata: quicksolve(HDFE, y)

exit
