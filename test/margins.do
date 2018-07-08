* Not a proper cscript yet...

clear all
* cls
sysuse auto


// --------------------------------------------------------------------------
// Initial examples
// --------------------------------------------------------------------------
qui areg price c.weight#c.weight c.weight#c.gear i.foreign, a(turn)
margins
margins foreign
margins, dydx(weight) atmeans

qui reghdfe price c.weight#c.weight c.weight#c.gear i.foreign, a(turn) keepsing
margins
margins foreign
margins, dydx(weight) atmeans


// --------------------------------------------------------------------------
// Examples from Richard DiSalvo
// --------------------------------------------------------------------------
qui areg price c.weight##c.headroom, a(turn)
margins, dydx(weight) at(headroom=(1.5(0.5)5))

qui reghdfe price c.weight##c.headroom, a(turn) keepsing
margins, dydx(weight) at(headroom=(1.5(0.5)5))

qui areg price weight headroom c.weight#c.headroom, a(turn)
margins, dydx(weight) at(headroom=(1.5(0.5)5))

qui reghdfe price weight headroom c.weight#c.headroom, a(turn) keepsing
margins, dydx(weight) at(headroom=(1.5(0.5)5))

qui areg price c.weight##i.foreign, a(turn)
margins, dydx(weight) at(foreign = 1)

qui reghdfe price c.weight##i.foreign, a(turn) keepsing
margins, dydx(weight) at(foreign = 1)


// --------------------------------------------------------------------------
// Examples from Jeff Pitblado (StataCorp)
// --------------------------------------------------------------------------
areg price c.weight##c.headroom, a(turn)
est store areg1
margins, dydx(weight) at(headroom=(1.5(0.5)5)) post
est store areg1_marg

reghdfe price c.weight##c.headroom, a(turn)
est store reghdfe1
margins, dydx(weight) at(headroom=(1.5(0.5)5)) post
est store reghdfe1_marg

* show equivalent model fits, and the same parameterization
est table areg1 reghdfe1, b se stat(rmse)

* show -margins- yields the same marginal effects and SE estimates
est table areg1_marg reghdfe1_marg, b se

areg price c.weight##i.foreign, a(turn)
est store areg2
margins for, dydx(weight) post
est store areg2_marg

reghdfe price c.weight##i.foreign, a(turn)
est store reghdfe2
margins for, dydx(weight) post
est store reghdfe2_marg

* show equivalent model fits, but a different parameterization because
* the c.weight#1.foreign coefficient was fit by -areg-, but
* c.weight#0.foreign was fit by -reghdfe-
est table areg2 reghdfe2, b se stat(rmse)

* show -margins- yields the same marginal effects and SE estimates
est table areg2_marg reghdfe2_marg, b se
