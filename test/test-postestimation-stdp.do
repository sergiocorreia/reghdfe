noi cscript "reghdfe postestimation: predict stdp" adofile reghdfe

* Setup
	discard
	clear all
	set more off
	* cls

// -------------------------------------------------------------------------------------------------
capture program drop FixAreg
program define FixAreg
	tempname b V el
	*matrix list e(V)
	matrix `b' = e(b)
	matrix `V' = e(V)
	local names : colnames `b'
	replace se_bench = se_bench ^ 2
	replace se = se ^ 2

	foreach var of local names {
		*di "`var'"
		matrix `el' = `V'["`var'", "_cons"]
		scalar `el' = `el'[1,1]
		qui replace se_bench = se_bench - (1+("`var'"!="_cons")) * (`var') * (`el')
	}
	su se*
end
// -------------------------------------------------------------------------------------------------

* Toy case
	sysuse auto, clear
	
	reghdfe price , a(turn) keepsing
	predict double se, stdp
	
	areg price , a(turn)
	predict double se_bench, stdp
	
	FixAreg
	_vassert se se_bench, tol(1e-8)

* Bivariate case
	sysuse auto, clear
	
	reghdfe price weight, a(turn) keepsing
	predict double se, stdp
	
	areg price weight, a(turn)
	predict double se_bench, stdp
	
	FixAreg
	_vassert se se_bench, tol(1e-8)

* Multivariate case
	sysuse auto, clear
	
	reghdfe price weight gear trunk, a(turn) keepsing
	predict double se, stdp
	
	areg price weight gear trunk, a(turn)
	predict double se_bench, stdp
	
	FixAreg
	_vassert se se_bench, tol(1e-8)

* Now with weights

	sysuse auto, clear

	reghdfe price weight gear [fw=turn], a(turn) keepsing
	predict double se, stdp

	areg price weight gear [fw=turn], a(turn)
	predict double se_bench, stdp
	
	FixAreg
	_vassert se se_bench, tol(1e-8)

* Now with p.weights

	sysuse auto, clear

	reghdfe price weight gear [aw=length], a(turn) keepsing
	predict double se, stdp

	areg price weight gear [aw=length], a(turn)
	predict double se_bench, stdp
	
	FixAreg
	_vassert se se_bench, tol(1e-8)
	
cd "D:/Github/reghdfe/test"
exit
