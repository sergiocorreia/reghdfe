* Put this in parse if the option is needed
cap qui which tuples
assert inlist() => el rc debe ser 0 o 111
si es != 0 , dar un msg de instalar ssc install tuples , y hacer un exit 111 o algo asi

clear all
cls
sysuse auto, clear
keep if rep<.

cd "D:/github/reghdfe/source/_reghdfe_absorb"

local indepvars weight length
* Compute D := X'X (the bread of the sandwich)
reg price `indepvars', mse1
* Store the residuals to compute the meat of the sandwich
predict double resid, resid

local clustervars turn rep foreign trunk

* Notation: V = DMD
matrix D = e(V)
matrix backupD = D
local K = rowsof(D)
local minus = `K'
matrix ans = J(`K', `K', 0)

local tmp : word count `indepvars'
di "(`tmp'+1==`minus')"
assert (`tmp'+1==`minus')

* This gives all the required combinations (ssc install tuples)
tuples `clustervars' // create locals i) ntuples, ii) tuple1 .. tuple#

forval i = 1/`ntuples' {
	local vars `tuple`i''
	local numvars : word count `vars'
	local sign = cond(mod(`numvars', 2), "+", "-") // + with odd number of variables, - with even

	GenerateID `vars', gen(group)
	
	if (`numvars'==1) {
		su group, mean
		di as text "`vars' has `r(max)' groups"		
	}
	
	* Compute the full sandwich
	_robust resid, variance(D) minus(`minus') cluster(group)
	di as error "`sign' `vars'"

	* Add it to the other sandwiches
	matrix ans = ans `sign' D
	matrix D = backupD
	drop group
}

matrix list ans

* If the VCV matrix is not positive-semidefinite, use the fix from
* Cameron, Gelbach & Miller - Robust Inference with Multi-way Clustering (JBES 2011)
* 1) Use eigendecomposition V = U Lambda U' where U are the eigenvectors and Lambda = diag(eigenvalues)
* 2) Replace negative eigenvalues into zero and obtain FixedLambda
* 3) Recover FixedV = U * FixedLambda * U'
* This will fail if V is not symmetric (we could use -mata makesymmetric- to deal with numerical precision errors)
mata:
	ans = st_matrix("ans")
	if (!issymmetric(ans)) exit(error(505))
	symeigensystem(ans, U=., lambda=.) 
	st_local("eigenfix", "0")
	if (min(lambda)<0) {
		lambda = lambda :* (lambda :>= 0)
		// ans = U * diag(lambda) * U'
		ans = quadcross(U', lambda, U')
		st_local("eigenfix", "1")
	}
	st_replacematrix("ans", ans)
end

matrix list ans

matrix b = e(b)
local n = e(N)
local dof = `n' - `K'
matrix asd = ans
ereturn post b asd // .00027968, dof(`dof')
ereturn display

cd "D:/tmp"
cgmreg price `indepvars', cluster(`clustervars')

di as text "ANS"
matrix list ans

di as text "CGMREG"
matrix list e(V)

di as text "DELTA"
matrix cgmreg = e(V)
matrix delta = ans - cgmreg
matrix list delta

di as text "RAWCOVMAT"
matrix list e(rawcovmat)

* TODO: sortpreserve b/c GenerateID messes up sorting
* TODO: Allow WEIGHTS with _robust
exit

REFS:

Samuel B. Thompson, Simple formulas for standard errors that cluster by both firm and time, Journal of Financial Economics, Volume 99, Issue 1, January 2011, Pages 1-10, ISSN 0304-405X, http://dx.doi.org/10.1016/j.jfineco.2010.08.016.
(http://www.sciencedirect.com/science/article/pii/S0304405X10001923)
Abstract: When estimating finance panel regressions, it is common practice to adjust standard errors for correlation either across firms or across time. These procedures are valid only if the residuals are correlated either across time or across firms, but not across both. This paper shows that it is very easy to calculate standard errors that are robust to simultaneous correlation along two dimensions, such as firms and time. The covariance estimator is equal to the estimator that clusters by firm, plus the estimator that clusters by time, minus the usual heteroskedasticity-robust ordinary least squares (OLS) covariance matrix. Any statistical package with a clustering command can be used to easily calculate these standard errors.
Keywords: Cluster standard errors; Panel data; Finance panel data



Save results in
e(clustvar) // NOT clustvarS
e(N_clust_VARNAME)
e(cluster_eigenfix) 0 / 1
e(N_clust) -> EL MINIMO

Tenemos e(N_hdfe) , asi q hacer e(N_clustvars) ??
e(N)


Como hacer el small-sample scaling?

i) Cada -D- es multiplicado por M/M-1 de ese D (i.e. r(max) dp del sum)
ii) Tb: solo hacer una vez al final, donde M = min(G_i) solo para las unidades, no las combinaciones de clusters

Aparte de eso ajustar segun K
