* How are robust standard error computed with 2sls?

* [Formula]
* VCV = q * V *S * V
* q = q = [ (N-1)/(N-K) ] * [ M / (M-1) ]
* V = X'Z (Z'Z)-1 Z'X
* S = Sum over M of U_k ' U_k where U_k are the score row vectors, thus
* S = Sum for all k over M of resid_k^2 Xhat_k' Xhat_k
* Xhat are PREDICTED not actual

* [Steps]
* Obtain Xhat from first stage (or from when betas are computed
* Run regression with -mse1- option, so we get V in e(V)
* Obtain resid "e"
* Calculate -q-
* Run _robust e, minus(..) cluster(..)
* Recall that minus is the K in N/N-K +-
* And dof(#) needs to be provided by ereturn post to compute t & F. Without those asymptotic distribs are assumed
* Basically, dof sets e(df_r)

clear all
sysuse auto
cls

*regress price weight length (weight displacement)
*ivregress 2sls price weight (length = displacement), small

reg length weight displacement
predict double lengthhat, xb

_regress price weight length  (weight displacement), mse1
predict double e, resid
*matrix D = e(V)
*matrix list D

rename length _length
rename lengthhat length
_robust e , minus(3) // v(D) 
rename length lengthhat
rename _length length
matrix list e(V)
estimates replay

ivreg2 price weight (length = displacement), robust small
matrix list e(V)

