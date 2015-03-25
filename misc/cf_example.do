// -------------------------------------------------------------------------------------------------
// Murphy-Topel: Asymptotic VCE estimator for two-step models
// -------------------------------------------------------------------------------------------------
/// link is (minus sign) derivative of RESID2 wrt beta1
/// Equivalently, is the derivative of the INDEX FUNCTION x2b2 wrt b1
// In the OLS case, it's b2[zhat] * x1

local Varlist 		string scalar
*local Integer 		real scalar
*local VarByFE 		real colvector // Should be levels*1
local Series		real colvector // Should be N*1
local Matrix		real matrix
*local SharedData 	external struct FixedEffect vector

clear all

mata:
mata set matastrict on
void function murphy_topel(`Varlist' x1_name, `Varlist' x2_name, `Varlist' s1_name, `Varlist' s2_name,
	`Varlist' link_name, string scalar v1_name, string scalar v2_name,
	`Varlist' touse, real scalar hasconstant) {

	`Matrix' X1, X2, V1, V2, C, R, MTVAR
	real colvector s1, s2, link


	printf("{txt}(calculating Murphy-Topel Variance)\n")
	st_view(X1, ., x1_name, touse)
	st_view(X2, ., x2_name, touse)
	st_view(s1, ., s1_name, touse)
	st_view(s2, ., s2_name, touse)
	st_view(link, ., link_name, touse)
	V1 = st_matrix(v1_name)
	V2 = st_matrix(v2_name)


	C = quadcross(X2, hasconstant, s2:^2 :* link, X1, hasconstant)
	R = quadcross(X2, hasconstant, s2 :* s1, X1, hasconstant)

	MTVAR = V2 + V2 * (C*V1*C' - R*V1*C' - C*V1*R') * V2
	//mm_matlist(MTVAR)
	//mm_matlist(invsym(MTVAR))
	st_replacematrix(v2_name, MTVAR)
}
end

capture program drop murphy_topel
program define murphy_topel, eclass
	syntax, x1(varlist numeric) x2(varlist numeric) ///
		s1(varname numeric) s2(varname numeric) ///
		link(varname numeric) ///
		var1(string) var2(string) ///
		sample(varname) ///
		[cluster(varname)]

	local hasconstant = 1

	if ("`cluster'"!="") {
		preserve
		local uniquevarlist `x1' `x2' `s1' `s2' `link'
		local uniquevarlist : list uniq uniquevarlist
		keep if `sample'
		collapse (sum) `uniquevarlist' (max) `sample', by(`cluster') fast
	}
	de

	* Will overwrite V2
	mata: murphy_topel("`x1'", "`x2'", "`s1'", "`s2'", "`link'", "`var1'", "`var2'", ///
		"`sample'", `hasconstant')

	if ("`cluster'"!="") {
		restore
	}

	ereturn repost V = `var2'
	*st_matrix(name)
	*st_replacematrix(name, X) --> REEMPLAZAR V2
	*st_view(X, ??, "var1 var2 ..", "touse")
	*usar tokens() pa las vars
end

* Playing with control functions
	clear
	cls
	set more off
	sysuse auto

* Variables
	local depvar price
	local endogvar head
	local exogvars weight length
	local instruments gear disp foreign rep

* Fake dataset
	clear
	set obs 1000000
	gen int group = 1 + int((_n-1)/1000)
	su group
	sort group
	cou if group==0
	cou if group==1
	cou if group==1000
	cou if group==1001
	by group: gen fe1 = rnormal() // * 0
	by group: replace fe1 = fe1[1]
	by group: gen fe2 = rnormal() // * 0
	by group: replace fe2 = fe2[1]

	foreach var in `exogvars' `instruments' error1 error2 {
		gen double `var' = rnormal()
	}
	gen `endogvar' = 1*gear + 1*disp + 1*foreign + 1*rep + 1*length + 15*error1 + 10*error2 + 20 * fe1
	gen `depvar' = 100*`endogvar' + 100*weight + 100*length + 1000*error2 + 20 * fe2

* OLS (wrong)
	regress `depvar' `endogvar' `exogvars'
	local tss = e(rss) + e(mss)

* Manual First Stage
	regress `endogvar' `exogvars' `instruments'
	matrix V1 = (e(df_r)/e(N)) * e(V)
	
	predict double s1, scores // resid
	local mse = e(rss)/e(N)
	replace s1 = s1 / `mse'
	
	predict double zhat, xb

* "Naive Covariance Matrix"
	regress `depvar' zhat `exogvars'
	*regress `depvar' zhat `exogvars', robust
	
* IV Version
	*regress `depvar' zhat `exogvars' // , robust
	matrix V2 = (e(df_r)/e(N)) * e(V)
	gen byte smpl = e(sample)
	
	predict double s2, scores // resid
	local mse = e(rss)/e(N)
	replace s2 = s2 / `mse'

	gen double link = _b[zhat]

* Compute Corrected VCE
	murphy_topel, x1(`exogvars' `instruments') x2(zhat `exogvars') s1(s1) s2(s2) link(link) ///
		var1(V1) var2(V2) sample(smpl) // DOESNT WORK AT ALL, code is wrong cluster(group)
	estadd local vcetype "Murphy Topel", replace
	regress // Replay correct

* Benchmark
	ivreg2 `depvar' `exogvars' (`endogvar' = `instruments') , small
	ivreg2 `depvar' `exogvars' (`endogvar' = `instruments') , small robust
	ivreg2 `depvar' `exogvars' (`endogvar' = `instruments') , small cluster(group)
	*ivregress 2sls `depvar' `exogvars' (`endogvar' = `instruments') , small robust




exit
/*
M-T Variance = V2 + V2 W V2 where
W = C V1 C' - R V1 C' - C V1 R'

C = sum(1..N) of deriv(logl2i, beta2) * deriv(logl2i, beta1)'
R = sum(1..N) of deriv(logl2i, beta2) * deriv(logl1i, beta1)'

C = sum s_22i s_21i' = S22 S21
R = sum s_22i s_11i' = S22 S11

Since
logl = .. - 0.5 log sigma^2 - 1/(2 sigma) resid^2 ; resid=y - x beta

Then, calling z the depvar of eqn1 and zhat its predicted value,
and ztilde the derivative of zhat wrt the index function xb
(which is 1 in a linear model)

Also, call v = resid / mse ; where mse = (resid'resid/n)

S22 = X2 v2
S11 = X1 v1
S21 = (-resid/sigma) * resid(Theta)' = -v2 * resid(Zhat) * Zhat(Theta) = v2 * _b[Zhat] * X1
    = X1 (v2 _b[Zhat])
WRONG!!: zhat v2 ztilde X1 = X1 (v2 zhat ztilde)

Cross notation:
cross(X, xc, w, Z, zc) = X'diag(w)Z
where X/Z are augmented on the right with a col of 1s if xc/zc!=0
Simplifies to cross(X,w,Z)

Then, assuming constant was excluded from Xs
C = cross(X2, 1, v2^2 zhat ztilde, X1, 1)
R = cross(X2, 1, v2 v1, X1, 1)
