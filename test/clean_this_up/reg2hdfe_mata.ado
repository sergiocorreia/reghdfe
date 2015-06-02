*! version 2.01 07oct2011
* Estimates linear regression model with two high dimensional fixed effects 
/*---------------------------------------------------------*/
/* Guimaraes & Portugal Algorithm */
/* Author: Paulo Guimaraes */
/*---------------------------------------------------------*/
/* Based on: */
/* Paulo Guimaraes and Pedro Portugal. "A Simple Feasible Alternative Procedure to Estimate Models with */
/* High-Dimensional Fixed Effects", Stata Journal, 10(4), 628-649, 2010. */

include fastmean.mata // MODIFIED

program reg2hdfe_mata, eclass
version 9.1
if replay() {
if ("`e(cmd)'"!="reg2hdfe") error 301
Display `0'
}
else Estimate `0'
cap drop _id2
end

program define Estimate, eclass
syntax varlist [if] [in], id1(str) id2(str)  ///
[TOL(real 0.000001) MAXiter(integer 0) ///
CHECK NODOTS SIMPLE fe1(str) fe2(str) cluster(str) GROUPid(str)  INdata(string) ///
OUTdata(string) IMProve(str) VERBose NOREGress PARAM1 PARAM2 OP1(integer 3) ///
OP2(integer 1) OP3(integer 10) OP4(integer 1000) OP5(real 0.001) AUTOFF NESTEDFE] // MODIFIED

*********************************************************************
* Checking syntax
*********************************************************************
tokenize `varlist'
local lhs `1'
mac shift
local rhs `*'

if ("`fe1'"!=""&"`fe2'"=="")|("`fe2'"!=""&"`fe1'"=="") {
di in red "Error: You must specify both options fe1 and fe2"
error 198
}

if "`param1'"=="param1"&"`param2'"=="param2" {
di in red "Error: Choose either param1 or param2"
error 198
}

if "`indata'"!=""&"`improve'"!="" {
di in red "Error: Indata option not valid with improve option"
error 198
}

if "`outdata'"!=""&"`improve'"!="" {
di in red "Error: Outdata option not valid with improve option"
error 198
}

if "`improve'"!=""&"`rhs'"!="" {
di in red "Error: Improve option can only be used with a single variable"
error 198
}

if "`improve'"!=""&"`fe1'"!="" {
di in red "Error: Can not estimate fixed effects with improve option"
error 198
}

if "`improve'"!=""&"`check'"!="" {
di in red "Error: Can not use check and improve option simultaneously"
error 198
}

if "`improve'"!=""&"`cluster'"!="" {
di in red "Error: Can not use cluster and improve option simultaneously"
error 198
}

if "`indata'"==""&"`improve'"==""&"`outdata'"=="" {
local standard "standard"
}

if `"`fe1'"'!=`""' confirm new var `fe1' 
if `"`fe2'"'!=`""' confirm new var `fe2'
if `"`groupid'"'!=`""' confirm new var `groupid'

capture drop __uid
**********************************************************************
* Define Initial Variables
**********************************************************************
tempvar clustervar

di in ye "=============================================================="

local dots `=cond("`nodots'"=="",1,0)'

***********************************************************************
* Do Main Loop
***********************************************************************

if "`improve'"!="" {
preserve 
tempvar varfe2
tempfile tmp1 tmp2
quietly {
if "`verbose'"!="" {        
noisily di "Reading `improve'_ids"
}
use `improve'_ids, clear
sort __uid
qui save `tmp1', replace
merge __uid using `improve'_`lhs'
sum _merge, meanonly
if r(min)<r(max) { 
di "There was an error merging `improve'_ids with `improve'_`lhs'"
di "Data sets do not match"
error 198
}
drop _merge
rename __t_`lhs' `lhs'
}
* Now we try to improve convergence
di in ye "Improving Convergence for Variable: `lhs'" 
di in red "Improve may not work if fixed effects are not specified in the same order as saved"
conf var _id2
iteralg "`lhs'" "`id1'" "`id2'" "__fe2_`lhs'" "`tol'" "`maxiter'" "`simple'" ///
"`verbose'" "`dots'" "`varfe2'" "`op1'" "`op2'" "`op3'" "`op4'" "`op5'" "`autoff'"
if "`outdata'"!="" {
outvars "__o_`lhs'" "`lhs'" "`varfe2'" "`outdata'" 
}
if "`check'"=="check" {
checkvars "__o_`lhs'" "`varfe2'" "`id1'"
}
restore
}     

if "`indata'"!="" {      
tempfile tmp1 tmp2 tmp3 readdata
quietly {
if "`verbose'"!="" {        
noisily di "Reading `indata'_ids"
}
use `indata'_ids, clear
sort __uid
qui save `tmp1', replace
if "`cluster'"!="" {
noisily di in yellow "The clustering variable is the one used when the data was created!!! "
if "`verbose'"!="" {        
noisily di in yellow "Adding `indata'_clustervar"
}
merge __uid using `indata'_clustervar
sum _merge, meanonly
if r(min)<r(max) { 
di "There was an error merging `indata'_ids with `indata'_clustervar "
di "Data sets do not match"
error 198
}
drop _merge
sort __uid
rename __clustervar `clustervar'
qui save `tmp1', replace
}
* Now read the original variables
foreach var in `varlist' {
if "`verbose'"!="" {        
noisily di "Adding original `indata'_`var'"
}
merge __uid using `indata'_`var'
sum _merge, meanonly
if r(min)<r(max) { 
di "There was an error merging `indata'_ids with `indata'_`var' "
di "Data sets do not match"
error 198
}
drop _merge
drop __fe2*
drop __t_*
sort __uid
qui save `tmp2', replace
}
foreach var in `lhs' `rhs' {
rename __o_`var' `var'
}
sum `lhs', meanonly
tempvar yy sy
gen double `yy'=(`lhs'-r(mean))^2
gen double `sy'=sum(`yy')
local tss=`sy'[_N]
drop `yy' `sy'
qui save `readdata'
use `tmp1', clear
foreach var in `varlist' {
if "`verbose'"!="" {        
noisily di "Adding transformed `indata'_`var'"
}
merge __uid using `indata'_`var'
sum _merge, meanonly
if r(min)<r(max) { 
di "There was an error merging `indata'_ids with `indata'_`var' "
di "Data sets do not match"
error 198
}
drop _merge
drop __fe2*
drop __o_*
sort __uid
qui save `tmp3', replace
}
        
foreach var in `lhs' `rhs' {
rename __t_`var' `var'
}
if "`verbose'"!="" {
noisily di "Done reading data"
}
}
}   

if "`outdata'"!=""|"`standard'"!="" {

************ Mark usable sample and store all data
mark __touse `if' `in'
if "`cluster'"!="" {
gen double `clustervar'=`cluster'
markout __touse `lhs' `rhs' `id1' `id2' `clustervar'
}
else {
markout __touse `lhs' `rhs' `id1' `id2'
}

tempfile origdata
if "`verbose'"!="" {
di in yellow "Saving original file"
}
gen long __uid = _n
sort __uid

* MODIFIED (just in case ID2 is not 1,2,3..G)
egen long _id2 = group(`id2') if __touse

qui save `origdata'

* Restrict data to usable sample
qui keep if __touse
if "`cluster'"!="" {
keep __uid `id1' `id2' `lhs' `rhs' `clustervar' _id2
}
else {
keep __uid `id1' `id2' `lhs' `rhs' _id2
}


*************************
sum `lhs', meanonly
tempvar yy sy
gen double `yy'=(`lhs'-r(mean))^2
gen double `sy'=sum(`yy')
local tss=`sy'[_N]
drop `yy' `sy'

***************************
if "`outdata'"!="" {
preserve
keep __uid `id1' `id2' _id2 // MODIFIED
order __uid `id1' `id2'
if "`verbose'"!="" {
di in yellow "Saving fixed effects variables... "
}
capture qui save `outdata'_ids
local err _rc
if `err'!=0 {
use `origdata', clear
drop __touse
error `err'
}
restore
if "`cluster'"!="" {
if "`verbose'"!="" {
di in yellow "Saving cluster variable... "
}
preserve
keep __uid `clustervar'
rename `clustervar' __clustervar
qui save `outdata'_clustervar
restore
}        
}                         
di in ye "Tolerance Level for Iterations: `tol' "
tempvar start2
tempname dum1
gen double `start2'=0
foreach var of varlist `varlist' {
di
di in ye "Transforming variable: `var' " 
gen double __o_`var' = `var'
conf var _id2
iteralg "`var'" "`id1'" "`id2'" "`start2'" "`tol'" "`maxiter'" "`simple'" ///
"`verbose'" "`dots'" "`dum1'" "`op1'" "`op2'" "`op3'" "`op4'" "`op5'" "`autoff'"
if "`outdata'"!="" {
outvars "__o_`var'" "`var'" "`dum1'" "`outdata'" 
}
checkvars "__o_`var'" "`dum1'" "`id1'"
drop __o_`var'
drop `dum1'
}
drop `start2'   
}  

if "`improve'"==""&"`noregress'"=="" {
if "`verbose'"!="" {        
di "Calculating degrees of freedom ... "
}

* Create group variable
tempvar group
qui __makegps, id1(`id1') id2(`id2') groupid(`group')

* Calculate Degrees of Freedom	
qui count
local N = r(N)
local k : word count `rhs'
sort `id1'
qui count if `id1'!=`id1'[_n-1]
local G1 = r(N)
sort `id2'
qui count if `id2'!=`id2'[_n-1]
local G2 = r(N)
sort `group'
qui count if `group'!=`group'[_n-1]
local M = r(N)
local kk = `k' + `G1' + `G2' - `M'
local dof = `N' - `kk'	

* Estimate the model
if "`verbose'"!="" {        
            di "Estimating Regression ... "
}

tempname name1 name2

if "`cluster'"=="" {
di
* Estimate Regression		
qui _regress `lhs' `rhs', nocons dof(`dof')
estimates store `name2'
local r=1-e(rss)/`tss'
ereturn scalar df_m = `kk'-1
ereturn scalar mss=`tss'-e(rss)
ereturn scalar r2=`r'
ereturn scalar r2_a=1-(e(rss)/e(df_r))/(`tss'/(e(N)-1))
ereturn scalar F=(`r'/(1-`r'))*(e(df_r)/(`kk'-1))
ereturn local cmdline "reg2hdfe `0'"
ereturn local cmd "reg2hdfe"
ereturn local predict ""
ereturn local estat_cmd ""
estimates store `name1'
}
else {
sort `clustervar'
qui count if `clustervar'!=`clustervar'[_n-1]
local Nclust = r(N)
qui _regress `lhs' `rhs', nocons mse1
estimates store `name2'
tempname b V
matrix `V'=e(V)
matrix `b'=e(b)
local rss=e(rss)
local r=1-`rss'/`tss'
local nobs=e(N)
tempvar res
predict double `res', residual


* MODIFIED
if "`cluster'"!="" & "`nestedfe'"!="" {
	*** Relevant? est scalar df_r = min(`df_r', `df_cl')
	
	* Replaced G1 with 1
	*local adj_kk = `k' + 1 + `G2' - `M'
	local adj_kk = `kk' - `G1' + 1
	local adj_dof = `G1' - 1
	
	_robust `res', v(`V') minus(`adj_kk') cluster(`clustervar')
	ereturn scalar Mgroups = `M'
	ereturn post `b' `V', depname(`lhs') obs(`nobs') dof(`adj_dof')
}
else {
	_robust `res', v(`V') minus(`kk') cluster(`clustervar')
	ereturn scalar Mgroups = `M'
	ereturn post `b' `V', depname(`lhs') obs(`nobs') dof(`kk')
}


ereturn local eclustvar "`cluster'"
ereturn local vce "cluster"
ereturn local vcetype "Robust"
ereturn local cmdline "reg2hdfe `0'"
ereturn local depvar "y"
ereturn local cmd "reg2hdfe"
ereturn scalar N_clust=`Nclust'
ereturn scalar r2=`r'
ereturn scalar rss=`rss'
ereturn scalar mss=`tss'-`rss'
estimates store `name1'

* MODIFIED
ereturn scalar G1 = `G1'
ereturn scalar G2 = `G2'
ereturn scalar M = `M'
assert "`G1'"!=""

di
}

if "`verbose'"!="" {        
di "Done with estimation ... "
}
}

if "`indata'"!="" {
use `readdata', clear
}
if "`indata'"==""&"`improve'"=="" {
use `origdata', clear
quietly keep if __touse==1
drop __touse
}
conf var _id2
* Compute Fixed Effects

if `"`fe1'"'!=`""' & `"`fe2'"'!=`""' { 
di in ye "Calculating Fixed Effects"
tempvar dum1 dum2 dum3
qui estimates restore `name2'
qui predict `dum1', res
qui gen double `dum2'=`dum1'     
*di in ye "Tolerance Level for Iterations: `tol'"
local dots `=cond("`nodots'"=="",1,0)'
tempvar start2
gen double `start2'=0
quietly {
conf var _id2
iteralg "`dum1'" "`id1'" "`id2'" "`start2'" "`tol'" "`maxiter'" "`simple'" ///
"`verbose'" "`dots'" "`dum3'" "`op1'" "`op2'" "`op3'" "`op4'" "`op5'" "`autoff'"
}
drop `start2'
qui replace `dum2'=`dum2'-`dum3', nopromote 
sort `id1' 
qui by `id1': g double `fe1' = sum(`dum2')/_n
qui by `id1': replace `fe1' = `fe1'[_N], nopromote
rename `dum3' `fe2'
di in ye "Done!!! "
* Implement parameterization 1
if "`param1'"=="param1" {
tempvar sumfe 
gen double `sumfe'=`fe1'+`fe2'
qui __makegps, id1(`id1') id2(`id2') groupid(`group')
sort `group' 
qui by `group': replace `dum1' = sum(`fe1')/_n, nopromote
qui by `group': replace `fe1' = `fe1'-`dum1'[_N], nopromote
qui replace `fe2'=`sumfe'-`fe1', nopromote
}
* Implement parameterization 2
if "`param2'"=="param2" {
tempvar sumfe 
gen double `sumfe'=`fe1'+`fe2'
qui __makegps, id1(`id1') id2(`id2') groupid(`group')
sort `group' 
qui by `group': replace `dum1' = sum(`fe2')/_n, nopromote
qui by `group': replace `fe2' = `fe2'-`dum1'[_N], nopromote
qui replace `fe1'=`sumfe'-`fe2', nopromote
}
Display `name1'
local nodisp nodisp
* Test Final Model 
if "`check'"=="check" {
qui _regress `lhs' `rhs' `fe1' `fe2' 
di
di in yellow "Checking if final model converged - Coefficients for fixed effects should equal 1"
di in yellow "Coefficient for `id1' --> "_b[`fe1']
di in yellow "Coefficient for `id2' --> "_b[`fe2']
}
}

if `"`groupid'"'!=`""' {
qui __makegps, id1(`id1') id2(`id2') groupid(`groupid')
label var `groupid' "Unique identifier for mobility groups"
}

if "`indata'"==""&"`improve'"=="" {
tempfile addvars
qui describe
if r(k)> 1 {
keep __uid `fe1' `fe2' `groupid'
sort __uid
qui save `addvars', replace 
use `origdata', clear
drop __touse
sort __uid
merge __uid using `addvars'
drop _merge
}
else {
use `origdata', clear
drop __touse
}
}

if "`nodisp'"!="nodisp" {
Display `name1'
}

capture drop __uid
di
end

program Display
args name
if "`name'"!="" {
qui estimates restore `name'
_coef_table_header, title( ********** Linear Regression with 2 High-Dimensional Fixed Effects ********** )
_coef_table
}
end

program define iteralg
args var id1 id2 start2 tol maxiter simple verbose dots jfe op1 op2 op3 op4 op5 autoff
conf var _id2
local count1=0
local count2=0
local count3=0
local count4=0
local count5=0
quietly {
recast double `var'
tempvar temp v1 v2 mean
gen double `v1'=0
gen double `temp'=0
gen double `v2'=`start2'
qui sum `v2', meanonly
if r(min)!=r(max) {
qui replace `var' = `var' + `v2', nopromote
}
if "`simple'"=="" {
tempvar v0 ym1 ym2 dum
gen double `v0'=0
gen double `ym1'=0
gen double `ym2'=0
}
local iter=1
local dif=1
local c1 "(`v2'>`v1'&`v1'>`v0'&((`v2'+`v0')<(2*`v1')))"
local c2 "(`v2'<`v1'&`v1'<`v0'&((`v2'+`v0')>(2*`v1')))"
capture drop `mean'
sort `id1'
by `id1': g double `mean' = sum(`var')/_n
qui by `id1': replace `var' = `var' - `mean'[_N], nopromote	
while abs(`dif')>`tol' & `iter'!=`maxiter'{
*while `dif'<1 & `iter'!=`maxiter'{

capture drop `mean'
sort `id1'
by `id1': g double `mean' = sum(`v2')/_n
qui by `id1': replace `mean' = `mean'[_N], nopromote				
if "`simple'"=="" {
capture drop `v0'
rename `v1' `v0'
}
else {
capture drop `v1'
}
rename `v2' `v1'

* MODIFIED
*sort `id2'
*by `id2': g double `v2' = sum(`var'+`mean')/_n 
*qui by `id2': replace `v2' = `v2'[_N], nopromote
gen double _tmp = `var'+`mean'
gen double `v2' = .
su `v2'
conf var _id2
mata: fastmean("_tmp", "`v2'", "_id2") // (var, newvar, groupid)
drop _tmp

if `iter'>`op1'&"`simple'"=="" {
capture drop `dum'
gen `dum'=`v1'+((`v2'-`v1')*(`v1'-`v0')/(2*`v1'-`v0'-`v2'))*(1-[(`v2'-`v1')/(`v1'-`v0')]^`op2')
replace `v2'=`dum' if (`c1'|`c2')&(`dum'<.), nopromote 
if mod(`iter',`op3')==0 {
replace `ym1'=`ym2', nopromote
replace `ym2'=`v1'+((`v2'-`v1')*(`v1'-`v0')/(2*`v1'-`v0'-`v2'))*(1-[(`v2'-`v1')/(`v1'-`v0')]^`op4')
replace `v2'=`ym2' if (`c1'|`c2')&(abs(`ym1'-`ym2')<`op5'), nopromote
}
}
qui replace `temp'=sum(reldif(`v2',`v1')), nopromote
local dif=`temp'[_N]/_N
*count if abs(`v2'-`v1')<`tol'
*local dif=r(N)/_N
if `dots' {
if `iter'==1 {
_dots 0, title(Iterations) reps(`maxiter')
}
_dots `iter' 0
}
if "`verbose'"!="" {
count if abs(`v2'-`v1')<`tol'
noisily di " `iter' - Dif --> " %-12.7g `dif' "  % fe below tolerance --> "  %-07.5f 100*r(N)/_N " OP2 -> "`op2'
*noisily di "`iter' --> % fe below tolerance --> "  %-08.5f 100*`dif' " OP2 -> "`op2'
} 
local count1=`count2'
local count2=`count3'
local count3=`count4'
local count4=`count5'
local count5=`dif'
if (`count5'>`count4')&(`count4'>`count3')&("`autoff'"=="")&("`simple'"=="")&(`op2'>1) {
local op2=1
local count5=0
}
if (`count5'<`count4')&(`count4'<`count3')&(`count3'<`count2')&(`count2'<`count1')&("`autoff'"=="")&("`simple'"=="") {
local op2=`op2'+1
local count5=.
}
local iter=`iter'+1
}
qui replace `var' = `var' - `v2' + `mean', nopromote
}
if `iter'==`maxiter' {
di
di in red "Maximum number of iterations reached"
di in red "Algorithm did not converge for variable `var'"
di in red "Last improvement: `dif'"
}			
else {			
di
di in yellow "Variable `var' converged after `iter' Iterations"
}
if "`jfe'"!="" {
gen double `jfe'=`v2'
}	
end

program define outvars
args orig var fe2 outdata
preserve
keep __uid `orig' `var' `fe2' 
sort __uid
rename `var' __t_`var'
rename `fe2' __fe2_`var'
qui save `outdata'_`var', replace
di in yellow " `var' was saved "
restore
end       

program define checkvars
args orig fe2 id1
tempvar fe1 dum2
gen double `dum2'=`orig'-`fe2'
sort `id1' 
by `id1': g double `fe1' = sum(`dum2')/_n
qui by `id1': replace `fe1' = `fe1'[_N], nopromote
qui _regress `orig' `fe1' `fe2' 
di in yellow "Checking if model converged - Coefficients for fixed effects should equal 1"
di in yellow "Coefficient for id1 --> "_b[`fe1']
di in yellow "Coefficient for id2 --> "_b[`fe2']
end

/* This routine is from Amine Ouazad's a2reg program */
/* It establishes the connected groups in the data */

*Find connected groups for normalization
capture program drop __makegps
program define __makegps
version 9.2
syntax [if] [in], id1(varname) id2(varname) groupid(name)
marksample touse
markout `touse' `id1' `id2'
confirm new variable `groupid'
sort `id1' `id2'
preserve
*Work with a subset of the data consisting of all id1-id2 combinations
keep if `touse'
collapse (sum) `touse', by(`id1' `id2')
sort `id1' `id2'
*Start by assigning the first id1 value to group 1, then iterate to fill this out
tempvar group newgroup1 newgroup2
gen double `group'=`id1'
local finished=0
local iter=1
while `finished'==0 {
quietly {
bysort `id2': egen double `newgroup1'=min(`group')
bysort `id1': egen double `newgroup2'=min(`newgroup1')
qui count if `newgroup2'~=`group'
local nchange=r(N)
local finished=(`nchange'==0)
replace `group'=`newgroup2'
 drop `newgroup1' `newgroup2'
}
di in yellow "On iteration `iter', changed `nchange' assignments"
local iter=`iter'+1
}
sort `group' `id1' `id2'
tempvar nobs complement
by `group': egen double `nobs'=sum(`touse')
replace `nobs'= -1*`nobs'
egen double `groupid'=group(`nobs' `group')
keep `id1' `id2' `groupid' // _id2 // MODIFIED
sort `id1' `id2'
tempfile gps
save `gps'
restore
tempvar mrg2group
merge `id1' `id2' using `gps', uniqusing _merge(`mrg2group')
assert `mrg2group'~=2
assert `groupid'<. if `mrg2group'==3
assert `groupid'==. if `mrg2group'==1
drop `mrg2group'
end  

