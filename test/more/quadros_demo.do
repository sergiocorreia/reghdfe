clear all

cap ado uninstall ftools
net install ftools, from("C:/git/ftools/src")
ftools, compile

cap ado uninstall fereg
net install fereg , from(C:/git/reghdfe/src)
fereg, compile

sysuse auto
fereg price weight, a(turn trunk) timeit verbose(3) dof(pairwise)
fereg price weight, a(trunk turn) timeit verbose(3) dof(pairwise)
*/
//fereg price weight, a(trunk turn foreign) timeit verbose(3) dof(pairwise) noprune



use  "C:\Dropbox\Projects\HDFE\datasets\quadros\QP_Sergio.dta", clear
//use if year <= 2005 & firm<=39000 using "C:\Dropbox\Projects\HDFE\datasets\quadros\QP_Sergio.dta", clear


set seed 4543534
keep if runiform()<0.10

//use if year <= 1991 & firm<=35000 using "C:\Dropbox\Projects\HDFE\datasets\quadros\QP_Sergio.dta", clear
//keep if mod(worker, 7)
//keep if mod(worker, 11)
//keep if mod(firm, 7)
//keep if mod(firm, 11)
//save prune4, replace
//asd
//use prune4, clear

/*
foreach seed of global megaseeds  {
	set seed `seed'
	keep if runiform()<0.80
}

cou
loc N = c(N)
set seed $megaseed
keep if runiform()<0.90
if (c(N)==`N') exit

gen y = 0
reghdfe y, a(firm worker) dof(none) tol(1e-1)
gen byte smpl = e(sample)
keep if smpl
drop smpl y
if (c(N)==0) exit
// contract worker firm
*/



set seed 1234
gen double x = rt(4)
gen double y = 2000 * rt(4) + 10 * x

//collapse (mean) y x, by(worker firm)
//expand 2

// egen long id1 = group(worker)
// egen long id2 = group(firm)
// sort id1 id2
// order id1 id2 y x
// keep id1 id2 y x
// saveold c:\git\fereg\test\prune, replace
// exit


//clear mata
set rmsg off
set trace off
set more off
//cls
//mata: mata desc using lreghdfe

replace firm = firm*100

mata: mata desc

//clear
//use c:\git\mini
//rename id1 firm
//rename id2 worker


timer clear
loc tol 1e-12
loc tol 1e-8

/*
reghdfe y x, a(firm worker) dof(none) tol(1e-1)
gen byte smpl = e(sample)
keep if smpl
drop smpl
*/

local absorb worker firm // year jobtitle firm  year
sort `absorb'
timer on 99
fereg   y x, a(`absorb') timeit verbose(3) dof(pairwise)  tol(`tol')  // accel(none) transform(cimmino)
timer off 99

matrix list e(b)
matrix b = e(b)
scalar b1 = b[1,1]

//timer list
//exit

timer on 98
fereg   y x, a(`absorb') timeit verbose(3) dof(pairwise)  tol(`tol')  noprune
//reghdfe y x, a(`absorb') timeit verbose(3) fast dof(pairwise) tol(`tol')
timer off 98
matrix list e(b)
matrix b = e(b)
scalar b2 = b[1,1]

/*reghdfe y x, a(firm worker, save) timeit verbose(3) fast dof(pairwise) tol(`tol') transform(cimm)
predict double resid2, resid
matrix list e(b)
matrix b = e(b)
scalar b3 = b[1,1]
*/

// reghdfe y x, a(firm worker) timeit verbose(3) fast dof(none) transform(cimm)

timer list
di r(t99)/r(t98) * 100 // 23%
di b1
di b2

global mufasa 0
if (abs(b1-b2)>0.008) & (c(N)<$initial_obs) {
	global mufasa 1
}
//di b3
asd
exit

