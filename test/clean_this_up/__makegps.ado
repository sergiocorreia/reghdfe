*************************************************************
* __MAKEGPS
*************************************************************
/* This routine is from Amine Ouazad's a2reg program */
/* It establishes the connected groups in the data */
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
keep `id1' `id2' `groupid'
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
