clear all

input long turn long trunk
1 1
2 1
2 2
3 2
4 2
4 3
5 2
5 3
6 3
6 4
6 5
6 5
6 5
6 5
end

li
set seed 123
gen weight = runiform()
gen price = runiform() * 50 + weight * 3
replace turn = turn * 100

cls


cap ado uninstall ftools
net install ftools, from("C:/git/ftools/src")
ftools, compile

cap ado uninstall fereg
net install fereg, from("C:/git/fereg/src")
fereg, compile


//sysuse auto, clear
//drop if turn==45 & trunk==13

/*
reghdfe price, a(turn trunk) 
gen byte smpl = e(sample)
keep if smpl
tab turn trunk, m
asd
*/
	
gen rep = int((_n-1)/7)
replace turn =  int(turn/500)
replace trunk = int(trunk/1.1)
gen z = 2

fereg price weight, a(turn rep trunk##c.z) v(3) keepsing 
fereg price weight, a(turn rep trunk##c.z) keepsing noprune
reghdfe price weight, a(turn rep trunk##c.z) keepsing


exit
reghdfe price weight, a(turn trunk) keepsing tol(1e-10)
