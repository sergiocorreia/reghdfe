clear all
pr drop _all

cap pr drop myreg
pr define myreg, eclass
	args N
	qui reg price weight
	matrix b = e(b)
	matrix V = e(V)
	ereturn clear
	local depvar = "price"
	local df_r 999
	ereturn post b 	V , dep(`depvar') obs(`N') dof(`df_r')
end

sysuse auto, clear
set more off
myreg 1000000
ereturn list

myreg 1000001
ereturn list

exit
