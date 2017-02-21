cls
clear all
discard
set more off
set trace off

cap ado uninstall ftools
net install ftools, from("C:/git/ftools/src")
ftools, check
ftools, compile

cap ado uninstall reghdfe
net install reghdfe , from(C:/git/reghdfe/src)
//set trace on
reghdfe, compile

clear all


cscript eret
sysuse auto
gen GEAR = gear
gen w = 10
loc w [aw=w]
loc w
loc cluster trunk turn
reghdfe price weight gear GEAR  `w', a(turn) keepsing old cluster(`cluster')
storedresults save old e()
reghdfe price weight gear GEAR  `w', a(turn) keepsing cluster(`cluster')
storedresults compare old e() , ///
	tol(1e-12) ///
	exclude(macros: vcesuite footnote estat_cmd predict dofadjustments subcmd cmdline marginsok scalar: savestages M1 K1 M1_exact M1_nested G1 unclustered_df_r F_absorb) 

exit


// whatabout number of iterations?!
