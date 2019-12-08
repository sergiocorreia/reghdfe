cap ado uninstall reghdfe
loc dir `c(pwd)'
loc dir : subinstr loc dir "test" "src"
di as text "DIR: `dir'"
net install reghdfe, from("`dir'")

cap findfile "lreghdfe.mlib"
if (!c(rc)) cap drop "`r(fn)'"


clear all
cls

loc path "C:\Git\asd\borrar"


sysuse auto
reghdfe price weight length, a(turn)
exit
set trace off
reghdfe price weight length, a(turn) verbose(1) parallel(2, tmp("C:\Git\asd\borrar"))



exit

mata:
fn = "C:\Git\asd\borrar\PARALLEL_286566349\data0.tmp"
fh = fopen(fn, "r")
//hdfe = fgetmatrix(fh)
//fclose(fh)
end
