clear all
cls

cd "../test"
do setup.do
cd "../benchmarks"

clear all


cls
sysuse auto
reghdfe price weight length headroom trunk displacement, a(turn) tech(map) parallel(2, tmp(C:\Git\asd\borrar) cores(3)) v(1) // noisy



exit

use "../test/part2/toy-patents-long", clear
reghdfe citations funding lab_size, a(inventor_id year) group(patent_id) individual(inventor_id) agg(average) precond(diagonal) tech(lsmr) parallel(2, tmp(C:\Git\asd\borrar) cores(2)) v(1)



reghdfe price weight length headroom trunk displacement, a(turn) tech(lsmr) parallel(2, tmp(C:\Git\asd\borrar) cores(3)) v(1) // noisy

reghdfe price weight length headroom trunk displacement, a(turn) parallel(2, tmp(C:\Git\asd\borrar) cores(3)) // quiet
reghdfe price weight length headroom trunk displacement, a(turn) parallel(2, tmp(C:\Git\asd\borrar) cores(3) verbose) // show verbose of parallel only

reghdfe price weight length headroom trunk displacement, a(turn)

exit
