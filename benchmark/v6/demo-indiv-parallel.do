clear all
cls

cd "../test"
do setup.do
cd "../benchmarks"

clear all

use "../test/part2/toy-patents-long", clear

cls
reghdfe citations funding lab_size, a(inventor_id year) group(patent_id) individual(inventor_id) agg(average) precond(diagonal) tech(lsmr) parallel(2, tmp(C:\Git\asd\borrar) cores(2) method(procexec)) v(1)

reghdfe citations funding lab_size, a(inventor_id year) group(patent_id) individual(inventor_id) agg(average) precond(diagonal) tech(lsmr)
exit
