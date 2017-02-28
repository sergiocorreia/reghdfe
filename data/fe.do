clear all
set more off
set maxvar 32767

/*
* Install
ssc install reg2hdfe
net install st0143_3.pkg
net get st0143_3.pkg
ssc install a2reg
*/

* (1) Ancilliary dataset from felsdvreg -> This is different
sysuse felsdvsimul

twfe y x1 x2 , id(i j) tol(1e-10)
areg y x1 x2 i.j , absorb(i)
reg2hdfe y x1 x2 , id1(i) id2(j) tol(1e-8)
a2reg y x1 x2, individual(i) unit(j)
felsdvreg y x1 x2 , ivar(i) jvar(j) peff(fe1) feff(fe2) xb(xb) res(res) mover(mover) mnum(mnum) pobs(pobs) group(group)
	drop fe1 fe2 xb res mover mnum pobs group

exit


* (2) Real example -> This works
use callao_tmp, clear
qui tab j1 , gen(fe2)
qui tab j2 , gen(fe3)
compress
	
areg y x1 x2 fe3*  , absorb(i)
reg2hdfe y x1 x2 , id1(i) id2(j2) tol(1e-8) maxiter(1000)
twfe y x1 x2 , id(i j2) tol(1e-10)
* DIES: a2reg y x1 x2, individual(i) unit(j2)

exit


