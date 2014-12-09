* Look for runtime errors in -reghdfe-

*cap ado uninstall reghdfe
*cd D:\Dropbox\Projects\stata\hdfe\code

* Config
	set more off
	clear all
	set trace off
	set tracedepth 5
	set traceexpand on
	cap cls
	set maxvar 5000
	timer clear
	
* Simple dataset
	sysuse auto
	gen int t = _n
	tsset t
	
* Battery of commands

* Absorb
	reghdfe price weight, a(foreign)
	reghdfe price weight, a(foreign turn)
	reghdfe price weight, a(foreign turn rep)
	reghdfe price weight, a(foreign#turn rep)
	reghdfe price weight, a(foreign#turn)
	reghdfe price weight, a(foreign#turn#rep)
	reghdfe price weight, a(foreign turn#rep)
	reghdfe price weight, a(foreign foreign#c.length)
	reghdfe price weight, a(foreign##c.length)
	
* Avge
	reghdfe price weight, a(foreign) avge(turn)
	reghdfe price weight, a(foreign) avge(turn#rep)
	cap reghdfe price weight, a(foreign) avge(turn#rep#c.length) // can't do this
		assert _rc!=0

* VCE
	reghdfe price weight, a(foreign) vce(robust)
	reghdfe price weight, a(foreign) vce(unadj)
	reghdfe price weight, a(foreign) vce(cluster t) // Should be same as unadj!

* Cluster
	reghdfe price weight, a(foreign) vce(cluster foreign)
		assert e(fe_nested_in_cluster)==1
	reghdfe price weight, a(foreign turn) vce(cluster foreign)
		assert e(fe_nested_in_cluster)==1
	reghdfe price weight, a(foreign#turn) vce(cluster foreign#turn)
		assert e(fe_nested_in_cluster)==1
	reghdfe price weight, a(foreign) vce(cluster foreign#turn)
		assert e(fe_nested_in_cluster)==0
		
* Save FEs
	reghdfe price weight, a(A=foreign##c.length B=turn) avge(C=mpg)
		conf var A A_slope B C
		drop A A_slope B C
		
* e(sample)
	reghdfe price weight, a(turn)
		cou if e(sample)
		assert r(N)==74
	reghdfe price weight, a(rep)
		cou if e(sample)
		assert r(N)==69

* IV - IVREG2
	reghdfe price (weight=length) , a(foreign) vce(unadj)
	reghdfe price (weight=length) , a(foreign) vce(unadj) dofminus(large)
	
	reghdfe price (weight=length) , a(foreign) vce(robust)
	reghdfe price (weight=length) , a(foreign) vce(cluster t) // Not the same as vce(unadj) due to V adjustements
	reghdfe price (weight=length) , a(foreign) vce(cluster t) dofminus(large)
	reghdfe price (weight=length) , a(mpg) vce(cluster mpg) ivsuite(ivreg2) dofminus(large)
	reghdfe price (weight=length) , a(mpg) vce(cluster mpg) ivsuite(ivreg2) dofminus(small)
	
* IV - IVREGRESS
	reghdfe price (weight=length) , a(foreign) vce(unadj)  ivsuite(ivregress)

* Time series and factors

	reghdfe S.price, a(foreign) verbose(3)
	reghdfe S.price i.L.turn, a(foreign) verbose(3)
	reghdfe L.price F(1/2).weight S.length i.turn, a(foreign)
	reghdfe L.price (weight=length) , a(foreign)
	reghdfe L.price i.turn (S.length = F(1/2).weight), a(foreign)
	
	replace t = .
	tsset t
	reghdfe price weight, a(foreign) // Should run since we are not using tsset; prevents reversion of M.Manacorda's bug
	replace t = _n


* Options
set trace off
	reghdfe price weight, a(foreign) verbose(3)
 	reghdfe price weight, a(foreign) fast // verbose(3)
	reghdfe price weight, a(A=foreign) fast
	drop A
	reghdfe price weight, a(foreign) nested
	reghdfe price weight, a(foreign) check
	reghdfe price weight, a(foreign) nested check
	reghdfe price weight, a(foreign) avge(turn) excludeself
	reghdfe price weight, a(foreign) group(G)
	reghdfe price weight, a(foreign turn) group(G)
	drop G

* IV Options
	reghdfe price (weight=length) , a(foreign) first
	estimates dir
	reghdfe price (weight=length) , a(foreign) first savefirst
	estimates dir
	reghdfe price (weight=length) , a(foreign) suboptions(nocollin) level(90)
	reghdfe price weight length , a(foreign) verbose(3)
	reghdfe price (weight=length) , a(foreign) suboptions(coeflegend) ivsuite(ivregress)
	
* Maxim options
	reghdfe price weight, a(foreign turn rep) tol(1e-12) maxiter(0) noaccel verbose(1)
	reghdfe price weight, a(foreign turn rep) accel_start(10) verbose(3)
	reghdfe price weight, a(foreign turn rep) accel_freq(2) verbose(3)
	
exit


These should be the same:
tab mpg, gen(FE)
cls
ivreg2 price FE* (weight=length) , partial(FE*) cluster(mpg) small

reghdfe price (weight=length) , a(mpg) vce(cluster mpg) ivsuite(ivreg2) dofminus(small) dofmethod(naive)
