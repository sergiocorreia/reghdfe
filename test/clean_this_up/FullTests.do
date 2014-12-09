// -------------------------------------------------------------
// Run all the tests related to the REGHDFE command
// -------------------------------------------------------------
// -Run- instead of -Do- for quick checkup

	set more off
	clear all
	set trace off
	set tracedepth 2
	set traceexpand on
	cls
	set maxvar 5000
	timer clear
	
	* adopath + "D:\Dropbox\Projects\stata\hdfe\tests"
	adopath + "D:\Dropbox\Projects\stata\hdfe\code\tests"
	adopath + "D:\Dropbox\Projects\stata\hdfe\code\_reghdfe"
	adopath + "D:\Dropbox\Projects\stata\hdfe\code\_common"
	cd "D:\Dropbox\Projects\stata\hdfe\code" // Else it won't find hdfe.mata
	mata: mata mlib index

	* Parameters
	local use0 "qui sysuse auto, clear"
	local use1 "qui use ../testdata/felsdvsimul, clear"
	local use2 "qui use ../testdata/callao_tmp, clear"
	local use3 "qui use ../testdata/big_lima, clear"
	
* Lists of test to run (b/c running them all takes a lot of time)
local TestsExcluded 
local TestSuite ConnectedGroups OLS1 OLS2 OLS3 OLS_Exhaustive DoF Interactions IV AvgE
global AllOk = 1

// -------------------------------------------------------------
// Programs
// -------------------------------------------------------------
cap pr drop test_connected
program test_connected
args id1 id2
	assert "`id1'"!=""
	assert "`id2'"!=""
	conf var `id1' `id2'
	timer on 1
	qui __makegps , id1(`id1') id2(`id2') groupid(ans1)
	timer off 1
	sort ans1
	qui count if ans1!=ans1[_n-1]
	local benchmark = r(N)
	
	* This is destructive
	*qui cd ../code
	mata: VERBOSE = 1
	timer on 2
	ConnectedGroups `id1' `id2', gen(ans3)
	timer off 2
	*qui cd ../tests
	*return list
	assert r(groups)==`benchmark'
	drop ans*
	noi di in ye  "  OK (id1=`id1', id2=`id2', #groups=`r(groups)')"
end
// -------------------------------------------------------------
cap pr drop prepare_one
program prepare_one
syntax [, CUT]
	matrix b = e(b)
	matrix V = e(V)
	scalar N = e(N)
	scalar F = e(F)
	if ("`cut'"!="") {
		matrix b = b[1, 1..2]
		matrix V = V[1..2, 1..2]
		
		noi testparm x*
		scalar F = r(F)
	}
end
// -------------------------------------------------------------
cap pr drop test_one
program test_one
syntax [, CUT]
	matrix bb = e(b)
	matrix VV = e(V)
	
	* Constant is meaningless in areg .. i.i , a(j)
	*di "<`cut'>"
	if ("`cut'"!="") {
		matrix bb = bb[1, 1..2]
		matrix VV = VV[1..2, 1..2]
	}
	*matrix list b
	*matrix list bb
	
	matrix eps = mreldif(b, bb)
	local eps = eps[1,1]
	assert `eps'<epsfloat()
	
	matrix eps = mreldif(V, VV)
	local eps = eps[1,1]
	assert `eps'<epsfloat()
	
	assert N==e(N)
	assert abs(F-e(F))<1e-4
	
	noi di in ye  "  OK (`e(cmdline)')"
end

// -------------------------------------------------------------
// Test -ConnectedGroups- (count number of mobility groups)
// -------------------------------------------------------------
if (`:list posof "ConnectedGroups" in TestSuite') {
	noi di in ye "Testing ConnectedGroups.ado"
	`use1'
	noi test_connected i j
	`use2'
	noi test_connected i j1
	`use2'
	noi test_connected i j2
	`use2'
	noi test_connected j1 j2
	`use2'
	noi test_connected j2 j1
	}
	
// -------------------------------------------------------------
// 1 HDFE, OLS
// -------------------------------------------------------------
if (`:list posof "OLS1" in TestSuite') {
	noi di in ye "Testing coefficients with one HDFE"
	
* Unadjusted errors
	noi di in ye "  Unadjusted VCE" _c
	`use1'
	areg y x*, absorb(i)
	prepare_one
	tsset i t
	xtreg y x*, fe
	noi test_one
	reghdfe y x* , a(i)
	noi test_one
	reghdfe y x* , a(i) vce(unadj)
	noi test_one
	
* Robust errors

	areg y x*, absorb(i) robust
	prepare_one
	reghdfe y x* , a(i) vce(robust)
	noi test_one

* Nested clustered errors
* (with xtreg, as areg/reg will have insufficient rank due to including FEs in K)

	xtreg y x*, fe robust
	prepare_one
	reghdfe y x* , a(i) vce(cluster i)
	noi test_one

* Non-nested clustered errors
	areg y x*, absorb(g) vce(cluster i) //  276 levels of j2, 5709 levels of j1
	prepare_one
	reghdfe y x* , a(g) vce(cluster i)
	noi test_one
	
* Larger sample
	`use2'
	areg y x*, absorb(i)
	prepare_one
	reghdfe y x* , a(i)
	noi test_one

	areg y x*, absorb(j1)
	prepare_one
	reghdfe y x* , a(j1)
	noi test_one

	areg y x*, absorb(j2)
	prepare_one
	reghdfe y x* , a(j2)
	noi test_one

* Non-nested clustered errors
	areg y x*, absorb(j2) vce(cluster j1) //  276 levels of j2, 5709 levels of j1
	prepare_one
	reghdfe y x* , a(j2) vce(cluster j1)
	noi test_one
	
}

// -------------------------------------------------------------
// 2 HDFE, OLS
// -------------------------------------------------------------
if (`:list posof "OLS2" in TestSuite') {
	noi di in ye "Testing coefficients with two HDFE"
	
	`use1'
	areg y x* ibn.j, absorb(i)
	prepare_one, cut
	reghdfe y x* , a(i j) tol(1e-12)
	test_one, cut
	reg2hdfe y x*, id1(i) id2(j) // Just for comparison
	
	areg y x* ibn.j, absorb(i) robust
	prepare_one, cut
	reghdfe y x* , a(i j) tol(1e-12) vce(rob)
	test_one, cut
	
	areg y x* ibn.g, absorb(t) cluster(i)
	prepare_one, cut
	reghdfe y x* , a(g t) tol(1e-12) vce(cluster i)
	test_one, cut
	
	// both xtreg and hdfe notice that -i- is both cluster and FE
	xtreg y x* ibn.t, fe cluster(i)
	matrix list e(V)
	prepare_one, cut
	reghdfe y x* , a(i t) tol(1e-12) vce(cluster i)
	ereturn list
	matrix list e(V)
	test_one, cut

* Add collinear regressors
	clonevar z1=x1
	clonevar z2=x2
	xtreg y x* z* ibn.t, fe cluster(i)
	prepare_one, cut
	reghdfe y x* z*, a(i t) tol(1e-12) vce(cluster i)
	test_one, cut
	drop z*
	
	`use2'
	reg2hdfe y x*, id1(i) id2(j1) // Slow
	* note that their FStat is different (we use a wald test like areg and xtreg)
	local x1 = _b[x1]
	local x2 = _b[x2]
	reghdfe y x* , a(i j1)
	assert abs(_b[x1]-`x1')<1e-4
	assert abs(_b[x2]-`x2')<1e-4
	noi di in ye  "  OK"
}

// -------------------------------------------------------------
// DoF Corner cases
// -------------------------------------------------------------
if (`:list posof "DoF" in TestSuite') {
	noi di in ye "Testing corner cases for DoF"

* FE t is nested within the cluster var
* This means no group of that FE spans more than one cluster
* Since one cluster = one superobs, we don't lose DoF for that FE since the FE estimates averages
* "within an obs"
	
	`use1'
	gen tt = int(t/2)
	*tab t tt
	areg y x*, absorb(t) cluster(tt)
	prepare_one, cut
	reghdfe y x* , a(t) tol(1e-10) vce(cluster tt)
	cap test_one, cut
	assert _rc // -areg- doesn't adjust for this
	noi di in ye  "  OK"
	
	* We can tell -hdfe- to avoid adjusting it
	reghdfe y x* , a(t) tol(1e-10) vce(cluster tt) dofmethod(naive)
	test_one, cut
	noi di in ye  "  OK"
	
* -g- nests -i-
	reghdfe y x*, a(i t g)
	assert e(dofmethod)=="bounds"
	assert e(K1)==20 & e(M1)==1
	assert e(K3)==5 & e(M3)==5 // Since FE3 is redundant given FE1
	noi di in ye  "  OK"

* Does it run with just LHS?
	reghdfe y, abs(i)
	assert e(F)==0
	noi di in ye  "  OK"
	reghdfe y, abs(i) vce(cluster i)
	assert e(F)==0 // VCV for constant will not be identified
	noi di in ye  "  OK"
	}

// -------------------------------------------------------------
// 3 HDFE, OLS
// -------------------------------------------------------------
if (`:list posof "OLS3" in TestSuite') {
	noi di in ye "Testing OLS 3 HDFE"
	`use1'
	
	areg y x* ibn.j ibn.t, absorb(i)
	prepare_one, cut
	reghdfe y x* , a(i j t) tol(1e-12)
	test_one, cut
	noi di in ye  "  OK"
}

// -------------------------------------------------------------
// OLS Interactions
// -------------------------------------------------------------
if (`:list posof "Interactions" in TestSuite') {
	noi di in ye "Testing i.foo#i.bar and i.foo#c.bar interactions"
	`use1'

* Between categories
	egen ID = group(g t)
	areg y x* , absorb(ID)
	prepare_one, cut
	reghdfe y x* , a(g#i.t)
	test_one, cut
	drop ID
	
	* Corner case where both x's are collinear with the FEs
	egen ID = group(i t)
	areg y x* , absorb(ID)
	local df_m = e(df_m)
	reghdfe y x* , a(i#i.t)
	assert `df_m'==0 & e(df_m)==0
	noi di in ye  "  OK"

* Between a category and a cont. var (i.foo#c.bar)
	gen z = uniform()
	reg y x* ibn.g#c.z
	prepare_one, cut
	gen byte c = 1
	reghdfe y x* , a(c g#c.z)
	test_one, cut
	drop z c

* Two interactions with cont vars!
	gen z1 = uniform()
	gen z2 = uniform()
	reg y x* ibn.g#c.z1 ibn.g#c.z2
	prepare_one, cut
	gen byte c = 1
	reghdfe y x* , a(c g#c.z1 g#c.z2)
	test_one, cut
	drop z1 z2 c


* Both
	gen z = uniform()
	reg y x* i.g#c.z i.t
	prepare_one, cut
	reghdfe y x* , a(t g#c.z)
	test_one, cut
	drop z
}
// -------------------------------------------------------------
// Exhaustive 1 HDFE OLS
// -------------------------------------------------------------
if (`:list posof "OLS_Exhaustive" in TestSuite') {
	noi di in ye "Testing IV"
	`use1'
	test_hdfe y x1 , abs(i)
	test_hdfe y x* , abs(g#t)
	* TODO: cluster, etc.
}

// -------------------------------------------------------------
// IV
// -------------------------------------------------------------
if (`:list posof "IV" in TestSuite') {
	noi di in ye "Testing IV"
	di as error "MANUAL: IV checks not automated"
	`use2'

	// Absorbing -i- is EXTREMELY slow on ivreg2 (999 levels) but else we won't spot the differences easily
	qui tab i, gen(FE) 
	drop FE1
	
	* PROBLEM: I WANT TO SHOW FIRST ONLY AND THEN THE REST!!!!!
	
	
	* TODO: Do automatic test
	
	* So far manually,
	* with -first- Uncentered TSS and UncenteredR2 are different
	
	* Massively slow
	*ivreg2 y x1 FE* (x2=j2) , small partial(FE*) first
	*hdfe y x1 (x2=j2) , absorb(i) tol(1e-10)

	cap drop FE*
	qui tab ubigeo, gen(FE)
	drop FE1
	gen z = uniform()

	* See if it runs with many FEs
	reghdfe y (x1 = x2) , absorb(i t ubigeo#c.z) tol(1e-10) ivsuite(ivregress) 

	* Check ivreg2
	ivreg2 y x1 FE* (x2=j2) , small partial(FE*) first
	reghdfe y x1 (x2=j2) , absorb(ubigeo) tol(1e-10)
	
	* Check ivregress
	ivregress 2sls y x1 FE* (x2=j2) , small
	test x1 x2
	reghdfe y x1 (x2=j2) , absorb(ubigeo) tol(1e-10) ivsuite(ivregress)

	* More endog
	ivregress 2sls y FE* (x1 x2= j1 j2) , small
	test x1 x2
	reghdfe y (x1 x2 = j1 j2) , absorb(ubigeo) tol(1e-10) ivsuite(ivregress)

	* Overid
	ivregress 2sls y x1 FE* (x2= j1 j2) , small
	test x1 x2
	reghdfe y x1 (x2 = j1 j2) , absorb(ubigeo) tol(1e-10) ivsuite(ivregress)
	
	* Robust
	ivregress 2sls y x1 FE* (x2= j1 j2) , small vce(robust)
	test x1 x2
	reghdfe y x1 (x2 = j1 j2) , absorb(ubigeo) tol(1e-10) ivsuite(ivregress) vce(robust)
	
	* Cluster
	ivregress 2sls y x1 FE* (x2= j1 j2) , small vce(cluster i)
	test x1 x2
	reghdfe y x1 (x2 = j1 j2) , absorb(ubigeo) tol(1e-10) ivsuite(ivregress) vce(cluster i)

}

// -------------------------------------------------------------
// AvgE
// -------------------------------------------------------------
if (`:list posof "AvgE" in TestSuite') {
	`use1'
	noi di in ye "Testing AvgE"
	set trace off
	reghdfe y x* , abs(i) avge(j)
	reghdfe y x* , abs(i j)
	reghdfe y x* , abs(i) avge(j t)

	* Not needed
	gen foo = x1
	reghdfe y x* foo, abs(i ) avge(j)
}


// -------------------------------------------------------------
// Advanced options, OLS
// -------------------------------------------------------------
* Nested F_absorb
	reghdfe y x*, abs(i t) nested
	
* Ensure that -check- works
	reghdfe y x*, abs(i t) check
	
* Save FEs
	reghdfe y x*, abs(FE1=i)
	drop FE1

exit
