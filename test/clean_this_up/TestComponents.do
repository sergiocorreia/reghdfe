//------------------------------------------------------------------------------
// Run all the tests related to the REGHDFE command
//------------------------------------------------------------------------------
// -Run- instead of -Do- for quick checkup

	set more off
	clear all
	set trace off
	set tracedepth 3
	set traceexpand on
	cap cls
	set maxvar 5000
	timer clear
	discard // pr drop _all
	
	cd "D:\Dropbox\Projects\stata\hdfe\code"
	*cd build
	*qui do build_mata
	*cd ..

//------------------------------------------------------------------------------
// _reghdfe_absorb/GenerateID.ado
//------------------------------------------------------------------------------
	adopath + "_reghdfe_absorb/"
	sysuse auto
	drop if missing(rep)
	egen byte alt1 = group(mpg)
	egen int alt2 = group(rep foreign)
	
	GenerateID mpg, replace
	rename mpg id1
	
	GenerateID rep foreign, gen(id2)
	
	assert id1==alt1
	assert id2==alt2
	
	noi di as result "[OK] " as text "_reghdfe_absorb/GenerateID.ado" _n

//------------------------------------------------------------------------------
// _reghdfe_absorb/AverageOthers.ado
//------------------------------------------------------------------------------
	* TODO

//------------------------------------------------------------------------------
// reghdfe_absorb
	// adopath + "."
	adopath - "_reghdfe_absorb/"
//------------------------------------------------------------------------------			
		clear all
		mata: VERBOSE=0
		sysuse auto , clear
		local vars price weight
	reghdfe_absorb, step(start) absorb(foreign turn#mpg) avge(rep#foreign) clustervar1(foreign#turn)
	drop if missing(rep)
	reghdfe_absorb, step(precompute) keepvars(`vars') depvar(price) excludeself
		local clustervar1 = r(clustervar1)
		di as input "<`clustervar1'>"
	reghdfe_absorb, step(demean) varlist(`vars') maxiter(5e3) tol(1e-9)
		de, all
		char list
		reg `vars'
	reghdfe_absorb, step(stop)
//------------------------------------------------------------------------------

		clear all
		local vars price weight length
		local fe turn
		sysuse auto, clear
	reghdfe_absorb, step(start) absorb(`fe')
		drop if missing(rep)
		preserve
		areg `vars', a(`fe') // benchmark
		mata: alt = st_matrix("e(b)")
	reghdfe_absorb, step(precompute) keepvars(`vars')
	reghdfe_absorb, step(demean) varlist(`vars') check(1)
		reg `vars'
		mata
			b = st_matrix("e(b)")
			(b \ alt)
			assert (sum(reldif(b, alt)) < 2e-12)
		end
		restore // back to untransformed
		predict double resid, resid
	reghdfe_absorb, step(demean) varlist(resid) check(1) save(1)
		reg `vars' __Z*__
		local cons = _b[_cons]
		replace __Z1__ = __Z1__ + `cons'
		reg `vars' ibn.`fe' , nocons
		mata: alt = st_matrix("e(b)")
		mata: alt = alt[| 3 \ length(alt) |]'
		collapse (mean) __Z1__, by(`fe') fast
		mata
			b = st_data(., "__Z1__")
			(b , alt)
			assert (sum(reldif(b, alt)) < 2e-12)
		end
	reghdfe_absorb, step(stop)
		
//------------------------------------------------------------------------------		
		local vars price weight length
		local fe rep#foreign mpg
		sysuse auto , clear
		areg `vars' rep#foreign, a(mpg)	
		local alt1 = _b[weight]
		local alt2 = _se[weight]
	reghdfe_absorb, step(start) absorb(`fe')
		drop if missing(rep)
	reghdfe_absorb, step(precompute) keepvars(`vars')
	reghdfe_absorb, step(demean) varlist(`vars')
		reg `vars'
		local val1 = _b[weight]
		local val2 = _se[weight]
	reghdfe_absorb, step(stop)
		mata: mata desc
	
		di "`val1' `alt1'"
		di "`val2' `alt2'" // This differs due to not making the DF adj
		assert abs(`val1' - `alt1') < 1e-6
		
//------------------------------------------------------------------------------			
		clear all
		mata: VERBOSE=0
		mata: mata desc
		local vars price weight length
		local fe rep#foreign#turn mpg#c.weight 
		sysuse auto , clear
	reghdfe_absorb, step(start) absorb(`fe')
		drop if missing(rep)
	reghdfe_absorb, step(precompute) keepvars(`vars')
	reghdfe_absorb, step(demean) varlist(`vars') maxiter(5e3) tol(1e-9)
		reg `vars'
	reghdfe_absorb, step(stop)
	
	* Collinear of course
	gen foo = _b[weight] * weight + _cons
	su foo
	assert r(sd)<0.001 


exit
	
//------------------------------------------------------------------------------
// Programs
//------------------------------------------------------------------------------
/*
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
*/

//------------------------------------------------------------------------------
// Test -ConnectedGroups- (count number of mobility groups)
//------------------------------------------------------------------------------
/*
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
*/

//------------------------------------------------------------------------------
// ExpandFactorVariables
//------------------------------------------------------------------------------
/*
	mata: VERBOSE = 3
	sysuse auto, clear
	ExpandFactorVariables price i.foreign
	
	gen const1 = 1
	gen const2 = 1
	ExpandFactorVariables i.const*
	return list
	
	gen t = _n
	tsset t
	ExpandFactorVariables L.price weight F(-1/1).turn
	return list
	
	sysuse auto, clear
	gen t = _n
	tsset t
	ExpandFactorVariables price i.turn i.trunk#foreign ibn.rep
	return list
	
	local vars i.mpg L(-1/1).(weight pric) foreign
	ExpandFactorVariables `vars', cache
	mata: mata desc
	mata: varlist_cache
	mata: asarray_keys(varlist_cache)
	fvunab mylist : `vars'
	foreach var in `mylist' {
		mata: st_local("ans",asarray(varlist_cache, "`var'"))
		di "<`var'>=<`ans'>"
	}
	
	* asarray(A,key,a) -- A[key] = a
	* asarray(A,key) -- A[key]
	* asarray_contains(A, key) -- key in A
	* asarray_elements(A) -- len(A)
	* asarray_keys(A) -- A.keys()
	* asarray_first(A) , asarray_next -> para loopear

	ExpandFactorVariables L1000.price
	
	darle tb
	L.i.rep
	ibn.rep
	ExpandFactorVariables i.rep if rep>=3
*/


	
	* TODO

//------------------------------------------------------------------------------
// EstimateDoF
//------------------------------------------------------------------------------
	
	sysuse auto, clear
	mata: VERBOSE = 3
	clonevar foo = turn
	reghdfe_absorb, step(start) absorb(turn foo#foreign)
	drop if missing(rep)
	
	local clustervar1 "foreign"
	set trace on
	set tracedepth 4
	reghdfe_absorb, step(precompute) keepvars(price weight) clustervar1(`clustervar1')
	local original_clustervar1 "`clustervar1'"
	local clustervar1 = r(clustervar1)
	EstimateDoF, dof_method(bounds) clustervar1(`clustervar1')
	
	asd
	
	EstimateDoF, dof_method(bounds) clustervar1(foreign)
	EstimateDoF, dof_method(bounds)
	
exit
