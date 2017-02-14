* ESTFE - Allow easy FE rows with estout
* See example at the end

capture program drop estfe
program define estfe
	syntax [anything(id="stored estimates" name=est_list)], [restore labels(string asis)]
	if ("`restore'"!="") Restore `est_list'
	else Add `est_list', labels(`labels')
end

capture program drop Add
program define Add, rclass
	syntax [anything(id="stored estimates" name=est_list)], [labels(string asis)]
	local dot .
	local hasdot : list dot in est_list
	local est_list : list est_list - dot

	if ("`est_list'"!="") {
		qui estimates dir `est_list'
		local models "`r(names)'"
	}

	if (`hasdot') {
		tempname hold
		estimates store `hold', nocopy
		local models `hold' `models'
	}

	foreach model of local models {
		AddOne `model' // injected `absvars'
		local fe_list `fe_list' `absvars'
	}
	local fe_list : list uniq fe_list

	if (`hasdot') {
		qui estimates restore `hold'
		estimates drop `hold'
	}

	local indicate_fe // This will contain our answer
	while (`"`labels'"'!="") {
		gettoken lhs labels : labels
		gettoken rhs labels : labels
		if "`rhs'"=="" {
			di as error "error: odd number of labels"
			error 123
		}

		foreach fe of local fe_list {
			local fixed_fe : subinstr local fe "0." "", all
			if ("`fixed_fe'"=="`lhs'") {
				local indicate_fe `"`indicate_fe' "`rhs'=`fe'" "'
				local fe_list : list fe_list - fe
				continue, break
			}
		}
	}

	* Parse remaining (w/out label)
	foreach fe of local fe_list {
		local fixed_fe : subinstr local fe "0." "", all
		local indicate_fe `"`indicate_fe' "`fixed_fe'=`fe'""'
	}

	return local indicate_fe `"`indicate_fe'"'
end

capture program drop AddOne
program define AddOne, eclass
	* From Ben Jann
	* See https://github.com/benjann/estout/issues/6
	* Requires erepost from SSC
	args model

	qui estimates restore `model'
	
	* Backup e(b) e(V)
	tempname b V new
	matrix `b' = e(b)
	matrix `V' = e(V)
	ereturn matrix b_backup = `b', copy
	ereturn matrix V_backup = `V', copy

	* Augment, reghdfe convention
	local K = colsof(`b')
	local G = e(N_hdfe_extended)
	local absvars "`e(extended_absvars)'"
	
	* Also allow areg
	if (`G'==.) local G = 1
	if ("`absvars'"=="") local absvars "`e(absvar)'"
	FixAbsvars `absvars'

	* Allow xtreg_fe, xtivreg_fe, etc
	if ("`absvars'"=="" & "`e(model)'"=="fe") local absvars "`e(ivar)'"

	matrix `new' = J(1, `G', 0)
	matrix colnames `new' = `absvars'
	matrix `b' = `b', `new'
	matrix `V' = (`V' , J(`K', `G', 0)) \ (J(`G', `K', 0), J(`G', `G', 0))

	erepost b=`b' V=`V', rename // Minor problem: removes "hidden" attribute
	estimates store `model', nocopy
	c_local absvars "`absvars'"
end

capture program drop FixAbsvars
program define FixAbsvars
	while ("`0'"!="") {
		gettoken absvar 0 : 0
		local newabsvar
		while ("`absvar'"!="") {
			gettoken part absvar : absvar, parse("# ")
			if (strpos("`part'", "#")==0 & strpos("`part'", "c.")==0) local part 0.`part'
			local newabsvar `newabsvar'`part'
		}
		local newabsvars `newabsvars' `newabsvar'
	}
	c_local absvars `newabsvars'
end

capture program drop Restore
program define Restore, eclass
	syntax [anything(id="stored estimates" name=est_list)]
	local dot .
	local hasdot : list dot in est_list
	local est_list : list est_list - dot

	if ("`est_list'"!="") {
		qui estimates dir `est_list'
		local models "`r(names)'"
	}

	if (`hasdot') {
		tempname hold
		estimates store `hold', nocopy
		local models `hold' `models'
	}

	foreach model of local models {
		qui estimates restore `model'
		tempname b V
		matrix `b' = e(b_backup)
		matrix `V' = e(V_backup)
		ereturn local b_backup
		ereturn local V_backup
		erepost b=`b' V=`V', rename
		estimates store `model', nocopy
	}

	if (`hasdot') {
		qui estimates restore `hold'
		estimates drop `hold'
	}
end

/*
* Setup
	pr drop _all
	set trace off
	clear all
	set more off
	sysuse auto

	bys turn: gen t = _n
	xtset turn t
	
* Run and store regressions
	reghdfe price weight, a(turn foreign#trunk##c.gear) keepsing
	estimates store model1, nocopy

	reghdfe price weight length, a(foreign turn) keepsing
	estimates store model2, nocopy
	
	areg price length, a(turn)
	estimates store model3, nocopy

	regress price weight
	estimates store model4, nocopy

	xtreg price gear, fe
	estimates store model5, nocopy

	xtreg price gear length, re
	*estimates store model6, nocopy

* Prepare estimates for -estout-
	estfe . model*, labels(turn "Turn FE" foreign#trunk "Foreign-Trunk FE" foreign#trunk#c.gear_ratio "Foreign-Trunk Gear Slope")
	return list
	
* Run estout/esttab
	esttab . model* , indicate("Length Controls=length" `r(indicate_fe)') varwidth(30)
		
* Return stored estimates to their previous state
	estfe . model*, restore

* Verify
	areg
	estimates dir _all
	estimates restore model2
	reghdfe
	di e(b_backup) // gives error if not restored
*/
