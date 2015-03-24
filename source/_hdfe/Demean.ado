cap pr drop Demean
program define Demean

	CheckCorrectOrder demean
	syntax , VARlist(varlist numeric) ///
		[TOLerance(real 1e-7) MAXITerations(integer 10000) ACCELerate(integer 1) /// See reghdfe.Parse
		CHECK(integer 0) SAVE_fe(integer 0) /// Runs regr of FEs
		NUM_fe(integer -1)] /// Regress only against the first Nth FEs (used in nested Fstats)
		[bad_loop_threshold(integer 1) stuck_threshold(real 5e-3) pause_length(integer 20) ///
		accel_freq(integer 3) accel_start(integer 6)] /// Advanced options

	assert inrange(`tolerance', 1e-20, 1) // However beyond 1e-16 we reach the limits of -double-
	assert inrange(`maxiterations',1,.)
	assert inlist(`accelerate',0,1)
	assert inlist(`check',0,1)
	assert inlist(`save_fe',0,1)
	assert inrange(`num_fe',1,100) | `num_fe'==-1 // -1 ==> Use all FEs

	assert `bad_loop_threshold'>0
	assert `stuck_threshold'>0 & `stuck_threshold'<=1
	assert `pause_length'>=0
	assert `accel_freq'>=0
	assert `accel_start'>0

	* We need to recast everything to -double- (-float- is not good enough)
	Debug, level(2) msg("(recasting variables as -double-)")
	recast double `varlist'

	* We can't save the FEs if there is more than one variable
	cap unab _ : `varlist', max(1)
	Assert (_rc==0 | `save_fe'==0) , rc(`=_rc') ///
		msg("hdfe.Demean: cannot save FEs of more than one variable at a time")

	tempvar resid
	local save = `save_fe' | `check' // check=1 implies save_fe=1
	local base_args `""`resid'", `tolerance', `maxiterations', `save', `accelerate', `num_fe'"'
	local adv_args `"`bad_loop_threshold', `stuck_threshold', `pause_length', `accel_freq', `accel_start'"'
	local args `"`base_args', `adv_args'"'
	Debug, level(3) msg(" - Structure of Mata calls: make_residual(" as result "{variable}" as text `", `args')"')

	Debug, level(2) tic(30)
	mata: st_local("weightexp", weightexp)
	
	foreach var of varlist `varlist' {
		cap drop __Z*__
		Assert !missing(`var'), msg("hdfe.Demean error: `var' has missing values and cannot be transformed")
		
		* Syntax: MAKE_RESIDUAL(var, newvar, tol, maxiter | , save=0 , accel=1, first_n=`num_fe')
		* Note: summarize doesn't allow pweight ( see http://www.stata.com/support/faqs/statistics/weights-and-summary-statistics/ )
		* Since we only want to compute means, replace with [aw]
		local tmpweightexp = subinstr("`weightexp'", "[pweight=", "[aweight=", 1)
		qui su `var' `tmpweightexp', mean
		char define `var'[mean] `r(mean)'
		mata: make_residual("`var'", `args')
		assert !missing(`resid')

		* Check that coefs are approximately 1
		if (`check') {
			unab _ : __Z*__, min(1)
			local backup = ("`e(cmd)'"!="")
			if (`backup') {
				tempname backup_results
				est store `backup_results', nocopy // nocopy needed to avoid having e(_estimates_name)
			}
			qui _regress `var' __Z*__
			local label : var label `var'
			if ("`label'"=="") local label `var'
			di as text "FE coefficients for `label':{col 36}" _continue
			foreach z of varlist __Z*__ {
				assert !missing(`z')
				di as text " `=string(_b[`z'], "%9.7f")'"  _continue
			}
			di
			
			if (`backup') qui est restore `backup_results'
			if (!`save_fe') cap drop __Z*__
		}

		* If the tol() is not high enough (e.g. 1e-14), we may fail to detect variables collinear with the absorbed categories
		* Again, we can't use pweight with summarize, but in this case it's just for debugging purposes so use [aw]
		qui su `resid' `tmpweightexp'
		local prettyvar `var'
		if (substr("`var'", 1, 2)=="__") local prettyvar : var label `var'
		if inrange(r(sd), 1e-20 , epsfloat()) di in ye "(warning: variable `prettyvar' is probably collinear, maybe try a tighter tolerance)"

		qui replace `var' = `resid' // This way I keep labels and so on
		drop `resid'
		Assert !missing(`var'), msg("REGHDFE.Demean: `var' has missing values after transformation")
	}
	Debug, level(2) toc(30) msg("(timer for calls to mata:make_residual)")
end

