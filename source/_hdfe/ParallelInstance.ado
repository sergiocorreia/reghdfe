cap pr drop ParallelInstance
program ParallelInstance
	syntax, core(integer) code(string asis)
	set more off
	assert inrange(`core',1,32)
	local path "`c(tmpdir)'reghdfe_`code'`c(dirsep)'"
	cd "`path'"
	set processors 1

	file open fh using "`core'_started.txt" , write text all
	file close _all

	cap noi {
		set linesize 120
		log using `core'_log.log, text

		mata: mata matuse "hdfe_mata.mo"
		mata: st_local("cores",strofreal(parallel_cores))
		assert `core' <= `cores'
		mata: st_local("usedta",parallel_dta)
		mata: st_local("vars",parallel_vars[`core'])
		mata: st_local("weightvar",weightvar)
		mata: st_local("opt",parallel_opt)
		Debug, msg(" - This is core `core'/`cores'")
		sleep 100
	
		local outfn "`core'_output.dta"
		conf new file "`outfn'"

		use `vars' `weightvar' using "`usedta'"
		de, full
		Demean, varlist(`vars') `opt'
		keep `vars'
		save `"`outfn'"'
		log close _all
	}

	local rc = _rc
	sleep 100

	if `rc'>0 {
		di in red "ERROR: `rc'"
		file open fh using "`core'_error.txt" , write text all
		file close _all
	}
	else {
		file open fh using "`core'_ok.txt" , write text all
		file close _all
	}

	file open fh using "`core'_done.txt" , write text all
	file close _all
	exit, STATA
end

