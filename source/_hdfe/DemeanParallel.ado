cap pr drop DemeanParallel
program define DemeanParallel
	* Notes:
	* First cluster is taking by this stata instance, to save HDD/memory/merge time
	* Also, this cluster should have more obs than the other ones so we let it have
	* the default number of processes
	* (the other start with 1 proc allowed, which should be fine)
	* Thus it will usually finish faster, to start waiting for the 2nd fastest  to merge

	CheckCorrectOrder demean
	syntax, VARlist(varlist numeric) FILEname(string) UID(varname numeric) CORES(integer) SELF(string) [*]

	local varlist : list uniq varlist
	local K : list sizeof varlist
	local cores = min(`cores',`K')
	local size = c(N) * c(width) / 2^30
	local wait = int(100 + 1000 * `size') // each gb wait 1 sec

	assert inlist("`self'", "reghdfe", "hdfe") // We will call `self', instance ...

	* Deal each variable like cards in Poker
	local core 1
	foreach var of local varlist {
		local varlist`core' `varlist`core'' `var'
		local ++core
		if (`core'>`cores') local core 1
	}

	* Folder name.. need some entropy.. use varlist + time
	mata: st_local("hash", strofreal(hash1("`varlist'"), "%20.0f"))
	local seed = real(subinstr(c(current_time),":","",.)) + `hash'
	local seed = mod(`seed',2^30) // Needs to be < 2^31-1
	set seed `seed'
	local code = string(int( uniform() * 1e6 ), "%08.0f")

	* Prepare
	* Note: On windows, tmpdir has a trailing / but on Linux it doesn't!
	local path "`c(tmpdir)'`c(dirsep)'hdfe_`code'`c(dirsep)'"
	* NOTE: Copy any changes to this line into ParallelInstance.ado

	Debug, level(1) msg(" - tempdir will be " as input "`path'")
	mata: parallel_cores = `cores'
	mata: parallel_dta = `"`filename'"'
	mata: parallel_vars = J(`cores',1,"")
	mata: parallel_opt = `"`options'"'
	mata: parallel_path = `"`path'"'
	forv i=1/`cores' {
		mata: parallel_vars[`i'] = "`varlist`i''"
	}

	local dropvarlist : list varlist - varlist1
	drop `dropvarlist' // basically, keeps UID and clustervar
	mata: st_global("hdfe_pwd",pwd())
	mkdir "`path'"
	qui cd "`path'"

	local objects VERBOSE G FEs betas prev_numstep parallel_* weightexp weightvar
	qui mata: mata matsave "`path'hdfe_mata.mo" `objects' , replace

	* Call -parallel-
	Debug, level(1) msg(" - running parallel instances")
	*qui mata: parallel_setstatadir("")
	qui parallel setclusters 1 // I just want the global PLL_DIR
	local binary `"$PLL_DIR"'
	Assert `"`binary'"'!="", msg("`self' error: after parallel, global PLL_DIR was empty")
	global PLL_DIR
	global PLL_CLUSTERS
	global PLL_STATA_PATH

	cap mata: st_local("VERBOSE",strofreal(VERBOSE))
	if (`VERBOSE'==0) local qui qui
	`qui' di as text _n "{dup 44:_}/ PARALLEL \{dup 44:_}"

	local flag = cond("`c(os)'"=="Windows", "/q", "-q")

	* Create instances
	forv i=2/`cores' {
		local cmd `"winexec `binary' `flag'  `self', instance core(`i') code(`code') "'
		Debug, level(1) msg(" - Executing " in ye `"`cmd' "')
		`cmd'
		Debug, level(1) msg(" - Sleeping `wait'ms")
		if (`i'!=`cores') sleep `wait'
	}
	Demean, varlist(`varlist1') `options' // core=1

	* Wait until all instances have started
	local timeout 20
	local elapsed 0
	forv i=2/`cores' {
		local ok 0
		while !`ok' {
			sleep 100
			local fn "`path'`i'_started.txt"
			cap conf file "`fn'"
			local rc = _rc
			if (`rc'==0) {
				local ok 1
				Debug, level(1) msg(" - process `i' started")
				erase "`fn'"
			}
			else {
				local elapsed = `elapsed' + 0.1
				Assert `elapsed'<`timeout', msg("Failed to start subprocess `i'")
			}
		}
		local remaining_cores `remaining_cores' `i' // Will contain remaining cores
	}

	* Wait for termination and merge
	while ("`remaining_cores'"!="") {
		foreach core of local remaining_cores {
			local donefn "`path'`core'_done.txt"
			local okfn "`path'`core'_ok.txt"
			local errorfn "`path'`core'_error.txt"
			local dtafn "`path'`core'_output.dta"
			local logfile "`path'`core'_log.log"


			cap conf file "`donefn'"
			local rc = _rc

			if (`rc'==0) {
				Debug, level(1) msg(" - process `core' finished")
				erase "`donefn'"
				cap conf file "`okfn'"
				if (`=_rc'>0) {
					type "`logfile'"
					//di as error "<`dtafn'> not found"
					Assert 0, msg("Call to subprocess `core' failed, see logfile")
				}

				erase "`okfn'"
				Debug, level(1) msg(" - Subprocess `core' done")
				local remaining_cores : list remaining_cores - core
				mata: st_local("VERBOSE",strofreal(VERBOSE))
				
				if (`VERBOSE'>=3) {
					type "`logfile'"
				}
				erase "`logfile'"

				* Merge file
				Debug, level(1) msg(" - Merging dta #`core'")
				merge 1:1 _n using "`dtafn'", nogen nolabel nonotes noreport sorted assert(match)
				erase "`dtafn'"
			}
			else {
				sleep 500 // increase this
			}
		}
	}

	* Cleanup
	qui cd "${hdfe_pwd}"
	erase "`path'hdfe_mata.mo"
	cap rmdir `"`path'"'
	`qui' di as text _n "{dup 44:_}\ PARALLEL /{dup 44:_}"

end
