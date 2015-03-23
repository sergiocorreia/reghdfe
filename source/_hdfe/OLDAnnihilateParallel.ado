
// The filename must have already been saved
cap pr drop AnnihilateParallelOld
program AnnihilateParallelOld
syntax, VARlist(varlist numeric) FILEname(string) UID(varname numeric) Numcores(integer) [*]



	* Save Mata structures
	local objects VERBOSE G FEs fe_cache prev_numstep

	tempfile hdfe_mata
	// local hdfe_mata = "borrar.mo"
	qui mata: mata matsave "`hdfe_mata'" `objects' , replace

	* Create temporary do-files
	forv i=1/`numcores' {
		
		* tempfile dofile`i'
		local code = string(int( (10*`i'+uniform() ) * 1e8 ), "%14.0f") // CORE + 8 RANDOM DIGITS
		local name`i' __reghdfe__`code'
		local dofile`i' `c(tmpdir)'`name`i''.tmp

		tempfile dataset`i'
		
		qui file open fh using "`dofile`i''", write text replace
		file write fh `" set more off "' _n
		file write fh `" sleep 10000 "' _n
		file write fh `" use `uid' `varlist`i'' using "`filename'" "' _n
		file write fh `" rename `uid' __uid__ "' _n

		file write fh `" cap qui do reghdfe_absorb.ado"' _n // Else the debug build will fail due to -clear mata-
		file write fh `" mata: mata matuse "`hdfe_mata'" "' _n

		file write fh `" reghdfe_absorb, step(demean) varlist(`varlist`i'') `options' "' _n

		file write fh `" rename __uid__ `uid' "' _n
		file write fh `" keep `uid' `varlist`i'' "' _n
		file write fh `" save "`dataset`i''" "' _n
		file write fh `" exit, clear "' _n _n
		file close fh
		// type "`dofile`i''"

		local cores `cores' `i' // Will contain remaining cores
	}

	* Misc
	drop `varlist' // basically, keeps UID and clustervar
	GetBinary
	local binary `r(binary)'
	cap shell erase __reghdfe__*.log

	* Call subprocesses
	forv i=1/`numcores' {
		sleep 1000 // increase this
		Debug, level(1)  msg(`" - Shell: winexec "`binary'" -e do "`dofile`i''" "')
		Debug, level(1)  msg(`" - (variables: `varlist`i'') "')
		winexec "`binary'" /e do "`dofile`i''"
	}

	* Wait until all instances have started
	forv i=1/`numcores' {
		local ok 0
		while !`ok' {
			sleep 100
			* If log file exists, then the Stata instance is already running
			local logfile `name`core''.log
			cap conf file "`logfile'"
			if (`=_rc'==0) local ok 1
			Debug, level(1) msg(" - process `i' started")
		}
	}

	* Wait for termination and merge
	while ("`cores'"!="") {
		foreach core of local cores {
			Debug, level(1) msg(" - sleeping")
			sleep 500 // increase this
			local logfile `name`core''.log
			local dofile  `dofile`core''
			local dtafile `dataset`core''

			* (This part assumes that the log files exist and stata is running)
			* Try to delete the do-file
			* If we are able to do so, Stata has finished
			* If we then find that the dta exists, then there was no error

			Debug, level(1) msg(" - trying to erase dofile `i'")
			cap erase "`dofile`core''"
			if (`=_rc'==0) {
				Debug, level(1) msg(" - dofile `i' succesfully erased")
				cap conf file "`dtafile'"
				if (`=_rc'!=0) {
					type "`logfile'"
					di as error "<`dtafile'> not found"
					Assert 0, msg("Call to subprocess `core' failed, see logfile above")
				}

				Debug, level(1) msg(" - Subprocess `core' done")
				local cores : list cores - core
				mata: st_local("VERBOSE",strofreal(VERBOSE))
				
				if (`VERBOSE'>=3) {
					type `name`core''.log
				}

				* Merge file
				Debug, level(1) msg(" - Merging dta #`core'")
				merge 1:1 `uid' using "`dtafile'", nogen nolabel nonotes noreport sorted assert(match)
				erase "`dtafile'"

			} // if could delete dofile
		} // for core
	} // while

	erase "`hdfe_mata'"
	cap shell erase __reghdfe__*.log

end

//------------------------------------------------------------------------------
// GetBinary
//------------------------------------------------------------------------------
// This is based on George Vega's -parallel- module
// See https://github.com/gvegayon/parallel
// Only works for windows atm
cap pr drop GetBinary
program GetBinary, rclass

	if ("${stata_executable}"!="") {

		cap conf file "${stata_executable}"
		Assert `=_rc'==0, msg(`"Stata executable not found (expected "${stata_executable}")"' _n `"please store the correct location in the global "stata_executable" "')

		return local binary "${stata_executable}"
	}

	local flv `c(flavor)'
	if ("`flv'"=="Small") local flv SM
	if ("`flv'"=="IC") local flv
	if (`c(MP)') local flv MP

	local ext ".exe"

	if (c(osdtl)!="" | c(bit)==64) local bit "-64" 
	local binary "`c(sysdir_stata)'Stata`flv'`bit'`ext'"
	
	cap conf file "`binary'"
	Assert `=_rc'==0, msg(`"Stata executable not found (expected "`binary'")"' _n `"please store the correct location in the global "stata_executable" "')

	return local binary "`binary'"
end
