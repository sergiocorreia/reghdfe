// -------------------------------------------------------------
// Simple debugging
// -------------------------------------------------------------
cap pr drop Debug
program define Debug

	syntax, [MSG(string asis) Level(integer 1) NEWline COLOR(string)] [tic(integer 0) toc(integer 0)]
	
	cap mata: st_local("VERBOSE",strofreal(VERBOSE)) // Ugly hack to avoid using a global
	if ("`VERBOSE'"=="") {
		di as result "Mata scalar -VERBOSE- not found, setting VERBOSE=3"
		local VERBOSE 3
		mata: VERBOSE = `VERBOSE'
	}


	assert "`VERBOSE'"!=""
	assert inrange(`level',0, 4)
	assert (`tic'>0) + (`toc'>0)<=1

	if ("`color'"=="") local color text
	assert inlist("`color'", "text", "res", "result", "error", "input")

	if (`VERBOSE'>=`level') {

		if (`tic'>0) {
			timer clear `tic'
			timer on `tic'
		}
		if (`toc'>0) {
			timer off `toc'
			qui timer list `toc'
			local time = r(t`toc')
			if (`time'<10) local time = string(`time'*1000, "%tcss.ss!s")
			else if (`time'<60) local time = string(`time'*1000, "%tcss!s")
			else if (`time'<3600) local time = string(`time'*1000, "%tc+mm!m! SS!s")
			else if (`time'<24*3600) local time = string(`time'*1000, "%tc+hH!h! mm!m! SS!s")
			timer clear `toc'
			local time `" as result " `time'""'
		}

		if (`"`msg'"'!="") di as `color' `msg'`time'
		if ("`newline'"!="") di
	}
end
