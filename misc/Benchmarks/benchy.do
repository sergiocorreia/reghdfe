* Track time through versions of the code

capture program drop RunOnce
program define RunOnce
	*cap qui reghdfe v3 v2 id4 id5 id6, absorb(g1 g2) vce(cluster g1) dofmethod(naive) tol(1e-6) fast
	*if (_rc==198) qui reghdfe v3 v2 id4 id5 id6, absorb(g1 g2) vce(cluster g1) dof(pair cont) tol(1e-6) fast
	qui reghdfe v3 v2 id4 id5 id6, absorb(g1 g2) vce(cluster g1) tol(1e-6) fast
end

capture program drop UseVersion
program define UseVersion
	args ver
	if "$lastversion"=="`ver'" {
		exit
	}

	if ("`ver'"!="working") {
		global lastversion `ver'
		cap ado uninstall reghdfe
		net from https://raw.githubusercontent.com/sergiocorreia/reghdfe/`ver'/package/
		net install reghdfe
	}
	else {
		rebuild_git reghdfe
	}
	discard
end

capture program drop UseFile
program define UseFile
	args fn
	local path "D:\Dropbox\Projects\stata\hdfe\testdata" // "D:\tmp"
	use "`path'/`fn'"
end

capture program drop RunMany
program define RunMany
	args rep
	RunOnce // Warm up (load ADOs into memory, etc)
	tic
	forval i = 1/`rep' {
		di as text "." _c
		RunOnce
	}
	qui toc, report
end

// -------------------------------------------------------------------------------------------------
set more off
cls
clear
// -------------------------------------------------------------------------------------------------


local fn 1e5 // 1e5 2e6 1e7 1e8
local reps 2
local versions 	working ///
				master /// current
				/// 67a346699f698c7e77e44c100eca775b47fea321 /// first on github OK
				b00d9b5e8425d5795f8f4183fed6981d440922c6 /// mar3 OK
				4e950be9d9a0f65749e8dfd9544cd51a7a7b0765 /// mar3 CAGADO
				/// 866f85551b77fe7fda2af0aafccbbf87f8a01987 /// last with _cons
				//

foreach ver of local versions {
	qui UseVersion `ver'
	UseFile `fn'
	RunMany `reps'
	local avg_time = `r(time)' / `reps'
	di as text `"[`ver']=[`avg_time']"'
}

exit

NOTES:
The new EstimateDoF created an intrinsic slowdown because:
we now check if there are redundant coefs in FE2 wrt FE1 *EVEN IF* FE1 is a cluster variable
