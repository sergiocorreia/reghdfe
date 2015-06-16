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
	* reghdfe // replay
	
	di
	reghdfe, version

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


local fn 2e6 // 1e5 2e6 1e7 1e8
local reps 2

* BUGBUG: Borrow a server with 8+ cores
forv n=1/8 {
	set processors `n'
	UseFile `fn'
	RunMany `reps'
	local avg_time = `r(time)' / `reps'
	di as text `"[`ver']=[`avg_time']"'
}

exit

Avg time per core:
1 13.7s
2 11.6	(85.16%)
3 11.2	(81.75%)
4 10.7	(77.92%)

Explanation:
- Precompute step has parallelizes poorly or not at all? (stays around 3.5s)
- map_solve() parallelizes poorly (7.1s to 5.7s)
