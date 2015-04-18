* Track time through versions of the code

capture program drop RunOnce
program define RunOnce
	qui reghdfe v3 v2 id4 id5 id6, absorb(g1 g2) vce(cluster g1) tol(1e-6) fast
end

capture program drop UseVersion
program define UseVersion
	args ver
	if "$lastversion"=="`ver'" {
		exit
	}
	global lastversion `ver'
	cap ado uninstall reghdfe
	net from https://raw.githubusercontent.com/sergiocorreia/reghdfe/`ver'/package/
	net install reghdfe
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
	di
end

local fn 2e6 // 1e5 2e6 1e7 1e8
local reps 1
local versions 	67a346699f698c7e77e44c100eca775b47fea321 /// first on github
				866f85551b77fe7fda2af0aafccbbf87f8a01987 /// last with _cons
				master // current

foreach ver of local versions {
	qui UseVersion `ver'
	UseFile `fn'
	RunMany `reps'
	local avg_time = `r(time)' / `reps'
	di as text `"[`ver']=[`avg_time']"'
}

exit
