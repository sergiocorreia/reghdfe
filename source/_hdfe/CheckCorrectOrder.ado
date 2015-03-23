capture program drop CheckCorrectOrder
program define CheckCorrectOrder
	args step

	local numstep = ("`step'"=="start") + 2*("`step'"=="precompute") + ///
		3*("`step'"=="demean") + 4*("`step'"=="save")
	Assert (`numstep'>0), msg("hdfe: -`step'- is an invalid step")

	cap mata: st_local("prev_numstep", strofreal(prev_numstep))
	if (_rc) local prev_numstep 0

	Assert (`numstep'==`prev_numstep'+1) | (`numstep'==3 & `prev_numstep'==3) ///
		, msg("hdfe: expected step `=`prev_numstep'+1' instead of step 	`numstep'")
	mata: prev_numstep = `numstep'
	Debug, msg(_n as text "{title:Running -hdfe- step `numstep'/5 (`step')}") level(3)
end
