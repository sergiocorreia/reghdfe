**** PREPARE

* Reload
	cap ado uninstall reghdfe
	loc dir `c(pwd)'
	loc dir : subinstr loc dir "test" "src"
	di as text "DIR: `dir'"
	net install reghdfe, from("`dir'")


* Clean slate
	set more off
	set trace off
	set varabbrev on
	do setup // cleanup, and reinstall reghdfe and dependencies
	cap cls // uncomment

**** RUN CERTIFICATION SCRIPTS

	run unadjusted
	run unadjusted-if
	run unadjusted-singular
	run robust
	run cluster
	run cluster-equals-robust
	run cluster-nested // Prevent regression of bug reported by Michael Wittry (see email)
	run cluster-strictly-nested
	run cluster-string
	run cluster-mwc

	run missing-absvar
	run options
	run omitted
	run alphas
	
	cap run cluster-after-expand
	if (c(rc)) di as error "TEST FAILED AS EXPECTED"
	if (c(rc)) di as error "This is a variant of the problem raised by Cattaneo et al 2017"
	
	run predict
	run predict-pweight
	run predict-fweight
	run predict-aweight
	run predict-slope
	run predict-slope-only
	run predict-slope-only-nocons
	
	run noabsorb
	run more-weights
	run pweight-cluster
	run fvstrip // Test cases with complex factor variable 
	run stdata
	run singleton-x
	run inconsistent-r2
	run singleton-and-fweights

	run bug_cluster

	* -prune- corner cases
	* run prune


	* run hdfe


exit
	





* OLD OR NOT-YET-REVIEWED TESTS:
	run test-fvvarlist
	run test-robust
	run test-zero
	run test-dof
	run test-cluster-nested 
	run test-vce-complex-bw
	run test-weights
	run test-mwc
	run test-ivreg2-nocons

	run test-weights-1clust
	run test-weights-2clust // 2 absvars and 2 clustervars

	run test-cont-interaction // fixes bug when the expansion of i.var#c.var  does not include base levels

	run test-gmm
	* Doesn't work: run test-cue
	run test-liml
	// run test-singletons
	run test-iv
	run test-slope

	run test-rank
	
	run test-postestimation-predict-pweight // predict after pweight
	run test-postestimation-stdp // predict stdp
	run test-postestimation-test
	run test-attach

	run test-cache

	run test-stages
	
	run test-hdfe // just tests that the syntax works, not for correctness

	* Don't run by default this as it's mostly a test about quipu.ado and not reghdfe.ado
	*run test-quipu

* Success
	di as text "No errors found!"
	clear
	exit
