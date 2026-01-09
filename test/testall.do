* ===========================================================================
* Run all certification tests
* ===========================================================================
	
	global tests_run 0


// --------------------------------------------------------------------------
// Set up Stata
// --------------------------------------------------------------------------

	which ftools // required


// --------------------------------------------------------------------------
// Reinstall reghdfe
// --------------------------------------------------------------------------

	do setup

// --------------------------------------------------------------------------
// Run certification scripts post v6
// --------------------------------------------------------------------------
	
	cd "./part2"

	* Quick tests:

* Disabled groupreg...
*run demo_groups
*run demo_teams
*run demo_indiv_preconditioner

	run test-aggregation
	run test-trivial-group
	run test-simple-absorb
	run test-error-no-obs
	run test-error-sample
	run test-error-missing-values
	run test-bug-sort // bug found by nconstantine when the data was not sorted by group_id

	* Test that we can save+load the HDFE object
	run test-parallel-save-load

	* Good quality tests:
	run test-alphas
	cd ".."


// --------------------------------------------------------------------------
// Run all certification scripts (pre v6)
// --------------------------------------------------------------------------
	
	cd "./part1"

	run report_constant // reghdfe,constant VS reghdfe,noconstant
	run alphas
	
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
	
	cap run cluster-after-expand
	if (c(rc)) di as error "TEST FAILED AS EXPECTED"
	if (c(rc)) di as error "This is a variant of the problem raised by Cattaneo et al 2017"
	
	run noabsorb
	run more-weights
	run fvstrip // Test cases with complex factor variable 
	
	run stdata1
	run stdata2
	run stdata3

	run singleton-x
	run inconsistent-r2
	run singleton-and-fweights
	run prune
	run memory // test -compact- and -pool(#)- options

	* Postestimation
	run predict
	run predict-pweight
	run predict-fweight
	run predict-aweight
	run predict-slope
	run predict-slope-only
	run predict-slope-only-nocons
	run pweight-cluster // also tests predict
	run estat-ic
	
	* Prevent specific bugs
	run bug_cluster
	run extreme_values // test numerical accuracy
	run bug_compact

	* Extra
	do margins // not yet cscript
	run groupvar // mobility group; not yet cscript

	* do lsmr // TODO


	cd ".."


set linesize 120
di as text "TESTS COMPLETED SUCCESSFULLY!"
exit





/*
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
*/

* Success
	di as text "No errors found!"
	clear
	exit
