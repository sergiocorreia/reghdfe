* Config
	log close _all
	clear all
	discard
	cap cls
	set more off
	cd "D:/Github/reghdfe/test"
	
* Run scripts
	*run test-minimal
	set trace off
	set tracedepth 5
	set varabbrev on

	run test-unadjusted
	run test-robust
	run test-cluster
	run test-ivreg2
	run test-cluster-same-as-absvar
	run test-vce-complex-bw
	run test-weights
	run test-mwc
	run test-ivreg2-nocons

	run test-cont-interaction

	run test-gmm
	* Doesn't work: run test-cue
	run test-liml

	run test-singletons
	run test-iv
	run test-slope
	
	run test-avar
	
	run test-postestimation-test
	run test-postestimation-predict

	run test-cores
	run test-attach
	run test-stages

	run test-cache
	run test-over

	run test-hdfe // just tests that the syntax works, not for correctness

	* Don't run by default this as it's mostly a test about quipu.ado and not reghdfe.ado
	*run test-quipu

* Success
	di as text "No errors found!"
	clear
	exit
