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
	
	run test-postestimation-test
	run test-postestimation-predict

	run test-cores
	run test-attach
	run test-stage

* Success
	di as text "No errors found!"
