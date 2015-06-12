* Config
	log close _all
	clear all
	discard
	cap cls
	set more off
	cd "D:/Github/reghdfe/test"

	rebuild_git reghdfe
	cap pr drop _all
	discard
	
* Run scripts
	adopath + "D:\Github\reghdfe\source\internal"
	adopath + "D:\Github\reghdfe\source\common"
	run test-absvars // Need to fix so we can find the path of the ADO
	adopath - "D:\Github\reghdfe\source\internal"
	adopath - "D:\Github\reghdfe\source\common"

	*run test-minimal
	set trace off
	set tracedepth 5
	set varabbrev on

	run test-unadjusted
	run test-robust
	run test-cluster
	run test-ivreg2

	run test-zero

	run test-avar

	run test-cluster-same-as-absvar
	run test-cluster-nested // Prevent regression of bug reported by Michael Wittry (see email)
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
	
	run test-postestimation-predict
	run test-postestimation-test
	run test-attach

	run test-cache

	* TODO NOW:
	run test-stages
	
	* TODO LATER:
	// run test-hdfe // just tests that the syntax works, not for correctness

	* Don't run by default this as it's mostly a test about quipu.ado and not reghdfe.ado
	*run test-quipu

* Success
	di as text "No errors found!"
	clear
	exit
