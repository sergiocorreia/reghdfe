* Config
	log close _all
	clear all
	discard
	cap cls
	
* Run scripts
	run test-weights
	* Asegurarme que falle con fw y noninteger!
* Success
	di as text "No errors found!"
