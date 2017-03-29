	log close _all
	clear all
	discard
	cap cls
	cap pr drop _all

* Dependencies
/*
	cap ado uninstall ftools
	net install ftools, from("C:/git/ftools/src")
	ftools, check
	ftools, compile
*/

* Uninstall reghdfe
	cap ado uninstall reghdfe

* Check that the .mata files don't have compile-time errors
* (inspect by hand for the word "unused")
	cd "../src"
	do reghdfe.mata
	cd "../test"

* Reininstall and compile reghdfe
* Note: this has to be hardcoded!
	net install reghdfe , from("c:/git/reghdfe/src/")
	reghdfe, compile
	noi di as text "reghdfe has been reinstalled from local, and recompiled"

exit
