noi cscript "reghdfe: prevent bug regression with cluster VCE" adofile reghdfe

* SEE:
* https://www.statalist.org/forums/forum/general-stata-discussion/general/1409974-reghdfe-4-x-standard-errors-and-categorical-variables

* Dataset
	sysuse auto, clear
	egen int turn_foreign = group(turn foreign)

* [TEST] Cluster

	* 1. Run benchmark
	reghdfe price weight, a(turn#trunk foreign) vce(cluster turn_foreign)
	storedresults save benchmark e()
	
	* 2. Run reghdfe with #
	reghdfe price weight, a(turn#trunk foreign) vce(cluster turn#foreign)
	storedresults compare benchmark e(), tol(1e-12) exclude(macro: clustvar1 clustvar cmdline)
	storedresults drop benchmark

exit
