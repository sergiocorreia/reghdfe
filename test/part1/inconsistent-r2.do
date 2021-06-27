noi cscript "reghdfe: R2 and other results were inconsistent with slope-only absvars" adofile reghdfe

* Prevent bug regression of https://github.com/sergiocorreia/reghdfe/issues/78

* Dataset
	sysuse auto

	reghdfe price weight gear length, a(turn trunk#c.head) //old
	storedresults save benchmark e()
	
	reghdfe price weight gear length, a(trunk#c.head turn) //old
	storedresults compare benchmark e(), tol(1e-11) exclude(macro: cmdline extended_absvars absvars matrix: dof_table)

storedresults drop benchmark
exit
