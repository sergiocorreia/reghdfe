capture program drop Wrapper_mwc
program define Wrapper_mwc
* This will compute an ols regression with 2+ clusters
syntax , depvar(varname) [indepvars(varlist) avgevars(varlist)] ///
	original_absvars(string) original_depvar(string) [original_indepvars(string) avge_targets(string)] ///
	vceoption(string asis) vcetype(string) ///
	kk(integer) ///
	[weightexp(string)] ///
	addconstant(integer) ///
	[SUBOPTions(string)] [*] // [*] are ignored!

	assert "`vcetype'"=="cluster"

end
