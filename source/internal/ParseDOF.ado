capture program drop ParseDOF
program define ParseDOF, sclass
	syntax, [ALL NONE] [PAIRwise TWO THREE] [CLusters] [CONTinuous]
	local opts `pairwise' `two' `three' `clusters' `continuous'
	local n : word count `opts'
	local first_opt : word 1 of `opt'

	opts_exclusive "`all' `none'" dofadjustments
	opts_exclusive "`pairwise' `two' `three'" dofadjustments
	opts_exclusive "`all' `first_opt'" dofadjustments
	opts_exclusive "`none' `first_opt'" dofadjustments

	if ("`none'" != "") local opts
	if ("`all'" != "") local opts pairwise clusters continuous

	if (`: list posof "three" in opts') {
		cap findfile group3hdfe.ado
		Assert !_rc , msg("error: -group3hdfe- not installed, please run {stata ssc install group3hdfe}")
	}

	sreturn local dofadjustments "`opts'"
end
