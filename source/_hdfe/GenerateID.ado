// -------------------------------------------------------------
// Faster alternative to -egen group-. MVs, IF, etc not allowed!
// -------------------------------------------------------------

cap pr drop GenerateID
program define GenerateID, sortpreserve
syntax varlist(numeric) , [REPLACE Generate(name)]

	assert ("`replace'"!="") + ("`generate'"!="") == 1
	// replace XOR generate, could also use -opts_exclusive -
	foreach var of varlist `varlist' {
		assert !missing(`var')
	}

	local numvars : word count `varlist'
	if ("`replace'"!="") assert `numvars'==1 // Can't replace more than one var!
	
	// Create ID
	tempvar new_id
	sort `varlist'
	by `varlist': gen long `new_id' = (_n==1)
	qui replace `new_id' = sum(`new_id')
	qui compress `new_id'
	assert !missing(`new_id')
	
	local name = "i." + subinstr("`varlist'", " ", "#i.", .)
	char `new_id'[name] `name'
	la var `new_id' "[ID] `name'"

	// Either replace or generate
	if ("`replace'"!="") {
		drop `varlist'
		rename `new_id' `varlist'
	}
	else {
		rename `new_id' `generate'
	}

end
