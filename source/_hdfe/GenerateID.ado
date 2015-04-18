// -------------------------------------------------------------
// Faster alternative to -egen group-. MVs, IF, etc not allowed!
// -------------------------------------------------------------

cap pr drop GenerateID
program define GenerateID, sortpreserve
syntax varlist(numeric) , [REPLACE Generate(name)] [CLUSTERVARS(namelist) NESTED]
assert ("`replace'"!="") + ("`generate'"!="") == 1

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

	// Could use these chars to speed up DropSingletons	and Wrapper_mwc
	*char `new_id'[obs] `c(N)' 
	*char `new_id'[id] 1 

	// Either replace or generate
	if ("`replace'"!="") {
		drop `varlist'
		rename `new_id' `varlist'
		local new_id `varlist' // I need to keep track of the variable for the clustervar part
	}
	else {
		rename `new_id' `generate'
		local new_id `generate'
	}

	// See if var. is nested within a clustervar
	local in_clustervar 0
	local is_clustervar 0

	if ("`clustervars'"!="") {
		
		* Check if clustervar===ID
		foreach clustervar of local clustervars {
			if ("`new_id'"=="`clustervar'") {
				local is_clustervar 1
				local nesting_clustervar "`clustervar'"
				continue, break
			}
		}
		
		* Check if ID is nested within cluster ("if two obs. belong to the same ID, they belong to the same cluster")
		if (!`is_clustervar' & "`nested'"!="") {
			tempvar same
			qui gen byte `same' = .
			foreach clustervar of local clustervars {

				* Avoid check if clustervar is another absvar
				* Reason: it would be stupid to have one absvar nested in another (same result as dropping nesting one)
				local clustervar_is_absvar = regexm("`clustervar'","__FE[0-9]+__")
				if (`clustervar_is_absvar') continue

				qui bys `new_id' (`clustervar'): replace `same' = (`clustervar'[1] == `clustervar'[_N])
				qui cou if (`same'==0)
				if r(N)==0 {
					local in_clustervar 1
					local nesting_clustervar "`clustervar'"
					continue, break
				}
			}
		}
	}

	char `new_id'[is_clustervar] `is_clustervar'
	char `new_id'[in_clustervar] `in_clustervar'
	char `new_id'[nesting_clustervar] `nesting_clustervar' 
end
