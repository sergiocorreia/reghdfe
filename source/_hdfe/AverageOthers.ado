// -------------------------------------------------------------
// AvgE: Average of all the other obs in a group, except each obs itself
// -------------------------------------------------------------
cap pr drop AverageOthers
program AverageOthers , sortpreserve
syntax varname , BY(varlist) Generate(name) [EXCLUDESELF]

* EXCLUDESELF: Excludes obs at hand when computing avg

***[EXCLUDE(varname)]
*** Do not use obs where `exclude'!=0 to compute the means, but do fill out these values

* Alternative:
* MeanOthers = MeanAll * N/(N-1) - X / (N-1) = (SumAll-X)/(N-1)
* Also, using mean() instead of total() would give less rounding errors

	sort `by'

	conf new var `generate'
	***if ("`exclude'"!="") local cond " if !`exclude'"
	
	* Calculate avg by group
	tempvar total count
	qui gen double `generate' = `varlist' `cond'
	
	* Sum
	*qui by `by' : egen double `generate' = mean(`var')
	qui by `by' : gen double `total' = sum(`generate')
	qui by `by' : replace `total' = `total'[_N]
	
	* Count
	qui by `by' : gen double `count' = sum(`generate'<.)
	qui by `by' : replace `count' = `count'[_N]
	
	* Substract itself
	if ("`excludeself'"!="") qui by `by' : replace `total' = `total' - `generate' if (`generate'<.)
	if ("`excludeself'"!="") qui by `by' : replace `count' = `count' - 1 if (`generate'<.)
	
	* Divide
	qui replace `generate' = `total' / `count'
	
	**qui by `by' : replace `generate' = `generate'[_N]
	
	* Adjust negative values b/c of rounding errors introduced by -excludeself- (risky)
	if ("`excludeself'"!="") {
		replace `generate' = 0 if inrange(`generate', -1e-8, 0)
		local note X
	}

	* Add label and chars
	local name = subinstr("`by'", " ", "_", .)
	char `generate'[avge_equation]  AvgE`note'
	char `generate'[name] `name'
	char `generate'[depvar] `varlist'
	la var `generate' "Avg`note'. of `varlist' by `by'"
end
