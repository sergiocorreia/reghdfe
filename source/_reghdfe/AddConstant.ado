capture program drop AddConstant
program define AddConstant
	syntax varlist(numeric)
	foreach var of local varlist {
		local mean : char `var'[mean]
		assert "`mean'"!=""
		assert !missing(`mean')
		qui replace `var' = `var' + `mean'
	}
end
