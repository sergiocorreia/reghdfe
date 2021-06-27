* ===========================================================================
* Create toy patent dataset loosely based on 'auto' dataset
* ===========================================================================

	clear all
	cls


// --------------------------------------------------------------------------
// Test dataset
// --------------------------------------------------------------------------

	set seed 1234
	sysuse auto, clear
	rename price lab_size
	rename weight funding
	rename turn year
	keep lab_size funding year

	replace lab_size = lab_size / 1000
	replace funding = funding / 1000

	*keep in 1/50 // 50 patents only
	gen int patent_id = _n
	gen int n = ceil(runiform()*10) // up to 10 inventors per parent
	expand n
	drop n

	sort patent_id, stable
	by patent_id: gen long inventor_id = ceil(runiform()*20) // up to 20 different inventors
	sort patent_id inventor_id, stable
	by patent_id inventor_id: drop if _N > 1 // drop obs with two inventors in one patent
	
	sort inventor_id, stable
	gen double ability = exp(rnormal())
	by inventor_id: replace ability = ability[1]

	gegen sum_ability = total(ability), by(patent_id)
	gen y = 0.2 * rnormal() + 0.1 * sum_ability + 3 * funding - 1 * lab_size
	sort patent_id, stable
	by patent_id: replace y = y[1]
	*gen citations = ceil(exp(y))
	*gen log_citations = log(1+citations)
	*assert !mi(citations)
	assert !mi(y)
	rename y citations
	drop ability sum_ability

	gen byte c = 1
	gen double w1 = 10 * runiform() // for aweights and pweights
	sort patent_id, stable
	by patent_id: replace w1 = w1[1]
	gen double w2 = ceil(w1) // for fweights

	sort inventor_id year, stable
	gen age = ceil(18+runiform()*30)
	by inventor_id: replace age = age[1] + year - year[1]
	su age


	sort patent_id inventor_id, stable

	preserve
		* Wide-data alternative
		tab inventor_id, gen(dummy)
		sort patent_id, stable
		* THIS IS WRONG: by patent_id: keep if _n == 1
		* THIS IS RIGHT:
		gcollapse (first) citations lab_size funding year c w1 w2 (sum) dummy*, by(patent_id) fast
		areg citations lab_size funding dummy*, a(year)
		save "toy-patents-wide", replace
	restore


	preserve
		* Wide-data alternative
		tab inventor_id, gen(dummy)
		sort patent_id, stable
		gcollapse (first) citations lab_size funding year c w1 w2 (sum) dummy* (count) n=inventor_id, by(patent_id) fast
		foreach dummy of varlist dummy* {
			replace `dummy' = `dummy' / n
		}
		areg citations lab_size funding dummy*, a(year)
		save "toy-patents-wide-average", replace
	restore


	save "toy-patents-long", replace
