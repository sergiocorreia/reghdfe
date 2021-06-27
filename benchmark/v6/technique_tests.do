
// Testing performance of LSMR, LSQR, MAP 

clear all
cls 

loc dir "/scratch/m1nbc00/flex_fe"
loc charts "`dir'/misc_charts"

// Testing specs
loc total_iter = 100
loc step = 5 

loc tol_start = 1
loc tol_step = 1
loc tol_stop = 14

use "/ofs/scratch2/m1txp01/moodies_clo/collateral_folder/collat_Oct_2020short"

gegen issuer_id = group(issuer_name)
keep y x1 issuer_id deal_id year error 

// trim it down a bit -- too slow otherwise 
// still keep unbalanced nature 
bys issuer_id year: keep if _n < 10

gen t = date(run_date, "YMD")
format t %td 
gen year = year(t)

loc fe1 "issuer_id"
loc fe2 "deal_id"
loc fe3 "year"

loc fes "`fe1' `fe2'#`fe3'"

gen double error  = runiform()
gen double x1 = rnormal()
gen double y = log(`fe1') + sin(`fe2'*`fe3') + 5*x1 + 2*error 

loc lhs "y"
loc rhs "x1"

loc cmd "groupreg `lhs' `rhs', a(`fes')"

// Testing relationship between tolerance & precision of beta hat

gen beta_hat = .
gen tech = ""
gen tol = .

loc i = 1
foreach tech in lsmr lsqr map {

	forvalues tol = `tol_start'(`tol_step')`tol_stop' {

		`cmd' tech(`tech') tol(1e-`tol')

		replace beta_hat = _b[`rhs'] if _n == `i'
		replace tech = "`tech'" if _n == `i'
		replace tol = `tol' if _n == `i'

		local ++i
	}
}

preserve
	
	keep beta_hat tech tol
	drop if mi(beta_hat)

	gen beta = beta_hat if tech == "lsmr" & tol == `tol_stop'
	egen true_beta = max(beta)

	gen precision = log(abs(beta_hat - true_beta))

	save "`charts'/tol.dta" // in case we need to revist these #s

	tw (scatter precision tol if tech == "lsmr", mcol(blue)) ///
		(scatter precision tol if tech == "lsqr", mcol(red)) ///
		(scatter precision tol if tech == "map", mcol(gs6)), ///
		scheme(s1mono) title("Precision vs. Tolerance by Method") ///
		ytitle("Log(|error|)") xtitle("Tolerance") ///
		legend(label(1 "LSMR") label(2 "LSQR") label(3 "MAP"))
	graph export "`charts'/tolerance.png", replace 

restore 

drop beta_hat tech tol

// Testing relationship between iterations & precision of beta hat

gen beta_hat = .
gen tech = ""
gen iter = .

loc i = 1
foreach tech in lsmr lsqr map {

	forvalues iter = 1(`step')`total_iter'+1 {

		`cmd' tech(`tech') maxiter(`iter') abort(0) tol(1e-14)

		replace beta_hat = _b[`rhs'] if _n == `i'
		replace tech = "`tech'" if _n == `i'
		replace iter = `iter' if _n == `i'

		local ++i
	}
}

preserve
	
	keep beta_hat tech iter
	drop if mi(beta_hat)

	gen beta = beta_hat if tech == "lsmr" & iter == 101
	egen true_beta = max(beta)

	gen precision = log(abs(beta_hat - true_beta))

	save "`charts'/iter.dta"

	tw (scatter precision iter if tech == "lsmr", mcol(blue)) ///
		(scatter precision iter if tech == "lsqr", mcol(red)) ///
		(scatter precision iter if tech == "map", mcol(gs6)), ///
		scheme(s1mono) title("Precision vs. Iterations by Method") ///
		ytitle("Log(|error|)") xtitle("Iterations") ///
		legend(label(1 "LSMR") label(2 "LSQR") label(3 "MAP"))
	graph export "`charts'/iterations.png", replace 

restore 
