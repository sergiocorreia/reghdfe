* ===========================================================================
* Testing how fast (in iterations and tolerance) do LSMR, LSQR, and MAP converge
* ===========================================================================
* In most cases there is no need to test beyond 1e-12 as results have pretty much converged by then


	clear all
	cls
	set varabbrev off
	set type double

	set seed 123

	global data_path	"C:/Dropbox/Projects/Flexible_FEs/mwfe-datasets"
	global output_path	"C:/Git/groupreg/benchmarks/output"

	* Dependencies
	which scheme-sergio.scheme // https://github.com/sergiocorreia/stata-schemes

	* Options
	loc total_iter = 2
	loc step = 5 

	loc tol_start = 1
	loc tol_step = 1
	loc tol_stop = 15

	loc tol_start = 4
	loc tol_step = 1
	loc tol_stop = 9


// --------------------------------------------------------------------------
// Programs
// --------------------------------------------------------------------------

program define LoadDirectors
	use id1 id2 using "$data_path/directors/dataset", clear
	// set seed 1234
	gen double x1 = rnormal()
	gen double x2 = uniform()
	gen double y = 10 + x1 + x2 + sin(id1) + log(10+id2) + 100 * rnormal()
	assert !missing(y)

	* Remove singletons
	reghdfe y, a(id1 id2) dof(none)
	gen byte sample = e(sample)
	keep if sample
	drop sample

	* Standardize results so betas are invariant on scale of y/x
	foreach var of varlist y x1 x2 {
		su `var'
		replace `var' = `var' / r(sd)
	}
end


// --------------------------------------------------------------------------
// Test relationship between tolerance() & the precision of beta hat
// --------------------------------------------------------------------------

	clear frames
	loc rhs x1
	frame create results double(iter tol beta_map beta_lsmr beta_lsqr)
	loc cmd "groupreg y x1 x2, a(id1 id2) dof(none)"

	forval iter = 1/`total_iter' {
		LoadDirectors
		forval tol = `tol_start'(`tol_step')`tol_stop' {
			foreach tech in map lsmr lsqr {
				di as text "iter=`iter' tol=`tol' tech=`tech'"
				qui `cmd' tech(`tech') tol(1e-`tol') keepsingletons verbose(-1)
				loc beta_`tech' = _b[`rhs']
			}
			frame post results (`iter') (`tol') (`beta_map') (`beta_lsmr') (`beta_lsqr')
		}
	}

	* Plot results
	frame change results
	frame drop default


	save "borrar", replace
	cls

	use "borrar", clear

	egen double beta_avg = rowmean(beta_*)
	gegen double true_beta = total(beta_avg), by(iter tol)

	* su beta_avg if tol == `tol_stop', mean
	* *assert r(min) == r(max)
	* gcollapse (mean)
	* loc true_beta = `r(mean)'

	foreach tech in map lsmr lsqr {
		gen double precision_`tech' = log(abs(beta_`tech' - true_beta))
	}

	gcollapse (mean) precision_*, by(tol) fast


	tw 	(con precision_lsmr tol, msize(small) mcol(blue) lcol(blue)) ///
		(con precision_lsqr tol, msize(small) mcol(red)  lcol(red)) ///
		(con precision_map  tol, msize(small) mcol(gs6)  lcol(gs6)), ///
		scheme(sergio) ///
		title("Precision vs. Tolerance by Method") ///
		subtitle("(Avg over `total_iter' repetitions)")
		ytitle("Log(|error|)") xtitle("Tolerance") ///
		legend(order(1 "LSMR" 2 "LSQR" 3 "MAP"))

	graph export "$output_path/figures/tolerance.png", replace width(600)
	graph export "$output_path/figures/tolerance.pdf", replace

exit

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
