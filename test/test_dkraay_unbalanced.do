* ==============================================================================
* Test Driscoll-Kraay Standard Errors with Unbalanced Panel
* ==============================================================================
* This do-file tests vce(dkraay) on a deeply imbalanced panel dataset
* and compares results across different implementations using esttab

cls
clear all
set more off

* Check for required packages
cap which estout
if _rc {
    di as error "estout/esttab not installed. Install with: ssc install estout"
    exit 199
}

* ------------------------------------------------------------------------------
* Setup: Create a deeply unbalanced panel dataset
* ------------------------------------------------------------------------------

clear
set seed 54321
set obs 500

* Generate unbalanced panel structure
* - 50 individuals with varying number of time periods
* - Some have only 2 periods, others have up to 20

gen int id = .
gen int time = .

* Create unbalanced structure
local obs = 1
forvalues i = 1/50 {
    * Random number of time periods for each individual (between 2 and 20)
    local n_periods = ceil(runiform() * 18) + 2
    
    * Random subset of time periods from 1-25
    forvalues t = 1/`n_periods' {
        if `obs' > 500 continue, break
        qui replace id = `i' in `obs'
        qui replace time = ceil(runiform() * 25) in `obs'
        local ++obs
    }
}

* Remove observations with missing id
drop if missing(id)

* Remove duplicate id-time combinations
duplicates drop id time, force

* Sort and verify panel structure
sort id time
xtset id time

di _n(2)
di as text "{hline 70}"
di as text "Panel Structure Summary (Unbalanced)"
di as text "{hline 70}"
xtdescribe
xtsum id time

* Show distribution of observations per individual
preserve
collapse (count) n_obs = time, by(id)
tabstat n_obs, stats(min p25 p50 p75 max mean sd) columns(statistics)
restore

* Generate variables
gen x1 = rnormal()
gen x2 = rnormal()

* Generate fixed effects
tempvar tmp_fe_id tmp_fe_time
gen `tmp_fe_id' = rnormal()
bysort id (time): replace `tmp_fe_id' = `tmp_fe_id'[1]
gen fe_id = `tmp_fe_id'

gen `tmp_fe_time' = rnormal()
bysort time (id): replace `tmp_fe_time' = `tmp_fe_time'[1]
gen fe_time = `tmp_fe_time'

* Generate outcome: y = 1 + 2*x1 + 0.5*x2 + fe_id + fe_time + epsilon
gen epsilon = rnormal()
gen y = 1 + 2*x1 + 0.5*x2 + fe_id + fe_time + epsilon

* Verify no missing values
assert !missing(y)
count
di as text "Total observations: " r(N)

* Re-establish panel structure
sort id time
xtset id time

di _n(2)
di as text "{hline 70}"
di as text "Number of time periods (for bandwidth selection): " 
di as result r(tmax) - r(tmin) + 1
di as text "{hline 70}"

* ==============================================================================
* Run regressions with different VCE methods
* ==============================================================================

estimates clear

* ------------------------------------------------------------------------------
* 1. reghdfe with different VCE options
* ------------------------------------------------------------------------------

di _n(2)
di as text "{hline 70}"
di as text "Running reghdfe estimations..."
di as text "{hline 70}"

* Robust SE
qui reghdfe y x1 x2, absorb(id) vce(robust)
estimates store reghdfe_robust

* Cluster on time
qui reghdfe y x1 x2, absorb(id) vce(cluster time)
estimates store reghdfe_cl_time

* Cluster on id
qui reghdfe y x1 x2, absorb(id) vce(cluster id)
estimates store reghdfe_cl_id

* Driscoll-Kraay with bandwidth 2 (1 lag)
qui reghdfe y x1 x2, absorb(id) vce(dkraay 2)
estimates store reghdfe_dk2

* Driscoll-Kraay with bandwidth 3 (2 lags)
qui reghdfe y x1 x2, absorb(id) vce(dkraay 3)
estimates store reghdfe_dk3

* Driscoll-Kraay with bandwidth 4 (3 lags)
qui reghdfe y x1 x2, absorb(id) vce(dkraay 4)
estimates store reghdfe_dk4

* Driscoll-Kraay with default bandwidth
qui reghdfe y x1 x2, absorb(id) vce(dkraay)
estimates store reghdfe_dk_default

* ------------------------------------------------------------------------------
* 2. xtscc (if installed)
* ------------------------------------------------------------------------------

cap which xtscc
if !_rc {
    di as text "Running xtscc estimations..."
    
    * xtscc with different lags
    qui xtscc y x1 x2, fe lag(1)
    estimates store xtscc_lag1
    
    qui xtscc y x1 x2, fe lag(2)
    estimates store xtscc_lag2
    
    qui xtscc y x1 x2, fe lag(3)
    estimates store xtscc_lag3
}
else {
    di as text "xtscc not installed - skipping"
}

* ------------------------------------------------------------------------------
* 3. ivreghdfe (if installed)
* ------------------------------------------------------------------------------

cap which ivreghdfe
if !_rc {
    di as text "Running ivreghdfe estimations..."
    
    * ivreghdfe with dkraay option
    qui ivreghdfe y x1 x2, absorb(id) dkraay(2)
    estimates store ivreghdfe_dk2
    
    qui ivreghdfe y x1 x2, absorb(id) dkraay(3)
    estimates store ivreghdfe_dk3
    
    qui ivreghdfe y x1 x2, absorb(id) dkraay(4)
    estimates store ivreghdfe_dk4
}
else {
    di as text "ivreghdfe not installed - skipping"
}

* ==============================================================================
* Display comparison tables using esttab
* ==============================================================================

di _n(3)
di as text "{hline 70}"
di as text "COMPARISON TABLES"
di as text "{hline 70}"

* ------------------------------------------------------------------------------
* Table 1: reghdfe - Different VCE methods
* ------------------------------------------------------------------------------

di _n(2)
di as text "{bf:Table 1: reghdfe with different VCE specifications}"
di as text "Coefficients should be identical; standard errors vary by method"
di as text "{hline 70}"

esttab reghdfe_robust reghdfe_cl_time reghdfe_cl_id reghdfe_dk2 reghdfe_dk3 reghdfe_dk4, ///
    se(%9.4f) b(%9.4f) ///
    mtitles("Robust" "Cl(time)" "Cl(id)" "DK(2)" "DK(3)" "DK(4)") ///
    title("reghdfe: Comparison of VCE methods") ///
    stats(N r2_within, fmt(%9.0fc %9.4f) labels("Observations" "Within R2")) ///
    star(* 0.10 ** 0.05 *** 0.01) ///
    note("DK(#) = Driscoll-Kraay with bandwidth=#")

* ------------------------------------------------------------------------------
* Table 2: Driscoll-Kraay across packages
* ------------------------------------------------------------------------------

* Check which estimates are available
local dk_estimates "reghdfe_dk2 reghdfe_dk3 reghdfe_dk4"

cap estimates describe xtscc_lag1
if !_rc {
    local dk_estimates "`dk_estimates' xtscc_lag1 xtscc_lag2 xtscc_lag3"
}

cap estimates describe ivreghdfe_dk2
if !_rc {
    local dk_estimates "`dk_estimates' ivreghdfe_dk2 ivreghdfe_dk3 ivreghdfe_dk4"
}

di _n(2)
di as text "{bf:Table 2: Driscoll-Kraay across different packages}"
di as text "Comparing reghdfe, xtscc, and ivreghdfe implementations"
di as text "{hline 70}"

* Build the esttab command dynamically based on available estimates
cap estimates describe xtscc_lag1
local has_xtscc = !_rc

cap estimates describe ivreghdfe_dk2
local has_ivreghdfe = !_rc

if `has_xtscc' & `has_ivreghdfe' {
    esttab reghdfe_dk2 xtscc_lag1 ivreghdfe_dk2 reghdfe_dk3 xtscc_lag2 ivreghdfe_dk3, ///
        se(%9.4f) b(%9.4f) ///
        mtitles("reghdfe" "xtscc" "ivreghdfe" "reghdfe" "xtscc" "ivreghdfe") ///
        title("Driscoll-Kraay: Cross-package comparison") ///
        mgroups("Bandwidth=2 (1 lag)" "Bandwidth=3 (2 lags)", pattern(1 0 0 1 0 0)) ///
        stats(N, fmt(%9.0fc) labels("Observations")) ///
        star(* 0.10 ** 0.05 *** 0.01)
}
else if `has_xtscc' {
    esttab reghdfe_dk2 xtscc_lag1 reghdfe_dk3 xtscc_lag2, ///
        se(%9.4f) b(%9.4f) ///
        mtitles("reghdfe" "xtscc" "reghdfe" "xtscc") ///
        title("Driscoll-Kraay: reghdfe vs xtscc") ///
        mgroups("Bandwidth=2 (1 lag)" "Bandwidth=3 (2 lags)", pattern(1 0 1 0)) ///
        stats(N, fmt(%9.0fc) labels("Observations")) ///
        star(* 0.10 ** 0.05 *** 0.01)
}
else if `has_ivreghdfe' {
    esttab reghdfe_dk2 ivreghdfe_dk2 reghdfe_dk3 ivreghdfe_dk3, ///
        se(%9.4f) b(%9.4f) ///
        mtitles("reghdfe" "ivreghdfe" "reghdfe" "ivreghdfe") ///
        title("Driscoll-Kraay: reghdfe vs ivreghdfe") ///
        mgroups("Bandwidth=2 (1 lag)" "Bandwidth=3 (2 lags)", pattern(1 0 1 0)) ///
        stats(N, fmt(%9.0fc) labels("Observations")) ///
        star(* 0.10 ** 0.05 *** 0.01)
}
else {
    esttab reghdfe_dk2 reghdfe_dk3 reghdfe_dk4 reghdfe_dk_default, ///
        se(%9.4f) b(%9.4f) ///
        mtitles("DK(2)" "DK(3)" "DK(4)" "DK(default)") ///
        title("Driscoll-Kraay: reghdfe with different bandwidths") ///
        stats(N, fmt(%9.0fc) labels("Observations")) ///
        star(* 0.10 ** 0.05 *** 0.01)
}

* ------------------------------------------------------------------------------
* Table 3: Detailed SE comparison (focus on x1 coefficient)
* ------------------------------------------------------------------------------

di _n(2)
di as text "{bf:Table 3: Standard Errors for x1 coefficient}"
di as text "Direct comparison of SE magnitudes across methods"
di as text "{hline 70}"

* Extract SEs manually for clearer comparison
tempname se_table
matrix `se_table' = J(7, 2, .)
matrix colnames `se_table' = "x1_SE" "x2_SE"
matrix rownames `se_table' = "Robust" "Cl(time)" "Cl(id)" "DK(2)" "DK(3)" "DK(4)" "DK(default)"

local row = 1
foreach est in reghdfe_robust reghdfe_cl_time reghdfe_cl_id reghdfe_dk2 reghdfe_dk3 reghdfe_dk4 reghdfe_dk_default {
    qui estimates restore `est'
    matrix `se_table'[`row', 1] = _se[x1]
    matrix `se_table'[`row', 2] = _se[x2]
    local ++row
}

matlist `se_table', format(%9.4f) title("Standard Errors by VCE Method")

* ==============================================================================
* Summary
* ==============================================================================

di _n(3)
di as text "{hline 70}"
di as text "TEST COMPLETED"
di as text "{hline 70}"
di as text ""
di as text "Key observations:"
di as text "1. Coefficients should be identical across all methods"
di as text "2. Standard errors vary by method and bandwidth"
di as text "3. reghdfe's vce(dkraay #) should match ivreghdfe's dkraay(#)"
di as text "4. xtscc uses lag(#) where # = bandwidth - 1"
di as text ""

estimates dir

