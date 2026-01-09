* ==============================================================================
* Test Driscoll-Kraay Standard Errors in reghdfe
* ==============================================================================
* This do-file compares reghdfe's vce(dkraay) option against other implementations

cls
clear all
set more off

* ------------------------------------------------------------------------------
* Setup: Create a panel dataset from auto
* ------------------------------------------------------------------------------

* Load auto and create a fake panel structure
sysuse auto, clear

* Create a panel structure: 
* - Use rep78 as time variable (has 5 levels: 1-5)
* - Use foreign as panel id (has 2 levels: 0, 1)
* We need to reshape or create fake panel; let's create one

* Alternative approach: create a proper panel dataset
clear
set seed 12345
set obs 200

* Generate panel structure
gen int id = ceil(_n / 10)      // 20 individuals
gen int time = mod(_n - 1, 10) + 1  // 10 time periods each
sort id time

* Generate variables
gen x1 = rnormal()
gen x2 = rnormal()

* Generate fixed effects properly
* fe_id: one random value per individual, constant across time
tempvar tmp_fe_id
gen `tmp_fe_id' = rnormal() if time == 1
bysort id (time): replace `tmp_fe_id' = `tmp_fe_id'[1]
gen fe_id = `tmp_fe_id'

* fe_time: one random value per time period, constant across individuals  
tempvar tmp_fe_time
gen `tmp_fe_time' = rnormal() if id == 1
bysort time (id): replace `tmp_fe_time' = `tmp_fe_time'[1]
gen fe_time = `tmp_fe_time'

* Re-sort by panel structure
sort id time

* Generate outcome with known DGP: y = 1 + 2*x1 + 0.5*x2 + fe_id + fe_time + epsilon
gen epsilon = rnormal()
gen y = 1 + 2*x1 + 0.5*x2 + fe_id + fe_time + epsilon

* Verify no missing values
assert !missing(y)
assert !missing(fe_id)
assert !missing(fe_time)

* Set as panel
xtset id time

di _n(3)
di as text "{hline 80}"
di as text "Dataset Summary"
di as text "{hline 80}"
describe
xtsum y x1 x2

* ==============================================================================
* Test 1: Basic reghdfe with vce(dkraay)
* ==============================================================================

di _n(3)
di as text "{hline 80}"
di as text "TEST 1: reghdfe with vce(dkraay) - default bandwidth"
di as text "{hline 80}"

set trace off
reghdfe y x1 x2, absorb(id) vce(dkraay)
estimates store reghdfe_dk_default


di _n(3)
di as text "{hline 80}"
di as text "TEST 2: reghdfe with vce(dkraay 2) - explicit bandwidth = 2"
di as text "{hline 80}"

reghdfe y x1 x2, absorb(id) vce(dkraay 2)
estimates store reghdfe_dk_2

di _n(3)
di as text "{hline 80}"
di as text "TEST 3: reghdfe with vce(dkraay 4) - explicit bandwidth = 4"
di as text "{hline 80}"

reghdfe y x1 x2, absorb(id) vce(dkraay 4)
estimates store reghdfe_dk_4

* ==============================================================================
* Test 2: Compare with xtscc (if installed)
* ==============================================================================

di _n(3)
di as text "{hline 80}"
di as text "COMPARISON: xtscc (Hoechle's implementation)"
di as text "{hline 80}"

cap which xtscc
if _rc {
    di as error "xtscc not installed. Install with: ssc install xtscc"
    di as text "Skipping xtscc comparison..."
}
else {
    * xtscc with default lag
    xtscc y x1 x2, fe
    estimates store xtscc_default
    
    * xtscc with lag(2)
    xtscc y x1 x2, fe lag(2)
    estimates store xtscc_lag2
    
    * xtscc with lag(4)
    xtscc y x1 x2, fe lag(4)
    estimates store xtscc_lag4
}

* ==============================================================================
* Test 3: Compare with ivreg2/ivreghdfe (if installed)
* ==============================================================================

di _n(3)
di as text "{hline 80}"
di as text "COMPARISON: ivreg2 with dkraay option"
di as text "{hline 80}"

cap which ivreg2
if _rc {
    di as error "ivreg2 not installed. Install with: ssc install ivreg2"
    di as text "Skipping ivreg2 comparison..."
}
else {
    * ivreg2 requires the bw() option for dkraay
    * Note: ivreg2's dkraay uses bw=#, which is bandwidth = lags + 1
    
    * With 2 lags (bw=3)
    ivreg2 y x1 x2 i.id, dkraay(3) small
    estimates store ivreg2_dk3
    
    * With 4 lags (bw=5)  
    ivreg2 y x1 x2 i.id, dkraay(5) small
    estimates store ivreg2_dk5
}

* ==============================================================================
* Test 4: Compare with ivreghdfe (if installed)
* ==============================================================================

di _n(3)
di as text "{hline 80}"
di as text "COMPARISON: ivreghdfe with dkraay option (if supported)"
di as text "{hline 80}"

cap which ivreghdfe
if _rc {
    di as error "ivreghdfe not installed"
    di as text "Skipping ivreghdfe comparison..."
}
else {
    * Note: Check if ivreghdfe supports dkraay
    cap ivreghdfe y x1 x2, absorb(id) cluster(time)
    if !_rc {
        estimates store ivreghdfe_cluster_time
    }
}

* ==============================================================================
* Test 5: Compare with reghdfe vce(cluster time)
* ==============================================================================

di _n(3)
di as text "{hline 80}"
di as text "COMPARISON: reghdfe with vce(cluster time)"
di as text "Note: Clustering on time is the no-lag version of Driscoll-Kraay"
di as text "{hline 80}"

reghdfe y x1 x2, absorb(id) vce(cluster time)
estimates store reghdfe_cluster_time

* ==============================================================================
* Summary: Compare all estimates
* ==============================================================================

di _n(3)
di as text "{hline 80}"
di as text "SUMMARY: Comparison of all estimates"
di as text "{hline 80}"

* Check which estimates exist and compare
estimates dir

di _n
di as text "Comparing coefficient estimates and standard errors:"
di as text "(All methods should give same coefficients, but different SEs)"
di _n

* Display comparison table
cap estimates table reghdfe_dk_default reghdfe_dk_2 reghdfe_dk_4 reghdfe_cluster_time, ///
    b(%9.4f) se(%9.4f) stats(N r2 df_r) ///
    title("reghdfe: Driscoll-Kraay vs Cluster on Time")

di _n
cap estimates table xtscc_default xtscc_lag2 xtscc_lag4, ///
    b(%9.4f) se(%9.4f) stats(N r2) ///
    title("xtscc: Various lag specifications")

di _n
cap estimates table ivreg2_dk3 ivreg2_dk5, ///
    b(%9.4f) se(%9.4f) stats(N) ///
    title("ivreg2: Driscoll-Kraay with different bandwidths")

* ==============================================================================
* Test 6: Edge cases
* ==============================================================================

di _n(3)
di as text "{hline 80}"
di as text "EDGE CASES"
di as text "{hline 80}"

* Test with bandwidth = 1 (0 lags, should be equivalent to clustering on time)
di as text _n "Test: vce(dkraay 1) - should match vce(cluster time)"
reghdfe y x1 x2, absorb(id) vce(dkraay 1)
estimates store reghdfe_dk_1

estimates table reghdfe_dk_1 reghdfe_cluster_time, b(%9.6f) se(%9.6f) ///
    title("Comparison: dkraay(1) vs cluster(time)")

* Test with two-way fixed effects
di as text _n "Test: Two-way FE with Driscoll-Kraay"
reghdfe y x1 x2, absorb(id time) vce(dkraay 2)
estimates store reghdfe_twoway_dk

* Test with weights
di as text _n "Test: With weights"
gen weight = runiform() + 0.5
reghdfe y x1 x2 [aw=weight], absorb(id) vce(dkraay 2)
estimates store reghdfe_dk_weighted

* ==============================================================================
* Cleanup
* ==============================================================================

di _n(3)
di as text "{hline 80}"
di as text "All tests completed!"
di as text "{hline 80}"

estimates dir

