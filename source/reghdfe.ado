*! reghdfe VERSION_NUMBER
*! Sergio Correia (sergio.correia@duke.edu)
* (built from multiple source files using build.py)

include _mata/reghdfe.mata

cap pr drop reghdfe
program define reghdfe
	local version `=clip(`c(version)', 11.2, 13.1)' // 11.2 minimum, 13+ preferred
	qui version `version'

	* Intercept version calls
	cap syntax, version
	local rc = _rc
	 if (`rc'==0) {
		Version
		exit
	}

	* Intercept multiprocessor/parallel calls
	cap syntax, instance [*]
	local rc = _rc
	 if (`rc'==0) {
		ParallelInstance, `options'
		exit
	}

	if replay() {
		if (`"`e(cmd)'"'!="reghdfe") error 301
		Replay `0'
	}
	else {
		* Estimate, and then clean up Mata in case of failure
		mata: st_global("reghdfe_pwd",pwd())
		Stop // clean leftovers for a possible [break]
		cap noi Estimate `0'
		if (_rc) {
			local rc = _rc
			Stop
			exit `rc'
		}
	}
end

* Note: Assert and Debug must go first
include "_common/Assert.ado"
include "_common/Debug.ado"
include "_common/Version.ado"

include "_mata/fix_psd.mata"

include "_reghdfe/Estimate.ado"
include "_reghdfe/Parse.ado"
include "_hdfe/DropSingletons.ado"
include "_reghdfe/ExpandFactorVariables.ado"
include "_reghdfe/FixVarnames.ado"
include "_reghdfe/Wrapper_regress.ado"
include "_reghdfe/Wrapper_mwc.ado"
include "_reghdfe/Wrapper_avar.ado"
include "_reghdfe/Wrapper_ivregress.ado"
include "_reghdfe/Wrapper_ivreg2.ado"
include "_reghdfe/AddConstant.ado"
include "_reghdfe/Attach.ado"
include "_reghdfe/Replay.ado"
include "_reghdfe/Header.ado"

include "_hdfe/ConnectedGroups.ado"
include "_hdfe/GenerateID.ado"
include "_hdfe/AverageOthers.ado"
include "_hdfe/EstimateDoF.ado"
include "_hdfe/Start.ado"
include "_hdfe/ParseOneAbsvar.ado"
include "_hdfe/Precompute.ado"
include "_hdfe/Demean.ado"
include "_hdfe/DemeanParallel.ado"
include "_hdfe/ParallelInstance.ado"
include "_hdfe/Save.ado"
include "_hdfe/Stop.ado"
include "_hdfe/CheckCorrectOrder.ado"
