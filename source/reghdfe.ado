*! reghdfe 1.4.420 20mar2015
*! Sergio Correia (sergio.correia@duke.edu)
* (built from multiple source files using build.py)

cap pr drop reghdfe
program define reghdfe
	local version `=clip(`c(version)', 11.2, 13.1)' // 11.2 minimum, 13+ preferred
	qui version `version'

	if replay() {
		if (`"`e(cmd)'"'!="reghdfe") error 301
		Replay `0'
	}
	else {
		* Estimate, and then clean up Mata in case of failure
		mata: st_global("reghdfe_pwd",pwd())
		cap noi Estimate `0'
		if (_rc) {
			local rc = _rc
			reghdfe_absorb, step(stop)
			exit `rc'
		}
	}
end

* Note: Assert and Debug must go first
include "_common/Assert.ado"
include "_common/Debug.ado"
include "_reghdfe_absorb/GenerateID.ado"

include "_mata/fix_psd.mata"
include "_reghdfe/Estimate.ado"
	include "_reghdfe/Parse.ado"
	include "_reghdfe/DropSingletons.ado"
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


