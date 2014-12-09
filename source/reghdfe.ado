*! reghdfe 1.2.0 10Aug2014
*! By Sergio Correia (sergio.correia@duke.edu)
* (built from multiple source files using push.py)

cap pr drop reghdfe
program define reghdfe
	local version `=clip(`c(version)', 11.2, 13.1)' // 11.2 minimum, 13+ preferred
	qui version `version'

	if replay() {
		if (`"`e(cmd)'"'!="reghdfe") error 301
		* Undocumented feature for debugging
		cap syntax, ALTernative
		if ("`alternative'"!="") {
			AlternativeCMD
			exit
		}
		else {	
			Replay `0'
		}
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

include "_reghdfe/Estimate.ado"
	include "_reghdfe/Parse.ado"
	include "_reghdfe/ExpandFactorVariables.ado"
	include "_reghdfe/EstimateDoF.ado"
		include "_reghdfe/ConnectedGroups.ado"
	include "_reghdfe/FixVarnames.ado"
	include "_reghdfe/Wrapper_regress.ado"
	include "_reghdfe/Wrapper_ivregress.ado"
	include "_reghdfe/Wrapper_ivreg2.ado"
include "_reghdfe/Replay.ado"
include "_reghdfe/AlternativeCMD.ado"
include "_common/Assert.ado"
include "_common/Debug.ado"
