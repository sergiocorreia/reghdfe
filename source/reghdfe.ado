
// Mata code is first, then main reghdfe.ado, then auxiliary .ado files
clear mata
include "mata/map.mata"

capture program drop reghdfe
program define reghdfe

* Set Stata version
	version `=clip(c(version), 11.2, 14.1)'

* Intercept old+version
	cap syntax, version old
	if !c(rc) {
		reghdfe_old, version
		exit
	}

* Intercept version
	cap syntax, version [*]
	if !c(rc) {
		Version , `options'
		exit
	}

* Intercept old
	cap syntax anything(everything) [fw aw pw/], [*] old
	if !c(rc) {
		di as error "(running historical version of reghdfe)"
		if ("`weight'"!="") local weightexp [`weight'=`exp']
		reghdfe_old `anything' `weightexp', `options'
		exit
	}

* Intercept cache(clear) (must be before replay)
	local cache
	cap syntax, CACHE(string)
	if ("`cache'"=="clear") {
		cap mata: mata drop HDFE_S // overwrites c(rc)
		cap mata: mata drop varlist_cache
		cap mata: mata drop tss_cache
		cap global updated_clustervars
		cap matrix drop reghdfe_statsmatrix
		exit
	}

* Intercept replay
	if replay() {
		if (`"`e(cmd)'"'!="reghdfe") error 301
		if ("`0'"=="") local comma ","
		Replay `comma' `0' stored // also replays stored regressions (first stages, reduced, etc.)
		exit
	}

* Intercept cache(save)
	local cache
	cap syntax anything(everything) [fw aw pw/], [*] CACHE(string)
	if (strpos("`cache'", "save")==1) {
		cap noi InnerSaveCache `0'
		if (c(rc)) {
			local rc = c(rc)
			cap mata: mata drop HDFE_S // overwrites c(rc)
			cap mata: mata drop varlist_cache
			cap mata: mata drop tss_cache
			global updated_clustervars
			cap matrix drop reghdfe_statsmatrix
			exit `rc'
		}
		exit
	}

* Intercept cache(use)
	local cache
	cap syntax anything(everything) [fw aw pw/], [*] CACHE(string)
	if ("`cache'"=="use") {
		InnerUseCache `0'
		exit
	}

* Finally, call Inner if not intercepted before
	local is_cache : char _dta[reghdfe_cache]
	Assert ("`is_cache'"!="1"), msg("reghdfe error: data transformed with -savecache- requires option -usecache-")
	cap noi Inner `0'
	if (c(rc)) {
		local rc = c(rc)
		cap mata: mata drop HDFE_S // overwrites c(rc)
		exit `rc'
	}
end

// -------------------------------------------------------------------------------------------------
include "common/Assert.ado"
include "common/Debug.ado"
include "common/Version.ado"
include "common/Tic.ado"
include "common/Toc.ado"
include "internal/Inner.ado"
	include "internal/Parse.ado"
		include "internal/ParseCache.ado"
		include "internal/ParseIV.ado"
		include "internal/ParseStages.ado"
		include "internal/ParseVCE.ado"
		include "internal/ParseAbsvars.ado"
		include "internal/ParseDOF.ado"
		include "internal/ParseImplicit.ado"
	include "internal/GenUID.ado"
	include "internal/Compact.ado"
		include "internal/ExpandFactorVariables.ado"
	include "internal/Prepare.ado"
	include "internal/Stats.ado"
	include "internal/JointTest.ado"
	include "internal/Wrapper_regress.ado"
		include "internal/RemoveCollinear.ado"
	include "internal/Wrapper_avar.ado"
	include "internal/Wrapper_mwc.ado"
	include "internal/Wrapper_ivreg2.ado"
	include "internal/Wrapper_ivregress.ado"
		include "internal/GenerateID.ado"
	include "internal/SaveFE.ado"
	include "internal/Post.ado"
		include "internal/FixVarnames.ado"
		include "internal/Subtitle.ado"
	include "internal/Attach.ado"
include "internal/Replay.ado"
	include "internal/Header.ado"
include "internal/InnerSaveCache.ado"
include "internal/InnerUseCache.ado"
// -------------------------------------------------------------------------------------------------
