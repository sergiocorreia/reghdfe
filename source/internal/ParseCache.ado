capture program drop ParseCache
program define ParseCache, sclass
	syntax, [CACHE(string)] [IFIN(string) ABSORB(string) VCE(string)] 
	if ("`cache'"!="") {
		local 0 `cache'
		syntax name(name=opt id="cache option"), [KEEPvars(varlist)]
		Assert inlist("`opt'", "save", "use"), msg("invalid cache option {cmd`opt'}") // -clear- is also a valid option but intercepted earlier
	}

	local savecache = ("`opt'"=="save")
	local usecache = ("`opt'"=="use")
	local is_cache : char _dta[reghdfe_cache]

	* Sanity checks on usecache
	if (`usecache') {
		local cache_obs : char _dta[cache_obs]
		local cache_absorb : char _dta[absorb]
		local cache_vce : char _dta[vce]

		Assert "`is_cache'"=="1" , msg("cache(use) requires a previous cache(save) operation")
		Assert `cache_obs'==`c(N)', msg("dataset cannot change after cache(save)")
		Assert "`cache_absorb'"=="`absorb'", msg("cached dataset has different absorb()")
		Assert "`ifin'"=="", msg("cannot use if/in with cache(use); data has already been transformed")
		Assert "`cache_vce'"=="`vce'", msg("cached dataset has a different vce()")
	}
	else {
		Assert "`is_cache'"!="1", msg("reghdfe error: data transformed with cache(save) requires cache(use)")
	}
	
	if (!`savecache') Assert "`keepvars'"=="", msg("reghdfe error: {cmd:keepvars()} suboption requires {cmd:cache(save)}")

	local keys savecache keepvars usecache
	foreach key of local keys {
		sreturn local `key' ``key''
	}
end
