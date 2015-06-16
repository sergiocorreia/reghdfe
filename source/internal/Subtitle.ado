capture program drop Subtitle
program define Subtitle, eclass
	* Fill e(title3/4/5) based on the info of the other e(..)

	if (inlist("`e(vcetype)'", "Robust", "Cluster")) local hacsubtitle1 "heteroskedasticity"
	if ("`e(kernel)'"!="" & "`e(clustvar)'"=="") local hacsubtitle3 "autocorrelation"
	if ("`e(kiefer)'"!="") local hacsubtitle3 "within-cluster autocorrelation (Kiefer)"
	if ("`hacsubtitle1'"!="" & "`hacsubtitle3'" != "") local hacsubtitle2 " and "
	local hacsubtitle "`hacsubtitle1'`hacsubtitle2'`hacsubtitle3'"
	if strlen("`hacsubtitle'")>30 {
		local hacsubtitle : subinstr local hacsubtitle "heteroskedasticity" "heterosk.", all word
		local hacsubtitle : subinstr local hacsubtitle "autocorrelation" "autocorr.", all word
	}
	if ("`hacsubtitle'"!="") {
		ereturn local title3 = "Statistics robust to `hacsubtitle'"
		
		if ("`e(kernel)'"!="") local notes " `notes' kernel=`e(kernel)'"
		if ("`e(bw)'"!="") local notes " `notes' bw=`e(bw)'"
		if ("`e(dkraay)'"!="") local notes " `notes' dkraay=`e(dkraay)'"
		local notes `notes' // remove initial space
		if ("`notes'"!="") ereturn local title4 = " (`notes')"
		if ("`notes'"!="") {
			if ("`_dta[_TSpanel]'"!="") local tsset panel=`_dta[_TSpanel]'
			if ("`_dta[_TStvar]'"!="") local tsset `tsset' time=`_dta[_TStvar]'
			local tsset `tsset'
			ereturn local title5 = " (`tsset')"
		}
	}
end
