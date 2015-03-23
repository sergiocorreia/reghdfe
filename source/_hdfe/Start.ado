cap pr drop Start
program define Start, rclass
	CheckCorrectOrder start
	syntax, Absorb(string) [AVGE(string)] [CLUSTERVARS(string)] [OVER(varname numeric)] [WEIGHT(string) WEIGHTVAR(varname numeric)]
	Assert !regexm("`absorb'","[*?-]"), ///
		msg("error: please avoid pattern matching in -absorb-")

	if ("`over'"!="") Assert "`avge'"=="", msg("-avge- needs to be empty if -over- is used")

	Assert inlist("`weight'", "", "fweight", "aweight", "pweight")

**** ABSORB PART ****

* First pass to get the true number of FEs
	local i 0
	Debug, level(3) msg(_n "Fixed effects:")
	foreach var of local absorb {
		ParseOneAbsvar, absvar(`var')
		local i = `i' + cond(r(is_bivariate), 2, 1)
		* Output: r(target) cvars ivars is_interaction is_cont_interaction is_bivariate
		Assert `i'>1 | "`r(cvars)'"=="" | `r(is_bivariate)', ///
			msg("error parsing absorb : first absvar cannot be continuous interaction" ///
			_n "solution: i) reorder absvars, ii) replace # with ##, iii) add a constant as first absvar (as a workaround)")

		if ("`over'"!="") {
			local ivars r(ivars)
			local dupe : list ivars & over
			Assert ("`dupe'"==""), msg("-over- cannot be part of any absvar")
		}
	}

	if ("`over'"!="") {
		local ++i // We'll add -over- as the first FE
		local pre_absorb `absorb'
		local absorb `over' `absorb'
	}

* Create vector of structures with the FEs
	Assert inrange(`i',1,100), msg("error: too many absorbed variables (do not include the dummies, just the variables)")
	Debug, msg(`"(`i' absorbed fixed `=plural(`i',"effect")': "' as result "`absorb'" as text ")")
	mata: weightexp = ""
	mata: weightvar = ""
	if ("`weightvar'"!="") {
		Debug, msg(`"(`weight': "' as result "`weightvar'" as text ")")
		mata: weightexp = "[`weight'=`weightvar']"
		mata: weightvar = "`weightvar'"
		**qui cou if `fweight'<=0 | `fweight'>=. | (`fweight'!=int(`fweight'))
		** Move this somewhere else.. else it will fail needlesly if some excluded obs. have missing weights
		**Assert (`r(N)'==0), msg("fweight -`fweight'- can only have strictly positive integers (no zero, negative, MVs, or reals)!")
	}
	mata: G = `i'
	mata: initialize()

* Second pass to save the values
	local i 0
	foreach var of local absorb {
		qui ParseOneAbsvar, absvar(`over_prefix'`var')
		local keepvars `keepvars' `r(ivars)' `r(cvars)'
		local varlabel = "i." + subinstr("`r(ivars)'", " ", "#i.", .)
		if (`r(is_cont_interaction)' & !`r(is_bivariate)') local varlabel "`varlabel'#c.`r(cvars)'"
		
		local args `" "`r(target)'", "`r(ivars)'", "`r(cvars)'", `r(is_interaction)', `r(is_cont_interaction)', `r(is_bivariate)', "`weightvar'" "'
		mata: add_fe(`++i', "`varlabel'", `args', 0)
		if (`r(is_bivariate)') {
			local varlabel "`varlabel'#c.`r(cvars)'"
			mata: add_fe(`++i', "`varlabel'", `args', 1)
		}

		if ("`over'"!="") local over_prefix "i.`over'#" // Not for the first one
	}
	local N_hdfe = `i'

	if ("`over'"!="") Debug, msg(`"absvars expanded due to over: `pre_absorb' -> `absorb'"')

**** AVGE PART ****

* First pass to get the true number of FEs
local N_avge = 0
if ("`avge'"!="") {
	local i 0
	foreach var of local avge {
		Debug, level(3) msg(_n "AvgE effects:")
		ParseOneAbsvar, absvar(`var')
		local ++i
		* Output: r(target) cvars ivars is_interaction is_bivariate
		Assert ("`r(cvars)'"=="" & `r(is_bivariate)'==0), ///
			msg("error parsing avge : continuous interactions not allowed")
	}

* Create vectors
	Assert inrange(`i',1,100), msg("error: too many avge variables (do not include the dummies, just the variables)")
	Debug, msg(`"(`i' avge `=plural(`i',"effect")': "' as result "`avge'" as text ")")
}

* Always save this to avoid not-found errors
	mata: avge_ivars = J(1, `i', "")
	mata: avge_target = J(1, `i', "")
	mata: avge_varlabel = J(1, `i', "")

* Second pass to save the values
if ("`avge'"!="") {
	local i 0
	foreach var of local avge {
		qui ParseOneAbsvar, absvar(`var')
		local ++i
		local varlabel = "i." + subinstr("`r(ivars)'", " ", "#i.", .)
		mata: avge_ivars[`i'] = "`r(ivars)'"
		mata: avge_target[`i'] = "`r(target)'"
		mata: avge_varlabel[`i'] = "`varlabel'"
		local keepvars `keepvars' `r(ivars)'
	}
	local N_avge = `i'
}
	mata: avge_num = `N_avge'

*** CLUSTER PART ****
* Create two string rowvectors, with the variables and ivars, and also add the ivars to keepvars
* EG: If clustervar1=foreign, absorb=foreign, then clustervar1 -> __FE1__
	mata: clustervars = tokens("`clustervars'")
	mata: clustervars_ivars = J(1, length(clustervars), "")
	mata: clustervars_original = J(1, length(clustervars), "")

	local i 0
	foreach var of local clustervars {
		local ++i
		Debug, level(3) msg(_n "Cluster by:")
		ParseOneAbsvar, absvar(`var')
		Assert "`r(cvars)'"=="", msg("clustervar cannot contain continuous interactions")
		local ivars = r(ivars)
		mata: clustervars_ivars[`i'] = "`ivars'"
		mata: clustervars_original[`i'] = invtokens( tokens(clustervars_ivars[`i']) , "#")
		local keepvars `keepvars' `ivars'
	}

**** Returns ****
	Debug, level(3) newline
	local keepvars : list uniq keepvars
	return local keepvars `keepvars'
	return scalar N_hdfe = `N_hdfe'
	return scalar N_avge = `N_avge'
end
