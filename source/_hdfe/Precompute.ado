cap pr drop Precompute
program define Precompute, rclass
	CheckCorrectOrder precompute
	syntax, KEEPvars(varlist) [DEPVAR(varname numeric) EXCLUDESELF] [TSVARS(varlist)] [OVER(varname numeric)]

**** AVGE PART ****
mata: st_local("N_avge", strofreal(avge_num))
if (`N_avge'>0) {
	forv g=1/`N_avge' {
		Assert ("`depvar'"!=""), msg("hdfe.Precompute error: depvar() required")
		mata: st_local("ivars", avge_ivars[`g'])
		mata: st_local("varlabel", avge_varlabel[`g'])
		mata: st_local("target", avge_target[`g'])
		local W __W`g'__

		local note = cond("`excludeself'"=="",""," (excluding obs. at hand)")
		local original_depvar = cond(substr("`depvar'",1,2)=="__", "`: var label `depvar''", "`depvar'")
		Debug, level(2) msg(" - computing AvgE(`original_depvar') wrt (`varlabel')`note'")

		* Syntax: by ... : AverageOthers varname , Generate(name) EXCLUDESELF
		qui AverageOthers `depvar', by(`ivars') gen(`W') `excludeself'
		char `W'[target] `target'
	}

	* Marked obs should have been previously deleted
	tempvar touse
	mark `touse'
	markout `touse' __W*__
	qui keep if `touse'
	drop `touse'
	local keepvars `keepvars' __W*__
}

	Assert c(N)>0, rc(2000) msg("Empty sample, check for missing values or an always-false if statement")
	Assert c(N)>1, rc(2001)

**** ABSORB PART ****
	mata: st_local("G", strofreal(G))
	mata: st_local("N_clustervars", strofreal(length(clustervars)))

	* 1. Clustervars
	* Corner case: if panelvar or timevar are clustervars, we can't touch them, because
	* i) the bandwidth calculation would be wrong if we have time holes, ii) -avar- and -ivreg2- will complain

	local clustervars
	forval i = 1/`N_clustervars' {
		mata: st_local("cluster_ivars", clustervars_ivars[`i'])

		* if clustervar is a panel/time var *AND* we are using a HAC VCE, then we can't touch it
		local is_tsvar : list cluster_ivars in tsvars

		local newname
		forv g=1/`G' {
			mata: fe2local(`g')
			if (`is_mock') continue
			
			local match : list cluster_ivars === ivars
			if (`match' & !`is_tsvar') local newname = "__FE`g'__"
		}

		* If clustervar is an interaction not found in absvars, create identifier variable
		local num_ivars : word count `cluster_ivars'
		if (`num_ivars'>1 & "`newname'"=="") {
			local newname __clustervar`i'__
			GenerateID `cluster_ivars',  gen(`newname')
		}
		if ("`newname'"!="") {
			mata: st_local("oldname", clustervars[`i'])
			mata: clustervars[`i'] = "`newname'"
			Debug, level(3) msg(" - clustervar `oldname' (" as result "`cluster_ivars'" as text ") -> " as result "`newname'")
		}
		else {
			mata: st_local("newname", clustervars[`i'])
		}
		local clustervars `clustervars' `newname'
	}

	* 2. Absvars

	* Get list of all cvars
	forv g=1/`G' {
		mata: fe2local(`g')
		local num_ivars : word count `ivars'
		local all_cvars : list all_cvars | cvars
	}
	local all_cvars : list uniq all_cvars

	* Create IDs for the absvars.
	* Will replace the varname except if i) is interaction so we can't, and ii) it's not interaction but the ivar is the cvar of something else
	* Also, if its in keepvars we can't replace it

	forv g=1/`G' {
		mata: fe2local(`g')
		if (`is_mock') continue
		local num_ivars : word count `ivars'
		local is_cvar : list ivars & all_cvars
		local is_cvar = "`is_cvar'"!=""
		local is_over = "`ivars'"=="`over'"

		local in_keepvars 0
		if (`num_ivars'==1) local in_keepvars : list ivars in keepvars

		if (`num_ivars'>1 | `is_cvar' | `in_keepvars' | `is_over') {
			GenerateID `ivars',  gen(__FE`g'__)
		}
		else {
			GenerateID `ivars' , replace
			rename `ivars' __FE`g'__
		}

		qui su __FE`g'__, mean
		local num_cats = r(max)
		Assert `num_cats'>0
		local name : char __FE`g'__[name]
		Debug, level(3) msg(as text " - absvar`g' " as result "`name'" as text " -> " as result "__FE`g'__")
		Debug, level(1) msg(as text " - absvar`g' " as result "`name'" as text " has " as result "`num_cats'" as text " categories")

	}

	* 3. Epilogue
	
	* Reduce dataset before preparing mata objects (which uses memory)
	keep `keepvars' `weightvar' `clustervars' `all_cvars' __FE*__

	* Fill in auxiliary Mata structures
	Debug, level(2) tic(20)
	mata: prepare()
	Debug, level(2) toc(20) msg("mata:prepare took")
end

