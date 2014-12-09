// -------------------------------------------------------------
// Calculate DoF lost from the FEs
// -------------------------------------------------------------
cap pr drop EstimateDoF
program define EstimateDoF
syntax, dofmethod(string) [clustervar1(string) group(name) uid(varname) groupdta(string)]
	
	Assert inlist("`dofmethod'", "bounds", "simple", "naive", "bootstrap")
	Assert "`dofmethod'"!="bootstrap" , msg("DoF Bootstrap: not yet implemented!")
	Debug, level(1) msg("(calculating degrees of freedom lost by the FEs)")
	Debug, level(2) msg(" - dofmethod: `dofmethod'")
	mata: st_local("G", strofreal(G))

* Start conservatively, assuming M`g'=1 (or =0 if interacted with cont. var)
	local g_list // List of FEs where M is still unknown
	forv g=1/`G' {
		reghdfe_absorb, fe2local(`g')
		// ivars cvars target varname varlabel is_interaction is_cont_interaction is_bivariate is_mock levels
		local M`g' = ("`cvars'"=="") | (`is_bivariate' & !`is_mock')
		* We know M1=1 and Mg=0 if g is has a continuous interaction. The others we don't
		if (`M`g''==1 & `levels'>1) local g_list `g_list' `g'
		// If the FE is just an intercept, it's redundant (useful with -over-)
	}


* (ADDENDUM) Look for nested within cluster for the categorical components of a FE with cont. interaction
	if ("`clustervar1'"!="" & "`dofmethod'"!="naive") {
		forv g=1/`G' {
			reghdfe_absorb, fe2local(`g')
			if (`M`g''==0) {
				local gg = `g' - `is_mock'
				cap _xtreg_chk_cl2 `clustervar1' __FE`gg'__
				assert inlist(_rc, 0, 498)
				if (!_rc) {
					local prettyname : char __FE`gg'__[name]
					Debug, msg("(Absorbed variable " as result "`prettyname'" ///
							as text " is nested within cluster " ///
							as result "`clustervar1'" as text ", adjusting DoF)")
					local M`g' `levels'
				}
			}
		}
	}

* Is the cluster var one of the FEs, or contains one of the FEs?
* Since the number of clusters is the effective "number of observations", 
* we shouldn't penalize the DoF for estimating means within a cluster ("within an obs")
	local g_cluster 0
	if ("`clustervar1'"!="" & "`dofmethod'"!="naive") {
	
* 1) See if it's exactly the same variable
		local regex = regexm("`clustervar1'","^__FE([0-9]+)__$")
		if (`regex') {
			local g_cluster = `=regexs(1)'
			local prettyname : char __FE`g_cluster'__[name]
			Debug, msg("(cluster variable " as result "`prettyname'" ///
				as text " is also an absorbed variable, adjusting DoF)")
		}
* 2) If that failed, see if one of the panels is nested within cluster
* i.e., if whenever two obs are in the same group (in the FE), 
* they are also in the same cluster
		else {
			foreach g of local g_list {
				reghdfe_absorb, fe2local(`g')
				cap _xtreg_chk_cl2 `clustervar1' __FE`g'__
				assert inlist(_rc, 0, 498)

				if (!_rc) {
					local g_cluster = `g'
					local prettyname : char __FE`g'__[name]
					Debug, msg("(Absorbed variable " as result "`prettyname'" ///
							as text " is nested within cluster " ///
							as result "`clustervar1'" as text ", adjusting DoF)")
					continue, break
				}
			}
		}

		if (`g_cluster'>0) {
			reghdfe_absorb, fe2local(`g_cluster')
			local M`g_cluster' = `levels'
			local g_list : list g_list - g_cluster
		}
	} // clustervar1

* Compute connected groups for the remaining FEs (except those with cont interactions)

	if ("`group'"!="") local group_option ", gen(`group')"
	local length : list sizeof g_list
	tokenize `g_list' // Saves g1 g2 .. into `1' `2' etc

	if (`length'==0) {
		// pass
	}
	else if (`length'==1) {
		local M`1' = 1
		local g_list // Only one FE besides cont. interaction and clustervar
	}
	else if (`length'==2) {
		Debug, level(1) msg(" - exact DoF computation using connected groups")
		ConnectedGroups __FE`1'__ __FE`2'__ `group_option'
		local M`1' = 1
		local M`2' = r(groups)
		local g_list
		local label "`: char __FE`1'__[name]' and `: char __FE`2'__[name]'"
	}
	else {
		ConnectedGroups __FE`1'__ __FE`2'__ `group_option'
		local M`1' = 1
		local M`2' = r(groups)
		local g_list : list g_list - 1
		local g_list : list g_list - 2

		* Get a conservative bound for M3 M4, etc:
		* Calculate the mobility group wrt M1 and M2 and use the max of both
		* Of course, this excludes effects in FE3 that are collinear with FE1 and FE2 but not FE1 or FE2 separately (and also excludes searching for FE4 vs FE3, etc)

		if ("`dofmethod'"=="bounds") {
			local msg_done = 0
			* Iterate over the remaining elements
			foreach g of local g_list {
				assert `1'!=`g' & `2'!=`g'
				ConnectedGroups __FE`1'__ __FE`g'__
				local candidate1 = r(groups)
				ConnectedGroups __FE`2'__ __FE`g'__
				local candidate2 = r(groups)
				local M`g' = max(`candidate1', `candidate2')

				reghdfe_absorb, fe2local(`g')
				Assert (`M`g''<=`levels'), msg("Mg should be at most Kg: `M`g''<=`levels'")
				
				* If Mg==Kg , the bound is exact and the FE is redundant!
				if (`M`g''==`levels') local g_list : list g_list - g
			}
			if ("`g_list'"!="") di in ye "(note: conservative DoF estimates; they do not account for all possible collinearities between FEs; dofmethod=`dofmethod')"
		}
		else if inlist("`dofmethod'","simple","naive") {
			di in ye "(note: DoF estimates ignore most possible collinearities between FEs, as dofmethod=`dofmethod')"
		}
		else {
			error 999
		}
	}

	local SumM 0
	local SumK 0
	forv g=1/`G' {
		reghdfe_absorb, fe2local(`g')
		// ivars cvars target varname varlabel is_interaction is_cont_interaction is_bivariate is_mock levels
		assert !missing(`M`g'') & !missing(`levels')
		local SumM = `SumM' + `M`g''
		local SumK = `SumK' + `levels'
		local is_exact = !`: list g in g_list'

		c_local M`g' `M`g''
		c_local K`g' `levels'
		c_local M`g'_exact `is_exact'
		Debug, level(2) msg(" - parameters of FE`g': K=`levels' M=`M`g'' is_exact=`is_exact'")
	}
	local NetSumK = `SumK' - `SumM'
	Debug, level(2) msg(" - DoF loss due to FEs: Sum(Kg)=`SumK', M:Sum(Mg)=`SumM' --> KK:=SumK-SumM=`NetSumK'")

* Save mobility group if needed
	c_local saved_group = 0
	if ("`group'"!="" & `length'>=2) {
		tempfile backup
		qui save "`backup'"
		
		keep `uid' `group'
		sort `uid'
		la var `group' "Mobility group between `label'"
		qui save "`groupdta'" // A tempfile from the caller program
		Debug, level(2) msg(" - mobility group saved")
		qui use "`backup'", clear
		cap erase "`backup'"
		c_local saved_group = 1
	}

	c_local M `SumM'
	c_local kk `NetSumK'
	local fe_nested_in_cluster = (`g_cluster'>0)
	c_local fe_nested_in_cluster `fe_nested_in_cluster'
end
