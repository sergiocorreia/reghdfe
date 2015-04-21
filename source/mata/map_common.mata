mata:
mata set matastrict on

real rowvector safe_divide(real rowvector numerator, real rowvector denominator, | real scalar epsi) {
	 // If the denominator goes below machine precision, the division explodes
	 if (args()<3) epsi = epsilon(1)
	return( numerator :/ colmax(denominator \ J(1,cols(denominator),epsi)) )
}

// -------------------------------------------------------------------------------------------------

void verbose2local(`Problem' S, string scalar loc) {
	st_local(loc, strofreal(S.verbose))
}

// -------------------------------------------------------------------------------------------------

void function store_uid(`Problem' S, `Varname' varname) {
	S.uid = st_data(., varname)
	assert_msg(rows(S.uid)==S.N, "assertion failed: rows(S.uid)==S.N")
}

void function drop_uid(`Problem' S) {
	S.uid = J(0,0,.)
}

// -------------------------------------------------------------------------------------------------

void function store_resid(`Problem' S, `Varname' varname) {
	S.resid = st_data(., varname)
	S.residname = varname
	assert_msg(rows(S.resid)==S.N, "assertion failed: rows(S.resid)==S.N")
}

void function resid2dta(`Problem' S) {
	st_store(., st_addvar("double", S.residname), S.resid)
	S.resid = J(0,0,.)
	S.residname = ""
}

// -------------------------------------------------------------------------------------------------

void function groupvar2dta(`Problem' S) {
	if (S.groupvar!="") {
		if (S.verbose>2) printf("{txt}    - Saving identifier for the first mobility group: {res}%s\n", S.groupvar)
		// WRONG: NEED TO INDEX BY UID
		st_store(S.uid, st_addvar(S.grouptype, S.groupvar), S.groupseries)

		S.groupseries = J(0,0,0)
		st_varlabel(S.groupvar, S.grouplabel)
	}
	else {
		if (S.verbose>=0) printf("{txt}Note: mobility group not saved as there was no need to compute it (e.g. there are less than two absvars not nested within a cluster)\n")
	}
}

// -------------------------------------------------------------------------------------------------

void function drop_ids(`Problem' S) {
	`Integer' g
	for (g=1;g<=S.G;g++) {
		if (!S.fes[g].is_clustervar & S.fes[g].target=="") st_dropvar(S.fes[g].idvarname)
	}
}

// -------------------------------------------------------------------------------------------------

void function esample2dta(`Problem' S, `Varname' esample) {
	assert(length(S.uid)>0)
	st_store(S.uid, st_addvar("byte", esample), J(rows(S.uid),1,1) )
}

// -------------------------------------------------------------------------------------------------

// Copy the fixed effect estimates (the alphas back into the original dataset)
void function alphas2dta(`Problem' S) {
	`Integer' g, i
	`Varlist' target
	`String' varlabel
	assert(S.will_save_fe==1)
	if (S.verbose>1) printf("{txt}    - Storing fixed effects in the original dataset\n")
	for (g=1; g<=S.G; g++) {
		target = S.fes[g].target
		if (length(target)>0) {
			st_store(S.uid, st_addvar("double", target), S.fes[g].alphas)
			S.fes[g].alphas = J(0,0,0)
			for (i=1; i<=length(target);i++) {
				varlabel = invtokens(S.fes[g].ivars, "#")
				if (i>1 | !S.fes[g].has_intercept) varlabel = varlabel + "#c." + S.fes[g].cvars[i-S.fes[g].has_intercept]
				st_varlabel(target[i], sprintf("[FE] %s", varlabel))
			}
		}
	}
}

end
