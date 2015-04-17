// -------------------------------------------------------------
// MATA functions and objects
// -------------------------------------------------------------
// struct FE: 				Container with data for each FE
// prepare(fe_varlist): 	Constructs the basic objects in order to obtain resid later on
// make_residual(varname): 	Obtains residual of varname wrt to the previously indicated FEs
// transform(x, g0, g1):	Transforms x from series (g=0) or means by FE (g>1) into either of those
// count_by_group(): 		Called by prepare

// -------------------------------------------------------------
// Shortcuts
// -------------------------------------------------------------
local Varlist 		string scalar
local Integer 		real scalar
local VarByFE 		real colvector // Should be levels*1
local Series		real colvector // Should be N*1
local Matrix		real matrix
local SharedData 	external struct FixedEffect vector

mata:
mata set matastrict on

// -------------------------------------------------------------
// Structures
// -------------------------------------------------------------
	// Notice that every FE costs around 3 -doubles- of memory (obs*8 bytes)
	struct FixedEffect {
		`Integer' levels // Distinct values of FE (recall its sorted 1..levels)
		`Integer' g // Identifier of this FE
		`Integer' is_interaction, is_cont_interaction, is_bivariate, is_mock
		`Integer' K // Number of params in a multivariate by-group regression (usually 1 or 2)
		`Matrix' v // The continuous var, if is_cont_interaction==1
		`Matrix' invxx // The inv(X'X) matrix where X is v + constant (constant last!)
		`VarByFE' count, sum_count // count=levels*1 with a count of each level. sum_count=Running sum of the above!
		`Series' group, indexfrom0, indexfrom1, sorted_weight // weight is freq. weight, sorted by indexfrom0
		`Varlist' varname, Z, target, ivars, cvars, varlabel, weightvar // weightvar is just the varname with the frequencies
	}

// -------------------------------------------------------------
// PREPARE DATA: Create almost-empty aux structures and fill them
// -------------------------------------------------------------

// Create the structures and shared variables
void function initialize() {
	`SharedData' FEs
	external `Integer' G // Number of FEs
	external `Matrix' betas

	assert(G>0 & G<=100)
	FEs = J(G, 1, FixedEffect()) // Use G=100 for now (should be enough)
	betas = 0
}

// Add basic data into each structure
void function add_fe(
		`Integer' g, `Varlist' varlabel, `Varlist' target, 
		`Varlist' ivars, `Varlist' cvars, `Integer' is_interaction,
		`Integer' is_cont_interaction, `Integer' is_bivariate,
		`Varlist' weightvar,
		`Integer' is_mock) {
	`SharedData' FEs
	assert(is_mock==0 | is_mock==1)
	FEs[g].g = g
	FEs[g].ivars = ivars
	FEs[g].cvars = cvars
	FEs[g].is_interaction = is_interaction
	FEs[g].is_bivariate = is_bivariate
	FEs[g].is_mock = is_mock
	FEs[g].K = 1 + length(tokens(cvars)) // 1 w/out cvars, 2 with bivariate, +2 with multivariate within the group
	FEs[g].is_cont_interaction = is_cont_interaction
	FEs[g].varname = "__FE" + strofreal(g-is_mock) + "__"
	FEs[g].Z = "__Z" + strofreal(g) + "__"
	FEs[g].varlabel = varlabel
	FEs[g].levels = -1 // Not yet filled
	FEs[g].target = target + (is_mock & target!="" ? "_slope" : "")
	if (FEs[g].K>2 & FEs[g].target!="") FEs[g].target = FEs[g].target + strofreal(FEs[g].K-1)
	FEs[g].weightvar = weightvar
}

// Dump data of one FE to locals
void function fe2local(`Integer' g) {
	`SharedData' FEs
	stata("local ivars " + FEs[g].ivars)
	stata("local cvars " + FEs[g].cvars)
	stata("local target " + FEs[g].target)
	stata("local varname " + FEs[g].varname)
	stata("local Z " + FEs[g].Z)
	stata("local varlabel " + FEs[g].varlabel)
	stata("local is_interaction " + strofreal(FEs[g].is_interaction))
	stata("local is_cont_interaction " + strofreal(FEs[g].is_cont_interaction))
	stata("local is_bivariate " + strofreal(FEs[g].is_bivariate))
	stata("local is_mock " + strofreal(FEs[g].is_mock))
	stata("local levels " + strofreal(FEs[g].levels))
	stata("local group_k " + strofreal(FEs[g].K))
	stata("local weightvar " + FEs[g].weightvar)
}

// Fill aux structures
void function prepare() {

	external `Integer'	G
	`SharedData' 		FEs
	external `Integer' 	VERBOSE
	`Integer'			g, obs, mem_used, is_weighted
	`Series'			group, weight
	`VarByFE'			count
	`Varlist'			weightvar

	// Setup
	assert(VERBOSE>=0 & VERBOSE<=5)
	obs = st_nobs()
	weightvar = ""
	is_weighted = (FEs[1].weightvar!="")
	if (is_weighted) {
		weightvar = FEs[1].weightvar
		weight = st_data(., weightvar)
	}
	else {
		weight = 1
	}

	// Main code
	for (g=1;g<=G;g++) {

		if (FEs[g].is_mock) {
			FEs[g].levels = FEs[g-1].levels
			continue
		}

		if (VERBOSE>1) printf("{txt}(preparing matrices for fixed effect {res}%s{txt})\n", FEs[g].varlabel)
		FEs[g].group = st_data(., FEs[g].varname)
		if (max(FEs[g].group)>1e8) _error("More than 100MM FEs found. Are you sure FE is in format 1..G?")
		if (min(FEs[g].group)!=1) _error("Minimum value for FE is not 1")
		FEs[g].indexfrom0 = order(FEs[g].group,1)
		if (g>1) FEs[g].indexfrom1 = FEs[1].group[FEs[g].indexfrom0]
		FEs[g].sorted_weight = is_weighted ? weight[FEs[g].indexfrom0, 1] : 0
		count = count_by_group(FEs[g].group, FEs[g].indexfrom0)
		FEs[g].sum_count = quadrunningsum(count)

		// If we have cont. interactions, use the above -count- only for -sum_count-
		// and use this one to get the denominators used in transform
		// Why? because -sum_count- is used to offset the submatrices, while -count- will be a denominator

		if (FEs[g].is_cont_interaction) {
			FEs[g].v = st_data(., FEs[g].cvars)

			if (!FEs[g].is_bivariate) {
				count = count_by_group(FEs[g].group, FEs[g].indexfrom0, FEs[g].v, FEs[g].sorted_weight)
			}
			else {
				count = count_by_group(FEs[g].group, FEs[g].indexfrom0, 0, FEs[g].sorted_weight)
				FEs[g].invxx = compute_invxx(FEs[g].v ,FEs[g].indexfrom0, FEs[g].sum_count, count, FEs[g].sorted_weight)
			}
		}
		else if (is_weighted) {
			count = count_by_group(FEs[g].group, FEs[g].indexfrom0, 0, FEs[g].sorted_weight)
		}
		
		FEs[g].count = count
		FEs[g].levels = rows(FEs[g].sum_count)

		if (VERBOSE>1) {
			mem_used = ( sizeof(FEs[g].group) + sizeof(FEs[g].indexfrom0) + sizeof(FEs[g].indexfrom1) + sizeof(FEs[g].sum_count) ) / 2^20
			printf(" - %3.1fMB used\n", mem_used)
		}
	}
}

// -------------------------------------------------------------
// TRANSFORM: Transform the dim of a vector by taking avgs
// -------------------------------------------------------------
// Type <`VarByFE'> but also works with <`Series'> (<real colvector> both)
// NOTE: g==0 means the raw Stata series, g==1 is collapsed by FE1, etc.
// (This function could benefit from a refactoring...)
//

`VarByFE' function transform(`VarByFE' indata, `Integer' g_from, `Integer' g_to) {
	`VarByFE'		outdata
	`SharedData'	FEs
	`Integer'		from_v, to_v // Whether the new or old FEs have cont. interactions
	`Integer'		is_weighted
	`Series'		indata_v
	`Integer'		from_bivariate, to_bivariate
	external `Matrix' betas

	assert(g_from>=0 & g_from<=100)
	assert(g_to>=0 & g_to<=100)
	assert(g_to!=g_from)
	is_weighted = (FEs[1].weightvar!="")

	// -v- is the possible cont. interaction of the FE
	to_v = from_v = 0 // There is no FEs[0]
	if (g_to>0) to_v = FEs[g_to].is_cont_interaction
	if (g_from>0) from_v = FEs[g_from].is_cont_interaction // & (!FEs[g_from].is_bivariate | FEs[g_from].is_mock)
	
	// g_from, g_to, from_v, to_v

	if (g_to>0) {
		if (FEs[g_to].is_mock==1) {
			_error("g_to shouldn't be the 2nd part of a ## interaction")
		}
	}

	if (g_from>0) {
		if (FEs[g_from].is_mock==1) {
			--g_from // Go back to non-mock
		}
	}


	if (g_to>0) {
		// If g_to is the first FE of a bivariate, do the regression and save the cache
		assert(FEs[g_to].is_bivariate==0 | FEs[g_to].is_bivariate==1)
		if (FEs[g_to].is_bivariate==1) {
			assert(!FEs[g_to].is_mock)
			// We have to deal with 0->g and 1->g transforms

			if (g_from>0 & !from_v) {
				assert(g_from==1)
				indata = indata[FEs[g_from].group, 1]
			}
		
			// -outdata- will have the predicted values (i.e. yhat = alpha + v1*beta1 + ..) by group, expanded
			outdata= regress_by_group(indata, FEs[g_to].v, FEs[g_to].indexfrom0, FEs[g_to].sum_count, FEs[g_to].count, FEs[g_to].invxx, FEs[g_to].sorted_weight, FEs[g_to].group)
			assert(rows(outdata)==st_nobs())
			assert(cols(outdata)==1)
			// b1 b2 .. alpha are in the external matrix -betas- , of size (K*levels,K)
			assert(rows(betas)==FEs[g_to].levels)
			assert(cols(betas)==FEs[g_to].K)
			return(outdata)
		}
	}

	if (!to_v & !from_v) {
		// sorted_weight should be of the same group as sum_count .. 
		if (g_to==0) {
			outdata = indata[FEs[g_from].group,1]
		}
		else if (g_from==0) {
			outdata = mean_by_group(indata, FEs[g_to].indexfrom0, FEs[g_to].sum_count, FEs[g_to].count, FEs[g_to].sorted_weight)
		}
		else if (g_from==1) {
			outdata = remean_by_group(indata, FEs[g_to].indexfrom1, FEs[g_to].sum_count, FEs[g_to].count, FEs[g_to].sorted_weight)
		}
		else {
			_error(sprintf("\nNot implemented"))
		}
	}
	else { // With cont interaction but not bivariate
		if (from_v & g_to==0) return(indata)
		
		if (from_v & g_to>0 & !to_v) {
			outdata = mean_by_group(indata, FEs[g_to].indexfrom0, FEs[g_to].sum_count, FEs[g_to].count, FEs[g_to].sorted_weight)
			return(outdata)
		}

		// In the remaining cases, the output has same dim as -v-
		
		if (g_from==0 | (from_v & to_v) ) {
			// Inline it?? BUGBUG TODO
			indata_v = FEs[g_to].v :* indata
		}
		else if (g_from==1 & !from_v) { // Either to_v==1 or g_to==0
			// Convert to a g==0 dimension
			indata_v = FEs[g_to].v :* indata[FEs[g_from].group, 1]
		}
		else {
			_error(sprintf("\nNot implemented"))
		}
		// Will have MVs for groups where all obs of V are zero
		outdata = editmissing(mean_by_group(indata_v, FEs[g_to].indexfrom0, FEs[g_to].sum_count, FEs[g_to].count, FEs[g_to].sorted_weight) , 0)
		outdata = FEs[g_to].v :* outdata[FEs[g_to].group, 1]
		assert(rows(outdata)==st_nobs())
	}
	return(outdata)
}


// -------------------------------------------------------------
// COUNT_BY_GROUP: alternative to mm_freq(group) (~10% of runtime)
// -------------------------------------------------------------
`VarByFE' function count_by_group(`Series' group, `Series' index, | `Series' v, `Series' sorted_weight)
{
	`Integer' 	levels, obs, i, count, g, is_cont_interaction, is_weighted, ww, vv
	`Series'	sorted_group, sorted_v
	`VarByFE'	ans
	sorted_group = group[index, 1]
	levels = max(group)
	obs = rows(group)
	ans = J(levels, 1, 0)
	if (levels>obs) _error("Error: more levels of FE than observations!")
	is_cont_interaction = (args()>=3 & length(v)>1 ) // -v- is the continuous interaction
	is_weighted = (args()>=4 & length(sorted_weight)>1 )
	if (is_cont_interaction) sorted_v = v[index]
	
	// -g- iterates over the values of the FE, -i- over observations
	// Should we have the -if- for v and weights inside or outside the hot loop?
	// If we put the -if- outside, we can use "count++" for the simple case
	// Check how large is the slowdown...
	count = 0
	for (i=g=1; i<=obs; i++) {
		if (g<sorted_group[i]) {
			ans[g++] = count
			count = 0
		}
		ww = is_weighted ? sorted_weight[i] : 1
		vv = is_cont_interaction ? sorted_v[i] ^ 2 : 1
		count = count + ww * vv
	}
	ans[g] = count // Last group
	if (!is_cont_interaction) assert( all(ans) ) // assert( all(ans:>0) )
	// BUGBUG -ans- may be zero for some cases!!!!!!!!!!!
	return(ans)
}

// -------------------------------------------------------------
// MEAN_BY_GROUP: Take a N*1 vector and save the avg by FE into a levels*1 vector
// -------------------------------------------------------------
`VarByFE' function mean_by_group(`Series' indata, `Series' index, `VarByFE' sum_count, `VarByFE' counts_to, `Series' sorted_weight)
{
	`Integer'	levels, i, j_lower, j_upper
	`Series'	sorted_indata
	`VarByFE'	outdata

	assert(rows(indata)==rows(index))
	levels = rows(sum_count)
	sorted_indata = indata[index] // Also very crucial to speed
	if (length(sorted_weight)>1) sorted_indata = sorted_indata :* sorted_weight
	outdata = J(levels, 1 , 0)
	
	// !! This is one of the most hot / crucial loops in the entire program
	j_lower = 1
	for (i=1; i<=levels; i++) {
		j_upper = sum_count[i]
		outdata[i] = quadcolsum(sorted_indata[| j_lower \ j_upper |])
		j_lower = j_upper + 1
	}
	outdata = outdata :/ counts_to
	return(outdata)
}

// -------------------------------------------------------------
// REMEAN_BY_GROUP: Transform one mean by group into another (of a diff group)
// -------------------------------------------------------------
`VarByFE' function remean_by_group(`VarByFE' indata, `Series' index, `VarByFE' sum_count, `VarByFE' counts_to, `Series' sorted_weight)
{
	`Integer'	levels_from, levels_to, i, j_lower, j_upper, obs
	`Series'	se_indata
	`VarByFE'	outdata

	obs = rows(index)
	levels_to = rows(sum_count)
	levels_from = rows(indata)
	assert(obs==st_nobs())
	assert(levels_from==max(index))
	//assert( all(indata:<.) )
	se_indata = indata[index, 1] // SE = Sorted & Expanded
	if (length(sorted_weight)>1) se_indata = se_indata :* sorted_weight
	outdata = J(levels_to, 1 , 0)
	// assert( all(se_indata:<.) )
	
	// !! This is one of the most hot / crucial loops in the entire program
	j_lower = 1
	for (i=1; i<=levels_to; i++) {
		j_upper = sum_count[i]
		outdata[i] = quadcolsum(se_indata[| j_lower \ j_upper |])
		j_lower = j_upper + 1
	}

	outdata = outdata :/ counts_to
	// mean() is much* slower than doing -quadcolsum- and then dividing by counts_to
	return(outdata)
}

// -------------------------------------------------------------
// REGRESS_BY_GROUP: Multivariate regression on constant and at least 1 var.
// -------------------------------------------------------------
// Returns block-column matrix with estimates by group; last estimate is the constant
// (This function hasn't been optimized very much)
`Matrix' function regress_by_group(`Series' y, `Matrix' x, `Series' index, 
	`VarByFE' offset, `VarByFE' count, `Matrix' invxx, `Series' sorted_weight, group)
{
	`Integer'			N, K, levels, is_weighted, j_lower, j_upper, i
	`Series'			predicted, tmp_y, tmp_w, sorted_y
	real colvector		b
	external `Matrix'   betas
	`Matrix'			tmp_x, tmp_invxx, sorted_x

	N = rows(x)
	K = 1 + cols(x)
	levels = rows(offset)
	is_weighted = length(sorted_weight)>1
	sorted_y = y[index,.]
	sorted_x = x[index,.]
	predicted = J(N, 1 , 0)
	betas = J(levels, K, 0)
	
	assert(rows(y)==N)
	assert(rows(index)==N)
	assert(rows(count)==levels)
	assert(rows(invxx)==levels*K & cols(invxx)==K)
	if (is_weighted) assert(rows(sorted_weight)==N)
	if (!is_weighted) assert(sorted_weight==0)
	
	j_lower = 1
	for (i=1; i<=levels; i++) {
		j_upper = offset[i]
		tmp_x = sorted_x[| j_lower , 1 \ j_upper , . |]
		tmp_y = sorted_y[| j_lower , 1 \ j_upper , . |]
		tmp_invxx = invxx[| 1+(i-1)*K , 1 \ i*K , . |]
		if (is_weighted) {
			tmp_w = sorted_weight[| j_lower , 1 \ j_upper , 1 |]
			b = tmp_invxx * quadcross(tmp_x, 1, tmp_w, tmp_y, 0)
		}
		else {
			b = tmp_invxx * quadcross(tmp_x, 1, tmp_y, 0)
		}
		betas[i, .] = b'
		//predicted[| j_lower , 1 \ j_upper , . |] = b[K] :+ tmp_x * b[|1 \ K-1|] // Doesn't work b/c its sorted by index, and I didn't save the reverse sort....
		j_lower = j_upper + 1
	}

	predicted = rowsum( (x , J(rows(predicted),1,1)) :* betas[group,.] )
	return(predicted)
}

// -------------------------------------------------------------
// COMPUTE_INVXX
// -------------------------------------------------------------
`Matrix' function compute_invxx(`Matrix' x, `Series' index, `VarByFE' offset, `VarByFE' count, `Series' sorted_weight)
{
	`Integer'	N, levels, K, is_weighted, j_lower, j_upper, i
	`Matrix'	ans, invxx, tmp_x, sorted_x
	`Series'	tmp_w

	N = rows(x)
	K = 1 + cols(x)
	levels = rows(offset)
	is_weighted = length(sorted_weight)>1
	sorted_x = x[index,.]
	ans = J(levels * K, K, 0)
	
	assert(rows(index)==N)
	assert(rows(count)==levels)
	if (is_weighted) assert(rows(sorted_weight)==N)
	if (!is_weighted) assert(sorted_weight==0)
	
	j_lower = 1
	for (i=1; i<=levels; i++) {
		j_upper = offset[i]
		tmp_x = sorted_x[| j_lower , 1 \ j_upper , . |]
		if (is_weighted) {
			tmp_w = sorted_weight[| j_lower , 1 \ j_upper , 1 |]
			invxx = invsym(quadcross(tmp_x,1,tmp_w,tmp_x,1))
		}
		else {
			invxx = invsym(quadcross(tmp_x,1,tmp_x,1))
		}
		ans[| 1+(i-1)*K , 1 \ i*K , . |] = invxx
		j_lower = j_upper + 1
	}
	return(ans)
}

// -------------------------------------------------------------
// MAKE_RESIDUAL: Take a variable and obtain its residual wrt many FEs
// -------------------------------------------------------------
// num_fe: Allows running nested models
void function make_residual(
	`Varlist' varname, `Varlist' resid_varname, 
	`Integer' tolerance, `Integer' max_iterations, | `Integer' save_fe, 
	`Integer' accelerate, `Integer' num_fe,
	`Integer' bad_loop_threshold, `Integer' stuck_threshold,
	`Integer' pause_length, `Integer' accel_freq, `Integer' accel_start)
{
	`SharedData' 		FEs
	external `Integer' 	VERBOSE
	external `Integer'	G
	
	// bad_loop_threshold, stuck_threshold, accel_freq, accel_start, pause_length
	`Integer'	update_error, converged, iter, accelerate_candidate, accelerated, mu, accelerate_norm
	`Integer'	eps, g, obs, stdev, levels, gstart, gextra, k // , _
	`Integer'	acceleration_countdown, old_error, oldest_error, bad_loop, improvement
	`Series' 	y, resid, ZZZ // ZZZ = sum of Zs except Z1
	`VarByFE'	P1y
	string scalar		code, msg
	pointer(`VarByFE') colvector	Deltas, oldDeltas, Zs, oldZs, Pytildes
	external `Matrix'	betas
	
	// Parse options
	assert(G==rows(FEs))
	obs = st_nobs()
	assert(VERBOSE>=0 & VERBOSE<=5)

	if (args()<5 | save_fe==.) save_fe = 0
	if (args()<6 | accelerate==.) accelerate = 1
	if (args()<7 | num_fe==-1) num_fe = G
	if (save_fe!=0 & save_fe!=1) _error("Option -save_fe- must be either 0 or 1")

	// See below for explanation
	if (args()<8 | bad_loop_threshold==-1) bad_loop_threshold = 1
	if (args()<9 | stuck_threshold==-1) stuck_threshold = 5e-3
	if (args()<10 | pause_length==-1) pause_length = 20
	if (args()<11 | accel_freq==-1) accel_freq = 3
	if (args()<12 | accel_start==-1) accel_start = 6
	// BUGBUG: These defaults are in triplicate: here, in Demean (?), and in reghdfe.parse

	// Should I expose these parameters?
	// bad_loop_threshold = 1 // If acceleration seems stuck X times in a row, pause it
	// stuck_threshold = 5e-3 // Call the improvement "slow" when it's less than e.g. 1%
	// pause_length = 20 // This is in terms of candidate accelerations, not iterations (i.e. x3)?
	// accel_freq = 3
	// accel_start = 6

	// Initialize vectors of pointers and others
	gstart = 1 + FEs[1].K // Usually 2, or if bivariate, 3
	Deltas = oldDeltas = Zs = oldZs = Pytildes = J(num_fe,1,NULL)
	ZZZ = J(obs, 1, 0) // oldZZZ = 
	for (g=gstart;g<=num_fe;g++) {
		if (FEs[g].is_cont_interaction) levels = obs
		else levels = FEs[g].levels

		if (FEs[g].is_mock) levels = 0 // Better than not initializing them
		
		Deltas[g] = &J(levels,1,.)
		oldDeltas[g] = &J(levels,1,.)
		Zs[g] = &J(levels,1,0) // Needs to start with 0s
		oldZs[g] = &J(levels,1,0) // Needs to start with 0s
		Pytildes[g] = &J(levels,1,.)
	}
	
	// Calculate P1*y and save in mata, then M1*y==ytilde and save in stata
	if (VERBOSE>0) {
		if (substr(varname, 1, 2)=="__") {
			msg = st_global(varname+"[fvrevar]")
			msg = msg + st_global(varname+"[tsrevar]")
			msg = msg + st_global(varname+"[avge]")
			if (msg=="") msg = "[residuals]"
		}
		else {
			msg = varname
		}
		printf("{txt}(computing residual of {res}%s{txt} with respect to %1.0f FE%s" + (VERBOSE==1 & num_fe>1? " " : ")\n"), msg, num_fe, num_fe==1? "" : "s")
		displayflush()
	}
	if (VERBOSE>1) printf("{txt} - Demeaning wrt FE1\n")

	st_view(y=., ., varname)
	stdev = sqrt(quadvariance(y))
	if (VERBOSE>1) printf("{txt} - Stdev of var is %12.6g\n",stdev)
	if (stdev<1e-8) stdev = 1 // Probably a constant, can't standardize
	
	if (FEs[1].is_cont_interaction & !FEs[1].is_bivariate) {
		_error("error: the first absvar cannot be an interaction with a cont. var")
	}

	P1y = transform(y, 0, 1)
	st_store(., st_addvar("double", resid_varname), transform(P1y, 1, 0)) // Store P1*y
	stata(sprintf(`"qui replace %s = %s - %s"', resid_varname, varname, resid_varname)) // ytilde===M1*y = y - P1*y
	stata(sprintf(`"la var %s "[reghdfe residuals of %s]" "', resid_varname, varname)) // Useful to know what is what when debugging

	if (num_fe<gstart) {
		// Save FEs
		if (save_fe!=0) {
			if (VERBOSE>1) printf("{txt} - Saving FE\n")
			if (gstart==2) {
				st_store(., st_addvar("double", FEs[1].Z), transform(P1y, 1, 0))
			}
			else {
				st_store(., st_addvar("double", FEs[1].Z), betas[FEs[1].group,FEs[1].K] )
				for (k=2; k<=FEs[1].K; k++) {
					st_store(., st_addvar("double", FEs[k].Z), FEs[1].v[.,k-1] :* betas[FEs[1].group,k-1] )
				}
			}
		}

		return
	}

	// Compute P2*ytilde, P3*ytilde and so on
	st_view(resid=., ., resid_varname) // BUGBUG? is using view too slow??
	for (g=gstart;g<=num_fe;g++) {
		if (FEs[g].is_mock) continue
		(*Pytildes[g]) = transform(resid, 0, g) :/ stdev // Standarize to get convergence independent of the scale (1000s, units, etc)
		assert(rows(*Pytildes[g])>0)
	}

	// --------------------------------
	if (VERBOSE>1) printf("{txt} - Starting iteration...\n")
	// --------------------------------
	converged = 0
	eps = epsilon(1) // Which one works better? sqrt(epsilon(1)) // 1e-8  ...  epsilon(1) ~ 2e-16
	old_error = oldest_error = bad_loop = acceleration_countdown = 0
	gextra = gstart + FEs[gstart].is_bivariate
	if (VERBOSE>1) timer_clear(40)
	for (iter=1; iter<=max_iterations; iter++) {
		// _ = _stata("parallel break")

		// Acceleration setup
		
		accelerated = 0
		accelerate_candidate = accelerate & (mod(iter, accel_freq)==1) & (iter>accel_start) // Optimal accel interval? // ==0 ?
		// If the FP is stuck, stop accelerating for a few periods
		code = accelerate_candidate ? "x" : "."
		if (accelerate_candidate==1 & acceleration_countdown>0) {
			--acceleration_countdown
			accelerate_candidate = 0
		}
		
		// Update Zs
		if (VERBOSE>1) timer_on(40)		
		for (g=gstart;g<=num_fe;g++) {
			if (FEs[g].is_mock) continue
			if (accelerate_candidate) (*oldDeltas[g]) = (*Deltas[g]) // Only update when needed

			// -reghdfe.ado- will spend most of its time in this line:
			if (FEs[g].is_bivariate) {
				(*Deltas[g]) = (*Pytildes[g]) + transform(transform(ZZZ,0,1), 1, g) - (num_fe>gextra? transform(ZZZ, 0, g) : (*Zs[g]) )
			}
			else {
				(*Deltas[g]) = (*Pytildes[g]) + transform(transform(ZZZ,0,1), 1, g) - (num_fe>gextra? transform(ZZZ, 0, g) : (*Zs[g]) )
			}


			(*Zs[g]) = (*Zs[g]) + (*Deltas[g])
			ZZZ = ZZZ + transform(*Deltas[g], g, 0)

			if (VERBOSE==5) {
				g
				(*Deltas[g]) // transform(transform(ZZZ,0,1), 1, g) - (*Zs[g])
				transform(transform(*Deltas[g],g,0),0,1)
				// 	(*Zs[g])
			}
		}
		
		if (VERBOSE>1) timer_off(40)
		// Optional: Acceleration
		// This is method 3 of Macleod (1986), a vector generalization of the Aitken-Steffensen method
		// Also: "when numerically computing the sequence.. stop..  when rounding errors become too 
		// important in the denominator, where the ^2 operation may cancel too many significant digits"

 		// Sometimes the iteration gets "stuck"; can we unstuck it with adding randomness in the accelerate decision?
		// There should be better ways too..
		
		if (accelerate_candidate) {
			mu = accelerate_norm = 0
			for (g=gstart;g<=num_fe;g++) {
				if (FEs[g].is_mock) continue
				mu = mu + quadcross( (*Deltas[g]) , (*Deltas[g]) - (*oldDeltas[g]) )
				accelerate_norm = accelerate_norm + norm((*Deltas[g]) - (*oldDeltas[g])) ^ 2
			}
			accelerate_norm = max((accelerate_norm, eps))
			//(iter, mu, accelerate_norm, mu/accelerate_norm, mean(*Zs[2]), mean(*Deltas[2]))
			mu = mu / accelerate_norm

			// Don't accelerate if mu is close to 0 (highly unlikely)
			if (abs(mu)>1e-6) {
				code = "a"
				accelerated = 1
				for (g=gstart;g<=num_fe;g++) {
					if (FEs[g].is_mock) continue
					(*Zs[g]) = (*Zs[g]) - mu :* (*Deltas[g])
					ZZZ = ZZZ - mu :* transform((*Deltas[g]), g, 0)
				}
			}
		} // accelerate_candidate
		
		// Reporting
		//update_error = (iter==1)? 1 : mean(reldif((*oldZs[g]), (*Zs[g])))
		update_error = 0
		for (g=gstart;g<=num_fe;g++) {
			if (FEs[g].is_mock) continue
			update_error = max(( update_error , mean(reldif( (*oldZs[g]) , (*Zs[g]) )) )) // max or mean?
			(*oldZs[g]) = (*Zs[g])
		}
		if (iter==1) update_error = 1
		//oldZZZ = ZZZ

		if ((VERBOSE>=2 & VERBOSE<=3 & mod(iter,1)==0) | (VERBOSE==1 & mod(iter,99)==0)) {
			printf(code)
			displayflush()
		}

		// Experimental: Pause acceleration when it seems stuck
		if (accelerated==1) {
			improvement = max(( (old_error-update_error)/update_error , (oldest_error-update_error)/update_error ))
			bad_loop = improvement < stuck_threshold ? bad_loop+1 : 0
			// bad_loop, improvement, update_error, old_error, oldest_error
			// Tolerate two problems (i.e. 6+2=8 iters) and then try to unstuck
			if (bad_loop>bad_loop_threshold) {
				bad_loop = 0
				if (VERBOSE==3) printf(" Fixed point iteration seems stuck, acceleration paused\n")
				acceleration_countdown = pause_length
			}
			assert(bad_loop<=3)	
			oldest_error = old_error
			old_error = update_error
		}

		if (VERBOSE>=2 & VERBOSE<=3 & mod(iter,99)==0) printf("%9.1f\n", update_error/tolerance)
		if (VERBOSE>=4) printf("%12.7e %1.0f \n", update_error, accelerate_candidate + accelerate_candidate*(acceleration_countdown==pause_length) ) // 0=Normal 1=Accel 2=BadAccel

		
		if ( (accelerated==0) & (update_error<tolerance) ) {
			converged = 1
			break
		}
	} // for

	if (VERBOSE>=2 & VERBOSE<=3 & mod(iter,99)!=0) printf("\n")
	if (!converged) {
		stata(sprintf(`"di as error "could not obtain resid of %s in %g iterations (last error=%e); to solve increase maxiter() or decrease tol().""', varname, max_iterations, update_error))
		exit(error(430))
	}
	if (VERBOSE>1) printf("{txt} - Converged in %g iterations (last error =%3.1e)\n", iter, update_error)
	if (VERBOSE==1) printf("{txt} converged in %g iterations, last error =%3.1e)\n", iter, update_error)
	if (VERBOSE>1) printf("{txt} - Saving output\n")

	// Recover Z1 = P1(y-ZZZ) where ZZZ=Z2+..+ZG
	Zs[1] = &transform(transform(y-stdev:*ZZZ, 0, 1), 1, 0)
	// Recover resid of y = y - ZZZ - Z1
	st_store(., resid_varname, y-stdev:*ZZZ-*Zs[1]) // BUGBUG if resid is just a vew, just do resid[.,.] = y-...
	
	// Save FEs
	if (save_fe!=0) {
		if (VERBOSE>1) printf("{txt} - Saving FEs\n")
		
		if (gstart==2) {
			st_store(., st_addvar("double", FEs[1].Z), *Zs[1])
		}
		else {
			st_store(., st_addvar("double", FEs[1].Z), betas[FEs[1].group,FEs[1].K] )
			for (k=2; k<=FEs[1].K; k++) {
				st_store(., st_addvar("double", FEs[k].Z), FEs[1].v[.,k-1] :* betas[FEs[1].group,k-1] )
			}
		}

		for (g=gstart;g<=num_fe;g++) {
			if (FEs[g].is_mock) continue

			if (!FEs[g].is_bivariate) {
				st_store(., st_addvar("double", FEs[g].Z), transform(stdev :* (*Zs[g]), g, 0))
			}
			else {
				(*Zs[g]) = transform((*Zs[g]), 0, g) // this saves -betas-
				st_store(., st_addvar("double", FEs[g].Z), stdev :* betas[FEs[g].group,FEs[g].K] )
				for (k=2; k<=FEs[g].K; k++) {
					st_store(., st_addvar("double", FEs[g+k-1].Z), stdev :* FEs[g].v[.,k-1] :* betas[FEs[g].group,k-1] )
				}
			}
		}

	}

	if (VERBOSE>1) {
		printf("{txt} - make_residual: inner loop took {res}%-6.2g \n\n", timer_value(40)[1])
		timer_clear(40)
		printf("")
	}

}

end
