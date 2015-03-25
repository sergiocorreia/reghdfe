*! hdfe 2.0.274 24mar2015
*! Sergio Correia (sergio.correia@duke.edu)
* (built from multiple source files using build.py)
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
	gstart = 1 + FEs[1].K
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

program define hdfe, rclass
	local version `=clip(`c(version)', 11.2, 13.1)' // 11.2 minimum, 13+ preferred
	qui version `version'

* Intercept version calls
	cap syntax, version
	local rc = _rc
	 if (`rc'==0) {
		Version
		exit
	}

* Intercept multiprocessor/parallel calls
	cap syntax, instance [*]
	local rc = _rc
	 if (`rc'==0) {
		ParallelInstance, `options'
		exit
	}

* Parse
	syntax varlist [if] [in] [fweight aweight pweight/] , Absorb(string) ///
		[PARTIAL(varlist numeric)] ///
		[CORES(integer 1)] ///
		[DROPSIngletons] ///
		[SAMPLE(name)] ///
		[GENerate(name)] [CLEAR] ///
		[CLUSTERVARs(string) Verbose(integer 0) TOLerance(real 1e-7) MAXITerations(integer 10000)] ///
		[SAVEFE] ///
		[*]

	Assert ("`generate'"!="") + ("`clear'"!="") == 1 , msg("hdfe error: you need to specify one and only one of the following options: clear generate(...)")

	opts_exclusive "(`partial') `savefe'"

	* Check that intersection(partial,varlist) = Null
	local intersection : list varlist & partial
	Assert "`intersection'"=="", msg("variables in varlist cannot appear in partial()")

	if ("`savefe'"!="") {
		local numvars : word count `varlist'
		Assert `numvars'==1 , msg("hdfe error: option savefe only allows one variable")
		local opt_savefe "save_fe(1)"
	}

	if ("`sample'"!="") conf new var `sample'

	if ("`weight'"!="") {
		local weightvar `exp'
		conf var `weightvar' // just allow simple weights
		local weighttype `weight'
		local weightequal =
	}

* Preserve if asked to
	if ("`generate'"!="") {

		* The stub must not exist!
		 cap ds `generate'*
		 Assert "`r(varlist)'"=="", msg("hdfe error: there are already variables that start with the stub `generate'")

		tempvar uid
		gen double `uid' = _n
		preserve
	}

* Clear previous errors
	Stop

* Time/panel variables
	cap conf var `_dta[_TStvar]'
	if (!_rc) local timevar `_dta[_TStvar]'
	cap conf var `_dta[_TSpanel]'
	if (!_rc) local panelvar `_dta[_TSpanel]'

* Set Verbosity
	mata: VERBOSE = `verbose' // Pick a number between 0 (quiet) and 4 (lots of debugging info)

* Parse: absorb, clusters, and weights
	Start, absorb(`absorb') clustervars(`clustervars') weight(`weighttype') weightvar(`weightvar')
	local absorb_keepvars = r(keepvars)
	local N_hdfe = r(N_hdfe)
	
* Keep relevant observations
	marksample touse, novar
	markout `touse' `varlist' `partial' `absorb_keepvars'
	qui keep if `touse'
	
* Keep relevant variables
	keep `varlist' `partial' `clustervars' `weightvar' `panelvar' `timevar' `uid' `absorb_keepvars'

* Drop singletons
	if ("`dropsingletons'"!="") DropSingletons, num_absvars(`N_hdfe')
	
* Construct Mata objects and auxiliary variables
	Precompute, ///
		keep(`varlist' `partial' `clustervars' `weightvar' `panelvar' `timevar' `uid') ///
		tsvars(`panelvar' `timevar')
	
* Compute e(df_a)
	EstimateDoF, dofadjustments(pairwise clusters continuous)
	* return list // what matters is r(kk) which will be e(df_a)
	local kk = r(kk)
	forval g = 1/`N_hdfe' {
		local df_a`g' = r(K`g') - r(M`g')
	}
	
* We don't need the FE variables (they are in mata objects now)
	*drop __FE*__

* Demean variables wrt to the fixed effects
	local opt varlist(`varlist' `partial') tol(`tolerance') maxiterations(`maxiterations') `options' `opt_savefe'
	if (`cores'>1) {
		DemeanParallel, `opt' self(hdfe) cores(`cores')
	}
	else {
		Demean, `opt'
	}

	if ("`savefe'"!="") {
		Save, original_depvar(`varlist')
		local saved_fe = r(keepvars)
	}

	return scalar df_a = `kk'
	return scalar N_hdfe = `N_hdfe'
	forv g=1/`N_hdfe' {
		mata: fe2local(`g') // copies Mata structure into locals
		* Will inject the following with c_local:
		* ivars cvars target varname varlabel is_interaction is_cont_interaction is_bivariate is_mock levels
		return local hdfe`g' = "`varlabel'"
		return scalar df_a`g' = `df_a`g'' // `levels'
	}
	
* Clean up Mata objects
	Stop

	if ("`partial'"!="") {
		tempvar resid
		_rmcoll `partial', forcedrop
		local partial = r(varlist)
		foreach var of local varlist {
			_regress `var' `partial' `weightexp' [`weighttype'`weightequal'`weightvar'], nohead notable
			_predict double `resid', resid
			qui replace `var' = `resid' // preserve labels
			drop `resid'
		}
		local numpartial : word count `partial'
		return scalar df_partial = `numpartial'
	}

	if ("`generate'"!="") {
		keep `varlist' `uid' `saved_fe'
		foreach var of local varlist {
			rename `var' `generate'`var'
		}

		tempfile output
		sort `uid'
		qui save "`output'"
		restore
		SafeMerge, uid(`uid') file("`output'") sample(`sample')
	}
end

* [SafeMerge: ADAPTED FROM THE ONE IN ESTIMATE.ADO]
* The idea of this program is to keep the sort order when doing the merges

program define SafeMerge, eclass sortpreserve
syntax, uid(varname numeric) file(string) [sample(string)]
	* Merging gives us e(sample) and the FEs / AvgEs
	if ("`sample'"!="") {
		tempvar smpl
		merge 1:1 `uid' using "`file'", assert(master match) nolabel nonotes noreport gen(`smpl')
		gen byte `sample' = (`smpl'==3)
		drop `smpl' // redundant
	}
	else {
		merge 1:1 `uid' using "`file'", assert(master match) nolabel nonotes noreport nogen
	}
end


// -------------------------------------------------------------
// Simple assertions
// -------------------------------------------------------------

program define Assert
    syntax anything(everything equalok) [, MSG(string asis) RC(integer 198)]
    if !(`anything') {
        di as error `msg'
        exit `rc'
    }
end


// -------------------------------------------------------------
// Simple debugging
// -------------------------------------------------------------

program define Debug

	syntax, [MSG(string asis) Level(integer 1) NEWline COLOR(string)] [tic(integer 0) toc(integer 0)]
	
	cap mata: st_local("VERBOSE",strofreal(VERBOSE)) // Ugly hack to avoid using a global
	if ("`VERBOSE'"=="") {
		di as result "Mata scalar -VERBOSE- not found, setting VERBOSE=3"
		local VERBOSE 3
		mata: VERBOSE = `VERBOSE'
	}


	assert "`VERBOSE'"!=""
	assert inrange(`level',0, 4)
	assert (`tic'>0) + (`toc'>0)<=1

	if ("`color'"=="") local color text
	assert inlist("`color'", "text", "res", "result", "error", "input")

	if (`VERBOSE'>=`level') {

		if (`tic'>0) {
			timer clear `tic'
			timer on `tic'
		}
		if (`toc'>0) {
			timer off `toc'
			qui timer list `toc'
			local time = r(t`toc')
			if (`time'<10) local time = string(`time'*1000, "%tcss.ss!s")
			else if (`time'<60) local time = string(`time'*1000, "%tcss!s")
			else if (`time'<3600) local time = string(`time'*1000, "%tc+mm!m! SS!s")
			else if (`time'<24*3600) local time = string(`time'*1000, "%tc+hH!h! mm!m! SS!s")
			timer clear `toc'
			local time `" as result " `time'""'
		}

		if (`"`msg'"'!="") di as `color' `msg'`time'
		if ("`newline'"!="") di
	}
end


// -------------------------------------------------------------
// Report HDFE/REGHDFE version
// -------------------------------------------------------------

program define Version, eclass
    local version "2.0.274 24mar2015"
    ereturn clear
    di as text "`version'"
    ereturn local version "`version'"

end



// -------------------------------------------------------------
// Faster alternative to -makegps-, but with some limitations
// -------------------------------------------------------------
* To avoid backuping the data, use option -clear-
* For simplicity, disallow -if- and -in- options

program ConnectedGroups, rclass
syntax varlist(min=2 max=2) [, GENerate(name) CLEAR]

* To avoid backuping the data, use option -clear-
* For simplicity, disallow -if- and -in- options

    if ("`generate'"!="") conf new var `generate'
    gettoken id1 id2 : varlist
    Debug, level(2) msg("    - computing connected groups between `id1' and`id2'")
    tempvar group copy

    tempfile backup
    if ("`clear'"=="") qui save "`backup'"
    keep `varlist'
    qui bys `varlist': keep if _n==1

    clonevar `group' = `id1'
    clonevar `copy' = `group'
    capture error 100 // We want an error
    while _rc {
        qui bys `id2' (`group'): replace `group' = `group'[1]
        qui bys `id1' (`group'): replace `group' = `group'[1]
        capture assert `copy'==`group'
        qui replace `copy' = `group'
    }

    assert !missing(`group')
    qui bys `group': replace `group' = (_n==1)
    qui replace `group' = sum(`group')
    
    su `group', mean
    local num_groups = r(max)
    
    if ("`generate'"!="") rename `group' `generate'
    
    if ("`clear'"=="") {
        if ("`generate'"!="") {
            tempfile groups
            qui compress
            la var `generate' "Mobility group for (`varlist')"
            qui save "`groups'"
            qui use "`backup'", clear
            qui merge m:1 `id1' `id2' using "`groups'" , assert(match) nogen
        }
        else {
            qui use "`backup'", clear
        }
    }
    
    return scalar groups=`num_groups'
end


// -------------------------------------------------------------
// Faster alternative to -egen group-. MVs, IF, etc not allowed!
// -------------------------------------------------------------

program define GenerateID, sortpreserve
syntax varlist(numeric) , [REPLACE Generate(name)]

	assert ("`replace'"!="") + ("`generate'"!="") == 1
	// replace XOR generate, could also use -opts_exclusive -
	foreach var of varlist `varlist' {
		assert !missing(`var')
	}

	local numvars : word count `varlist'
	if ("`replace'"!="") assert `numvars'==1 // Can't replace more than one var!
	
	// Create ID
	tempvar new_id
	sort `varlist'
	by `varlist': gen long `new_id' = (_n==1)
	qui replace `new_id' = sum(`new_id')
	qui compress `new_id'
	assert !missing(`new_id')
	
	local name = "i." + subinstr("`varlist'", " ", "#i.", .)
	char `new_id'[name] `name'
	la var `new_id' "[ID] `name'"

	// Either replace or generate
	if ("`replace'"!="") {
		drop `varlist'
		rename `new_id' `varlist'
	}
	else {
		rename `new_id' `generate'
	}

end


// -------------------------------------------------------------
// AvgE: Average of all the other obs in a group, except each obs itself
// -------------------------------------------------------------

program define AverageOthers , sortpreserve
syntax varname , BY(varlist) Generate(name) [EXCLUDESELF]

* EXCLUDESELF: Excludes obs at hand when computing avg

***[EXCLUDE(varname)]
*** Do not use obs where `exclude'!=0 to compute the means, but do fill out these values

* Alternative:
* MeanOthers = MeanAll * N/(N-1) - X / (N-1) = (SumAll-X)/(N-1)
* Also, using mean() instead of total() would give less rounding errors

	sort `by'

	conf new var `generate'
	***if ("`exclude'"!="") local cond " if !`exclude'"
	
	* Calculate avg by group
	tempvar total count
	qui gen double `generate' = `varlist' `cond'
	
	* Sum
	*qui by `by' : egen double `generate' = mean(`var')
	qui by `by' : gen double `total' = sum(`generate')
	qui by `by' : replace `total' = `total'[_N]
	
	* Count
	qui by `by' : gen double `count' = sum(`generate'<.)
	qui by `by' : replace `count' = `count'[_N]
	
	* Substract itself
	if ("`excludeself'"!="") qui by `by' : replace `total' = `total' - `generate' if (`generate'<.)
	if ("`excludeself'"!="") qui by `by' : replace `count' = `count' - 1 if (`generate'<.)
	
	* Divide
	qui replace `generate' = `total' / `count'
	
	**qui by `by' : replace `generate' = `generate'[_N]
	
	* Adjust negative values b/c of rounding errors introduced by -excludeself- (risky)
	if ("`excludeself'"!="") {
		replace `generate' = 0 if inrange(`generate', -1e-8, 0)
		local note X
	}

	* Add label and chars
	local name = subinstr("`by'", " ", "_", .)
	char `generate'[avge_equation]  AvgE`note'
	char `generate'[name] `name'
	char `generate'[depvar] `varlist'
	la var `generate' "Avg`note'. of `varlist' by `by'"
end


// -------------------------------------------------------------------------------------------------
// Calculate the degrees of freedom lost due to the absorbed fixed effects
// -------------------------------------------------------------------------------------------------
/*
	In general, we can't know the exact number of DoF lost because we don't know when multiple FEs are collinear
	When we have two pure FEs, we can use an existing algorithm, but besides that we'll just use an upper (conservative) bound

	Features:
	 - Save the first mobility group if asked
	 - Within the pure FEs, we can use the existing algorithm pairwise (FE1 vs FE2, FE3, .., FE2 vs FE3, ..)
	 - If there are n pure FEs, that means the algo gets called n! times, which may be kinda slow
	 - With FEs interacted with continuous variables, we can't do this, but can do two things:
		a) With i.a#c.b , whenever b==0 for all values of a group (of -a-), add one redundant
		b) With i.a##c.b, do the same but whenever b==CONSTANT (so not just zero)
     - With clusters, it gets trickier but in summary you don't need to penalize DoF for params that only exist within a cluster. This happens:
		a) if absvar==clustervar
		b) if absvar is nested within a clustervar. EG: if we do vce(cluster state), and -absorb(district)- or -absorb(state#year)
		c) With cont. interactions, e.g. absorb(i.state##c.year) vce(cluster year), then i) state FE is redundant, but ii) also state#c.year
		   The reason is that at the param for each "fixed slope" is shared only within a state

	Procedure:
	 - Go through all FEs and see if i) they share the same ivars as any clusters, and if not, ii) if they are nested within clusters
	 - For each pure FE in the list, run the algorithm pairwise, BUT DO NOT RUN IT BEETWEEN TWO PAIRS OF redundant
	   (since the redundants are on the left, we just need to check the rightmost FE for whether it was tagged)
	 - For the ones with cont interactions, do either of the two tests depending on the case

	Misc:
	 - There are two places where DoFs enter in the results:
		a) When computing e(V), we do a small sample adjustment (seen in Stata documentation as the -q-)
		   Instead of doing V*q with q = N/(N-k), we use q = N / (N-k-kk), so THE PURPOSE OF THIS PROGRAM IS TO COMPUTE "kk"
		   This kk will be used to adjust V and also stored in e(df_a)
		   With clusters, q = (N-1) / (N-k-kk) * M / (M-1)
		   With multiway clustering, we use the smallest N_clust as our M
	    b) In the DoF of the F and t tests (not when doing chi/normal)
	       When there are clusters, note that e(df_r) is M-1 instead of N-1-k
	       Again, here we want to use the smallest M with multiway clustering

	Inputs: +-+- if we just use -fe2local- we can avoid passing stuff around when building subroutines
	 - We need the current name of the absvars and clustervars (remember a#b is replaced by something different)
	 - Do a conf var at this point to be SURE that we didn't mess up before
	 - We need the ivars and cvars in a list
	 - For the c. interactions, we need to know if they are bivariate or univariate
	 - SOLN -> mata: fe2local(`g')  ; from mata: ivars_clustervar`i' (needed???) , and G
	 - Thus, do we really needed the syntax part??
	 - fe2local saves: ivars cvars target varname varlabel is_interaction is_cont_interaction is_bivariate is_mock levels // Z group_k weightvar

	DOF Syntax:
	 DOFadjustments(none | all | CLUSTERs | PAIRwise | FIRSTpair | CONTinuous)
	 dof() = dof(all) = dof(cluster pairwise continuous)
	 dof(none) -> do nothing; all Ms = 0 
	 dof(first) dof(first cluster) dof(cluster) dof(continuous)

	For this to work, the program MUST be modular
*/

program define EstimateDoF, rclass
syntax, [DOFadjustments(string) group(name) uid(varname) groupdta(string)]
	
	* Parse list of adjustments/tricks to do
	Debug, level(1) msg("(calculating degrees of freedom lost due to the FEs)")
	local adjustement_list firstpairs pairwise clusters continuous
	* This allows doing things like <if (`adj_clusters') ..>
	Debug, level(2) msg(`" - Adjustments:"')
	foreach adj of local adjustement_list {
		local adj_`adj' : list posof "`adj'" in dofadjustments
		Debug, level(2) msg(`"    - `adj' {col 18}{res} `=cond(`adj_`adj'',"yes","no")'"')
	}

	* Assert that the clustervars exist
	mata: st_local("clustervars", invtokens(clustervars))
	conf variable `clustervars', exact

	mata: st_local("G", strofreal(G))
	mata: st_local("N_clustervars", strofreal(length(clustervars)))

	if ("`group'"!="") {
		Assert (`adj_firstpairs' | `adj_pairwise'), msg("Cannot save connected groups without options pairwise or firstpair")
	}

	* Remember: fe2local stores the following:
	* ivars cvars target varname varlabel is_interaction is_cont_interaction is_bivariate is_mock levels

* Starting point assumes no redundant parameters
	forv g=1/`G' {
		mata: fe2local(`g')
		local redundant`g' = 0 // will be 1 if we don't penalize at all for this absvar (i.e. if it's nested with cluster or collinear with another absvar)
		local is_slope`g' = ("`cvars'"!="") & (!`is_bivariate' | `is_mock') // two cases: i.a#c.b , i.a##c.b (which expands to <i.a i.a#c.b> and we want the second part)
		local M`g' = !`is_slope`g'' // Start with 0 with cont. interaction, 1 w/out cont interaction

		*For each FE, only know exactly parameters are redundant in a few cases:
		*i) nested in cluster, ii) first pure FE, iii) second pure FE if checked with connected groups
		local exact`g' 0
		local drop`g' = !(`is_bivariate' & `is_mock')
	}

* Check if an absvar is a clustervar or is nested in a clustervar
* We *always* check if absvar is a clustervar, to prevent deleting its __FE__ variable by mistake
* But we only update the DoF if `adj_clusters' is true.

	local M_due_to_nested 0 // Redundant DoFs due to nesting within clusters
	if (`N_clustervars'>0) {
		mata: st_local("clustervars", invtokens(clustervars))
		forv g=1/`G' {
			mata: fe2local(`g')
			local gg = `g' - `is_mock'
			local absvar_in_clustervar 0 // 1 if absvar is nested in a clustervar
			
			* Trick: if the absvar is also a clustervar, then its name will be __FE*__
			local absvar_is_clustervar : list varname in clustervars
			if (`adj_clusters' & `absvar_is_clustervar') {
				Debug, level(1) msg("(categorical variable " as result "`varlabel'"as text " is also a cluster variable, so it doesn't count towards DoF)")
			}
			else if (`adj_clusters') {
				forval i = 1/`N_clustervars' {
					mata: st_local("clustervar", clustervars[`i'])
					mata: st_local("clustervar_original", clustervars_original[`i'])
					cap _xtreg_chk_cl2 `clustervar' __FE`gg'__
					assert inlist(_rc, 0, 498)
					if (!_rc) {
						Debug, level(1) msg("(categorical variable " as result "`varlabel'" as text " is nested within cluster variable " as result "`clustervar_original'" as text ", so it doesn't count towards DoF)")
						continue, break
					}
				}
			}

			if (`absvar_is_clustervar') local drop`g' 0

			if ( `adj_clusters' & (`absvar_is_clustervar' | `absvar_in_clustervar') ) {
				local M`g' = `levels'
				local redundant`g' 1
				local exact`g' 1
				local M_due_to_nested = `M_due_to_nested' + `levels' - 1
			}
		} // end for over absvars
	} // end cluster adjustment

* Just indicate the first pure FE that is not nested in a cluster
	forv g=1/`G' {
		if (!`is_slope`g'' & !`redundant`g'') {
			local exact`g' 1
			continue, break
		}
	}

* Compute connected groups for the remaining FEs (except those with cont interactions)

	local dof_exact 0 // if this code never runs, it's not exact
	if (`adj_firstpairs' | `adj_pairwise') {
		Debug, level(3) msg(" - Calculating connected groups for DoF estimation")
		local dof_exact 1
		local i_comparison 0
		forv g=1/`G' {
			if (`is_slope`g'') local dof_exact 0 // We may not get all redundant vars with cont. interactions
			if (`is_slope`g'') continue
			local start_h = `g' + 1
			forv h=`start_h'/`G' {

				if (`is_slope`h'' | `redundant`h'') continue
				local ++i_comparison
				if (`i_comparison'>1) local dof_exact 0 // Only exact with one comparison
				if (`i_comparison'>1 & `adj_firstpairs') continue // -firstpairs- will only run the first comparison
				if (`i_comparison'==1) local exact`h' 1

				* ConnectedGroups does destructive operations and thus backups the dta by default
				* This is very slow with huge datasets and e.g. 4 FEs (up to 3*2*1=6 saves).
				* As a soln, use the -clear- opt and save before. Rule:
				* - Save the cache on the first comparison, OR if we are saving the connected group, on the second

				if (`i_comparison'==1 & "`group'"!="") {
					ConnectedGroups __FE`g'__ __FE`h'__ , gen(`group')
				}
				else if (`i_comparison'==1 & "`group'"=="") | (`i_comparison'==2 & "`group'"!="") {
					tempfile backup
					qui save "`backup'"
					ConnectedGroups __FE`g'__ __FE`h'__ , clear
					qui use "`backup'", clear
				}
				else {
					ConnectedGroups __FE`g'__ __FE`h'__ , clear
					qui use "`backup'", clear
				}

				local candidate = r(groups)
				local M`h' = max(`M`h'', `candidate')
			}
		}
	} // end connected group comparisons

* Adjustment with cont. interactions
	if (`adj_continuous') {
		forv g=1/`G' {
			mata: fe2local(`g')
			if (!`is_slope`g'') continue
			CheckZerosByGroup, fe(`varname') cvars(`cvars') anyconstant(`is_mock')
			local M`g' = r(redundant)
		}
	}

	if (`dof_exact') {
		Debug, level(1) msg(" - DoF computation is exact")
	}
	else {
		Debug, level(1) msg(" - DoF computation not exact; DoF may be higher than reported")	
	}

	local SumM 0
	local SumK 0
	Debug, level(2) msg(" - Results of DoF adjustments:")
	forv g=1/`G' {
		mata: fe2local(`g')
		assert !missing(`M`g'') & !missing(`levels')
		local SumM = `SumM' + `M`g''
		local SumK = `SumK' + `levels'

		return scalar M`g' = `M`g''
		return scalar K`g' = `levels'
		return scalar M`g'_exact = `exact`g''
		return scalar drop`g' = `drop`g''
		Debug, level(2) msg("   - FE`g' ({res}`varlabel'{txt}): {col 40}K=`levels' {col 50}M=`M`g'' {col 60}is_exact=`exact`g''")
	}
	return scalar M = `SumM'
	local NetSumK = `SumK' - `SumM'
	Debug, level(2) msg(" - DoF loss due to FEs: Sum(Kg)=`SumK', M:Sum(Mg)=`SumM' --> KK:=SumK-SumM=`NetSumK'")
	return scalar kk = `NetSumK'

* Save mobility group if needed
	local saved_group = 0
	if ("`group'"!="") {
		conf var `group'
		tempfile backup
		qui save "`backup'"
		
		keep `uid' `group'
		sort `uid'
		la var `group' "Mobility group between `label'"
		qui save "`groupdta'" // A tempfile from the caller program
		Debug, level(2) msg(" - mobility group saved")
		qui use "`backup'", clear
		cap erase "`backup'"
		local saved_group = 1
	}
	return scalar saved_group = `saved_group'
	return scalar M_due_to_nested = `M_due_to_nested'
end

program define CheckZerosByGroup, rclass sortpreserve
syntax, fe(varname numeric) cvars(varname numeric) anyconstant(integer)
	tempvar redundant
	assert inlist(`anyconstant', 0, 1)
	if (`anyconstant') {
		qui bys `fe' (`cvars'): gen byte `redundant' = (`cvars'[1]==`cvars'[_N]) if (_n==1)
	}
	else {
		qui bys `fe' (`cvars'): gen byte `redundant' = (`cvars'[1]==0 & `cvars'[_N]==0) if (_n==1)
	}
	qui cou if `redundant'==1
	return scalar redundant = r(N)
end

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

program define ParseOneAbsvar, rclass
	syntax, ABSVAR(string)

	Assert !strpos("`absvar'","###"), msg("error parsing <`absvar'> : ### is invalid")
	Assert regexm("`absvar'", "^[a-zA-Z0-9_=.#]+$"), msg("error parsing <`absvar'> : illegal characters ")
	Assert !regexm("`absvar'", "##([^c]|(c[^.]))"), msg("error parsing <`absvar'> : expected c. after ##")
	local original_absvar `absvar'

* Split at equal sign
	local equalsign = strpos("`absvar'","=")
	local target = substr("`absvar'",1,`equalsign'-1)
	local absvar = substr("`absvar'",`equalsign'+1, .)
	if ("`target'"!="") conf new var `target'

	local is_interaction = strpos("`absvar'", "#")>0
	local is_bivariate = strpos("`absvar'", "##")>0

* Split interactions
	mata: st_local("vars", subinstr("`absvar'", "#", " ") )
	foreach var of local vars {

		local dot = strpos("`var'", ".")
		local root = substr("`var'", `dot'+1, .)
		unab root : `root' , max(1)
		conf numeric var `root'
		
		local prefix = substr("`var'", 1, `dot'-1)
		local prefix = lower( cond("`prefix'"=="", "i", "`prefix'") ) // -i.- is default prefix

		Assert inlist("`prefix'", "i", "c") , msg("error parsing <`absvar'><`var'> : only i. and c. allowed, not `prefix'.")
		Assert !strpos("`root'", ".") , msg("error parsing <`absvar'><`var'> : no time series operators allowed")
		
		if ("`prefix'"=="i") {
			local ivars `ivars' `root'
		}
		else {
			Assert "`cvars'"=="", msg("error: can't have more than one continuous variable in the interaction")
			local cvars `cvars' `root'
		}
	}
	local tab  "        "
	Debug, level(3) msg(as text "    Parsing " as result "`original_absvar'")
	Debug, level(3) msg(as text "`tab'ivars = " as result "`ivars'")
	if ("`cvars'"!="") Debug, level(3) msg(as text "`tab'cvars = " as result "`cvars'")
	if ("`target'"!="") Debug, level(3) msg(as text "`tab'target = " as result "`target'")
	Debug, level(3) msg(as text "`tab'is_interaction = " as result "`is_interaction'")
	Debug, level(3) msg(as text "`tab'is_bivariate = " as result "`is_bivariate'")
	// Debug, level(3) newline

	return scalar is_interaction = `is_interaction'
	return scalar is_cont_interaction = `is_interaction' & ("`cvars'"!="")
	return scalar is_bivariate = `is_bivariate'
	if ("`target'"!="") return local target "`target'"
	if ("`cvars'"!="") return local cvars "`cvars'"
	return local ivars "`ivars'"
end

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

program define Demean

	CheckCorrectOrder demean
	syntax , VARlist(varlist numeric) ///
		[TOLerance(real 1e-7) MAXITerations(integer 10000) ACCELerate(integer 1) /// See reghdfe.Parse
		CHECK(integer 0) SAVE_fe(integer 0) /// Runs regr of FEs
		NUM_fe(integer -1)] /// Regress only against the first Nth FEs (used in nested Fstats)
		[bad_loop_threshold(integer 1) stuck_threshold(real 5e-3) pause_length(integer 20) ///
		accel_freq(integer 3) accel_start(integer 6)] /// Advanced options

	assert inrange(`tolerance', 1e-20, 1) // However beyond 1e-16 we reach the limits of -double-
	assert inrange(`maxiterations',1,.)
	assert inlist(`accelerate',0,1)
	assert inlist(`check',0,1)
	assert inlist(`save_fe',0,1)
	assert inrange(`num_fe',1,100) | `num_fe'==-1 // -1 ==> Use all FEs

	assert `bad_loop_threshold'>0
	assert `stuck_threshold'>0 & `stuck_threshold'<=1
	assert `pause_length'>=0
	assert `accel_freq'>=0
	assert `accel_start'>0

	* We need to recast everything to -double- (-float- is not good enough)
	Debug, level(2) msg("(recasting variables as -double-)")
	recast double `varlist'

	* We can't save the FEs if there is more than one variable
	cap unab _ : `varlist', max(1)
	Assert (_rc==0 | `save_fe'==0) , rc(`=_rc') ///
		msg("hdfe.Demean: cannot save FEs of more than one variable at a time")

	tempvar resid
	local save = `save_fe' | `check' // check=1 implies save_fe=1
	local base_args `""`resid'", `tolerance', `maxiterations', `save', `accelerate', `num_fe'"'
	local adv_args `"`bad_loop_threshold', `stuck_threshold', `pause_length', `accel_freq', `accel_start'"'
	local args `"`base_args', `adv_args'"'
	Debug, level(3) msg(" - Structure of Mata calls: make_residual(" as result "{variable}" as text `", `args')"')

	Debug, level(2) tic(30)
	mata: st_local("weightexp", weightexp)
	
	foreach var of varlist `varlist' {
		cap drop __Z*__
		Assert !missing(`var'), msg("hdfe.Demean error: `var' has missing values and cannot be transformed")
		
		* Syntax: MAKE_RESIDUAL(var, newvar, tol, maxiter | , save=0 , accel=1, first_n=`num_fe')
		* Note: summarize doesn't allow pweight ( see http://www.stata.com/support/faqs/statistics/weights-and-summary-statistics/ )
		* Since we only want to compute means, replace with [aw]
		local tmpweightexp = subinstr("`weightexp'", "[pweight=", "[aweight=", 1)
		qui su `var' `tmpweightexp', mean
		char define `var'[mean] `r(mean)'
		mata: make_residual("`var'", `args')
		assert !missing(`resid')

		* Check that coefs are approximately 1
		if (`check') {
			unab _ : __Z*__, min(1)
			local backup = ("`e(cmd)'"!="")
			if (`backup') {
				tempname backup_results
				est store `backup_results', nocopy // nocopy needed to avoid having e(_estimates_name)
			}
			qui _regress `var' __Z*__
			local label : var label `var'
			if ("`label'"=="") local label `var'
			di as text "FE coefficients for `label':{col 36}" _continue
			foreach z of varlist __Z*__ {
				assert !missing(`z')
				di as text " `=string(_b[`z'], "%9.7f")'"  _continue
			}
			di
			
			if (`backup') qui est restore `backup_results'
			if (!`save_fe') cap drop __Z*__
		}

		* If the tol() is not high enough (e.g. 1e-14), we may fail to detect variables collinear with the absorbed categories
		* Again, we can't use pweight with summarize, but in this case it's just for debugging purposes so use [aw]
		qui su `resid' `tmpweightexp'
		local prettyvar `var'
		if (substr("`var'", 1, 2)=="__") local prettyvar : var label `var'
		if inrange(r(sd), 1e-20 , epsfloat()) di in ye "(warning: variable `prettyvar' is probably collinear, maybe try a tighter tolerance)"

		qui replace `var' = `resid' // This way I keep labels and so on
		drop `resid'
		Assert !missing(`var'), msg("REGHDFE.Demean: `var' has missing values after transformation")
	}
	Debug, level(2) toc(30) msg("(timer for calls to mata:make_residual)")
end

program define DemeanParallel
	* Notes:
	* First cluster is taking by this stata instance, to save HDD/memory/merge time
	* Also, this cluster should have more obs than the other ones so we let it have
	* the default number of processes
	* (the other start with 1 proc allowed, which should be fine)
	* Thus it will usually finish faster, to start waiting for the 2nd fastest  to merge

	CheckCorrectOrder demean
	syntax, VARlist(varlist numeric) FILEname(string) UID(varname numeric) CORES(integer) SELF(string) [*]

	local varlist : list uniq varlist
	local K : list sizeof varlist
	local cores = min(`cores',`K')
	local size = c(N) * c(width) / 2^30
	local wait = int(100 + 1000 * `size') // each gb wait 1 sec

	assert inlist("`self'", "reghdfe", "hdfe") // We will call `self', instance ...

	* Deal each variable like cards in Poker
	local core 1
	foreach var of local varlist {
		local varlist`core' `varlist`core'' `var'
		local ++core
		if (`core'>`cores') local core 1
	}

	* Folder name.. need some entropy.. use varlist + time
	mata: st_local("hash", strofreal(hash1("`varlist'"), "%20.0f"))
	local seed = real(subinstr(c(current_time),":","",.)) + `hash'
	local seed = mod(`seed',2^30) // Needs to be < 2^31-1
	set seed `seed'
	local code = string(int( uniform() * 1e6 ), "%08.0f")

	* Prepare
	* Note: On windows, tmpdir has a trailing / but on Linux it doesn't!
	local path "`c(tmpdir)'`c(dirsep)'hdfe_`code'`c(dirsep)'"
	* NOTE: Copy any changes to this line into ParallelInstance.ado

	Debug, level(1) msg(" - tempdir will be " as input "`path'")
	mata: parallel_cores = `cores'
	mata: parallel_dta = `"`filename'"'
	mata: parallel_vars = J(`cores',1,"")
	mata: parallel_opt = `"`options'"'
	mata: parallel_path = `"`path'"'
	forv i=1/`cores' {
		mata: parallel_vars[`i'] = "`varlist`i''"
	}

	local dropvarlist : list varlist - varlist1
	drop `dropvarlist' // basically, keeps UID and clustervar
	mata: st_global("hdfe_pwd",pwd())
	mkdir "`path'"
	qui cd "`path'"

	local objects VERBOSE G FEs betas prev_numstep parallel_* weightexp weightvar
	qui mata: mata matsave "`path'hdfe_mata.mo" `objects' , replace

	* Call -parallel-
	Debug, level(1) msg(" - running parallel instances")
	*qui mata: parallel_setstatadir("")
	qui parallel setclusters 1 // I just want the global PLL_DIR
	local binary `"$PLL_DIR"'
	Assert `"`binary'"'!="", msg("`self' error: after parallel, global PLL_DIR was empty")
	global PLL_DIR
	global PLL_CLUSTERS
	global PLL_STATA_PATH

	cap mata: st_local("VERBOSE",strofreal(VERBOSE))
	if (`VERBOSE'==0) local qui qui
	`qui' di as text _n "{dup 44:_}/ PARALLEL \{dup 44:_}"

	local flag = cond("`c(os)'"=="Windows", "/q", "-q")

	* Create instances
	forv i=2/`cores' {
		local cmd `"winexec `binary' `flag'  `self', instance core(`i') code(`code') "'
		Debug, level(1) msg(" - Executing " in ye `"`cmd' "')
		`cmd'
		Debug, level(1) msg(" - Sleeping `wait'ms")
		if (`i'!=`cores') sleep `wait'
	}
	Demean, varlist(`varlist1') `options' // core=1

	* Wait until all instances have started
	local timeout 20
	local elapsed 0
	forv i=2/`cores' {
		local ok 0
		while !`ok' {
			sleep 100
			local fn "`path'`i'_started.txt"
			cap conf file "`fn'"
			local rc = _rc
			if (`rc'==0) {
				local ok 1
				Debug, level(1) msg(" - process `i' started")
				erase "`fn'"
			}
			else {
				local elapsed = `elapsed' + 0.1
				Assert `elapsed'<`timeout', msg("Failed to start subprocess `i'")
			}
		}
		local remaining_cores `remaining_cores' `i' // Will contain remaining cores
	}

	* Wait for termination and merge
	while ("`remaining_cores'"!="") {
		foreach core of local remaining_cores {
			local donefn "`path'`core'_done.txt"
			local okfn "`path'`core'_ok.txt"
			local errorfn "`path'`core'_error.txt"
			local dtafn "`path'`core'_output.dta"
			local logfile "`path'`core'_log.log"


			cap conf file "`donefn'"
			local rc = _rc

			if (`rc'==0) {
				Debug, level(1) msg(" - process `core' finished")
				erase "`donefn'"
				cap conf file "`okfn'"
				if (`=_rc'>0) {
					type "`logfile'"
					//di as error "<`dtafn'> not found"
					Assert 0, msg("Call to subprocess `core' failed, see logfile")
				}

				erase "`okfn'"
				Debug, level(1) msg(" - Subprocess `core' done")
				local remaining_cores : list remaining_cores - core
				mata: st_local("VERBOSE",strofreal(VERBOSE))
				
				if (`VERBOSE'>=3) {
					type "`logfile'"
				}
				erase "`logfile'"

				* Merge file
				Debug, level(1) msg(" - Merging dta #`core'")
				merge 1:1 _n using "`dtafn'", nogen nolabel nonotes noreport sorted assert(match)
				erase "`dtafn'"
			}
			else {
				sleep 500 // increase this
			}
		}
	}

	* Cleanup
	qui cd "${hdfe_pwd}"
	erase "`path'hdfe_mata.mo"
	cap rmdir `"`path'"'
	`qui' di as text _n "{dup 44:_}\ PARALLEL /{dup 44:_}"

end

program define ParallelInstance
	syntax, core(integer) code(string asis)
	set more off
	assert inrange(`core',1,32)
	local path "`c(tmpdir)'`c(dirsep)'hdfe_`code'`c(dirsep)'"
	cd "`path'"
	set processors 1

	file open fh using "`core'_started.txt" , write text all
	file close _all

	cap noi {
		set linesize 120
		log using `core'_log.log, text

		mata: mata matuse "hdfe_mata.mo"
		mata: st_local("cores",strofreal(parallel_cores))
		assert `core' <= `cores'
		mata: st_local("usedta",parallel_dta)
		mata: st_local("vars",parallel_vars[`core'])
		mata: st_local("weightvar",weightvar)
		mata: st_local("opt",parallel_opt)
		Debug, msg(" - This is core `core'/`cores'")
		sleep 100
	
		local outfn "`core'_output.dta"
		conf new file "`outfn'"

		use `vars' `weightvar' using "`usedta'"
		de, full
		Demean, varlist(`vars') `opt'
		keep `vars'
		save `"`outfn'"'
		log close _all
	}

	local rc = _rc
	sleep 100

	if `rc'>0 {
		di in red "ERROR: `rc'"
		file open fh using "`core'_error.txt" , write text all
		file close _all
	}
	else {
		file open fh using "`core'_ok.txt" , write text all
		file close _all
	}

	file open fh using "`core'_done.txt" , write text all
	file close _all
	exit, STATA
end

program define Save, rclass
	* Run this after -Demean .. , save_fe(1)-
	* For each FE, if it has a -target-, add label, chars, and demean or divide
	CheckCorrectOrder save
	syntax , original_depvar(string)

	mata: st_local("G", strofreal(G))
	mata: st_local("weightexp", weightexp)
	forv g=1/`G' {

		// ivars cvars target varname varlabel is_interaction is_cont_interaction is_bivariate is_mock levels
		mata: fe2local(`g')
		if ("`target'"=="") continue

		* Rename, add label and chars
		rename __Z`g'__ `target'
		local label `varlabel'
		la var `target' "Fixed effects of `label' on `original_depvar'"
		char `target'[label] `label'
		char `target'[levels] `levels'

		* Substract mean, or divide by cvar (fixing division by zero errors)
		if ("`cvars'"!="" & !(`is_bivariate' & !`is_mock')) {
			char `target'[cvars] `cvars'
			qui replace `target' = cond(abs(`cvars')<epsfloat(), 0,  `target'/`cvars')
			// BUGBUG BUGBUG float(`target'/`cvars')) -> this makes them have the same FE but loses precision!
		}
		else {
			qui su `target' `weightexp', mean
			qui replace `target' = `target' - r(mean)
			// BUGBUG BUGBUG -> WHAT WAS THE PROBLEM WITH THIS?
		}

		local keepvars `keepvars' `target'
	}

	cap drop __Z*__
	return local keepvars " `keepvars'" // the space prevents MVs
end

program define Stop
	cap mata: mata drop prev_numstep // Created at step 1
	cap mata: mata drop VERBOSE // Created before step 1
	cap mata: mata drop G // Num of absorbed FEs
	cap mata: mata drop FEs // Main Mata structure
	cap mata: mata drop betas // Temporary matrices used to store bi/multivariate regr coefs
	cap mata: mata drop varlist_cache // Hash table with the names of the precomputed residuals
	cap mata: mata drop avge_* // Drop AvgE structures
	cap mata: mata drop weightexp weightvar

	cap mata: mata drop clustervars
	cap mata: mata drop clustervars_ivars
	cap mata: mata drop clustervars_original

	if ("${hdfe_pwd}"!="") {
		qui cd "${hdfe_pwd}"
		global hdfe_pwd
	}

	* PARALLEL SPECIFIC CLEANUP
	cap mata: st_local("path", parallel_path)
	if ("`path'"!="") {
		mata: st_local("cores", strofreal(parallel_cores))
		assert "`cores'"!=""
		local path "`path'"
		cap erase `"`path'hdfe_mata.mo"'
		forv core=1/`cores' {
			cap erase `"`path'`core'_done.txt"'
			cap erase `"`path'`core'_ok.txt"'
			cap erase `"`path'`core'_error.txt"'
			cap erase `"`path'`core'_output.dta"'
			cap erase `"`path'`core'_log.log"'
		}
		cap rmdir `"`path'"'
		cap mata: mata drop parallel_cores
		cap mata: mata drop parallel_dta
		cap mata: mata drop parallel_vars
		cap mata: mata drop parallel_opt
		cap mata: mata drop parallel_path
	}
end

program define CheckCorrectOrder
	args step

	local numstep = ("`step'"=="start") + 2*("`step'"=="precompute") + ///
		3*("`step'"=="demean") + 4*("`step'"=="save")
	Assert (`numstep'>0), msg("hdfe: -`step'- is an invalid step")

	cap mata: st_local("prev_numstep", strofreal(prev_numstep))
	if (_rc) local prev_numstep 0

	Assert (`numstep'==`prev_numstep'+1) | (`numstep'==3 & `prev_numstep'==3) ///
		, msg("hdfe: expected step `=`prev_numstep'+1' instead of step 	`numstep'")
	mata: prev_numstep = `numstep'
	Debug, msg(_n as text "{title:Running -hdfe- step `numstep'/5 (`step')}") level(3)
end


// -------------------------------------------------------------
// Iteratively drop singletons for each absvar
// -------------------------------------------------------------
* This could be done iteratively, dropping singletons for each absvar until no progress is made.
* However, that would be extremely slow for a modest gain

program define DropSingletons, sortpreserve
syntax, num_absvars(integer)

	forv g=1/`num_absvars' {
		mata: fe2local(`g')
		* ivars cvars target varname varlabel is_interaction is_cont_interaction is_bivariate is_mock levels
		
		* It's either redudant (the second part of i.a##c.b) or tricky (a simple i.a#c.b) to discard singletons with cont. interactions
		local is_slope =  ("`cvars'"!="") & (!`is_bivariate' | `is_mock')
		if (`is_slope') continue

		local N_old = c(N)
		qui bys `ivars': drop if _N==1
		local N_new = c(N)
		local N_dropped = (`N_old' - `N_new')
		if (`N_dropped'>0) Debug, level(0) msg("(dropped `N_dropped' singleton observations for absvar " as result "`varlabel'" as text ")")
	}
end

