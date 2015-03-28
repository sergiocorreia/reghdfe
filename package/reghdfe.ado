*! reghdfe 2.0.294 27mar2015
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

program define reghdfe
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

	if replay() {
		if (`"`e(cmd)'"'!="reghdfe") error 301
		Replay `0'
	}
	else {
		* Estimate, and then clean up Mata in case of failure
		mata: st_global("reghdfe_pwd",pwd())
		Stop // clean leftovers for a possible [break]
		cap noi Estimate `0'
		if (_rc) {
			local rc = _rc
			Stop
			exit `rc'
		}
	}
end

* Note: Assert and Debug must go first

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
    local version "2.0.294 27mar2015"
    ereturn clear
    di as text "`version'"
    ereturn local version "`version'"

    di as text _n "Dependencies installed?"
    local dependencies ivreg2 avar tuples parallel
    foreach dependency of local dependencies {
    	cap findfile `dependency'.ado
    	if (_rc) {
    		di as text "{lalign 20:- `dependency'}" as error "not"
    	}
    	else {
    		di as text "{lalign 20:- `dependency'}" as result "yes"
    	}
    }

end



mata:
mata set matastrict on

// -------------------------------------------------------------------------------------------------
// Fix nonpositive VCV; called from Wrapper_mwc.ado 
// -------------------------------------------------------------------------------------------------
void function fix_psd(string scalar Vname) {
	real matrix V, U, lambda

	V = st_matrix(Vname)
	if (!issymmetric(V)) exit(error(505))
	symeigensystem(V, U=., lambda=.)
	st_local("eigenfix", "0")
	if (min(lambda)<0) {
		lambda = lambda :* (lambda :>= 0)
		// V = U * diag(lambda) * U'
		V = quadcross(U', lambda, U')
		st_local("eigenfix", "1")
	}
	st_replacematrix(Vname, V)
}

end



// -------------------------------------------------------------------------------------------------
// Transform data and run the regression
// -------------------------------------------------------------------------------------------------

program define Estimate, eclass

/* Notation of created variables
	__FE1__        		Fixed effect categories
	__Z1__         		Fixed effect coefficients (estimates)
	__clustervar1__		Categories for the clusters that had to be generated
	__W1__         		AvgE transformed variables (avg of depvar by category)
*/

// PART I - PREPARE DATASET FOR REGRESSION

* 1) Parse main options
	Parse `0' // save all arguments into locals (verbose>=3 shows them)
	local sets depvar indepvars endogvars instruments // depvar MUST be first

* 2) Parse identifiers (absorb variables, avge, clustervar)
	Start, absorb(`absorb') over(`over') avge(`avge') clustervars(`clustervars') weight(`weight') weightvar(`weightvar')
	* Note: In this step, it doesn't matter if the weight is FW or AW
	local N_hdfe = r(N_hdfe)
	local N_avge = r(N_avge)
	local RAW_N = c(N)
	local RAW_K = c(k)
	local absorb_keepvars = r(keepvars) // Vars used in hdfe,avge,cluster
	
	qui de, simple
	local old_mem = string(r(width) * r(N)  / 2^20, "%6.2f") // This is just for debugging; measured in MBs

* 3) Preserve
if ("`usecache'"!="") {
	local uid __uid__
	if ("`over'"!="") {
		gettoken ifword ifexp : if
		expr_query `ifexp'
		local vars_in_if = r(varnames)
		Assert `: list over in vars_in_if', msg("Error: since you are using over(`over'), you need to include {it:`over'}=={it:value} to your -if- condition.")
		
		cap local regex = regexm("`if'", "(^| )`over'==([0-9.e-]+)")
		Assert `regex', rc(0) msg("Warning: {it:`over'}=={it:value} not found in -if- (perhaps was abbreviated); e(over_value) and e(over_label) will not be stored.")
		cap local over_value = regexs(2)
		cap local over_label : label (__uid__) `over_value'
	}
}
else {
	tempvar uid
	local uid_type = cond(`RAW_N'>c(maxlong), "double", "long")
	gen `uid_type' `uid' = _n // Useful for later merges
	la var `uid' "[UID]" // So I can recognize it in -describe-
}

	if (`savingcache') {
		cap drop __uid__
		rename `uid' __uid__
		local uid __uid__
		local handshake = int(uniform()*1e8)
		char __uid__[handshake] `handshake'
		char __uid__[tolerance] `tolerance'
		char __uid__[maxiterations] `maxiterations'
		if ("`over'"!="") {
			local label : value label `over'
			label value __uid__ `label', nofix // Trick, attach label to __uid__
		}
	}

	preserve
	Debug, msg("(dataset preserved)") level(2)

* 4) Drop unused variables
	if ("`vceextra'"!="") local tsvars `panelvar' `timevar' // We need to keep these when using an autoco-robust VCE
	local exp "= `weightvar'"
	marksample touse, novar // Uses -if- , -in- ; -weight-? and -exp- ; can't drop any var until this
	keep `uid' `touse' `timevar' `panelvar' `absorb_keepvars' `basevars' `over' `weightvar' `tsvars'

* 5) Expand factor and time-series variables (this *must* happen before precompute is called!)
	local expandedvars
	foreach set of local sets {
		local varlist ``set''
		if ("`varlist'"=="") continue
		local original_`set' `varlist'
		* the -if- prevents creating dummies for categories that have been excluded
		ExpandFactorVariables `varlist' if `touse', setname(`set')
		local `set' "`r(varlist)'"
		local expandedvars `expandedvars' ``set''
	} 

* 6) Drop unused basevars and tsset vars (usually no longer needed)
	keep `uid' `touse' `absorb_keepvars' `expandedvars' `over' `weightvar' `tsvars'

* 7) Drop all observations with missing values (before creating the FE ids!)
	markout `touse' `expandedvars'
	markout `touse' `expandedvars' `absorb_keepvars'
	qui keep if `touse'
	if ("`dropsingletons'"!="") DropSingletons, num_absvars(`N_hdfe')
	Assert c(N)>0, rc(2000) msg("Empty sample, check for missing values or an always-false if statement")
	drop `touse'
	if ("`over'"!="" & `savingcache') qui levelsof `over', local(over_levels)

* 8) Fill Mata structures, create FE identifiers, avge vars and clustervars if needed
	Precompute, keep(`uid' `expandedvars' `tsvars') depvar("`depvar'") `excludeself' tsvars(`tsvars') over(`over')
	Debug, level(2) msg("(dataset compacted: observations " as result "`RAW_N' -> `c(N)'" as text " ; variables " as result "`RAW_K' -> `c(k)'" as text ")")
	local avgevars = cond("`avge'"=="", "", "__W*__")
	local vars `expandedvars' `avgevars'

	* qui compress `expandedvars' // will recast to -double- later on
	qui de, simple
	local new_mem = string(r(width) * r(N) / 2^20, "%6.2f")
	Debug, level(2) msg("(dataset compacted, c(memory): " as result "`old_mem'" as text "M -> " as result "`new_mem'" as text "M)")

* 9) Check that weights have acceptable values
if ("`weightvar'"!="") {
	local require_integer = ("`weight'"=="fweight")
	local num_type = cond(`require_integer', "integers", "reals")

	local basenote "weight -`weightvar'- can only contain strictly positive `num_type', but"
	qui cou if `weightvar'<0
	Assert (`r(N)'==0), msg("`basenote' `r(N)' negative values were found!")
	qui cou if `weightvar'==0
	Assert (`r(N)'==0), msg("`basenote' `r(N)' zero values were found!")
	qui cou if `weightvar'>=.
	Assert (`r(N)'==0), msg("`basenote' `r(N)' missing values were found!")
	if (`require_integer') {
		qui cou if mod(`weightvar',1)
		Assert (`r(N)'==0), msg("`basenote' `r(N)' non-integer values were found!")
	}
}

* 10) Save the statistics we need before transforming the variables
if (`savingcache') {
	cap drop __FE*__
	cap drop __clustervar*__
}
else {
	* Compute TSS of untransformed depvar
	local tmpweightexp = subinstr("`weightexp'", "[pweight=", "[aweight=", 1)
	qui su `depvar' `tmpweightexp' // BUGBUG: Is this correct?!
	local tss = r(Var)*(r(N)-1)
	assert `tss'<.

	if (`: list posof "first" in stages') {
		foreach var of varlist `endogvars' {
			qui su `var' `tmpweightexp' // BUGBUG: Is this correct?!
			local tss_`var' = r(Var)*(r(N)-1)
		}
	}

* 11) Calculate the degrees of freedom lost due to the FEs
	if ("`group'"!="") {
		tempfile groupdta
		local opt group(`group') groupdta(`groupdta') uid(`uid')
	}
	EstimateDoF, dofadjustments(`dofadjustments') `opt'
	local kk = r(kk) // FEs that were not found to be redundant (= total FEs - redundant FEs)
	local M = r(M) // FEs found to be redundant
	local saved_group = r(saved_group)
	local M_due_to_nested = r(M_due_to_nested)

	Assert `kk'<.
	Assert `M'>=0 & `M'<.
	assert inlist(`saved_group', 0, 1)

	forv g=1/`N_hdfe' {
		local M`g' = r(M`g')
		local K`g' = r(K`g')
		local M`g'_exact = r(M`g'_exact)
		local M`g'_nested = r(M`g'_nested)

		assert inlist(`M`g'_exact',0,1) // 1 or 0 whether M`g' was calculated exactly or not
		assert `M`g''<. & `K`g''<.
		assert `M`g''>=0 & `K`g''>=0
		assert inlist(r(drop`g'), 0, 1)

		* Drop IDs for the absorbed FEs (except if its the clustervar)
		* Useful b/c regr. w/cluster takes a lot of memory
		if (r(drop`g')==1) drop __FE`g'__
	}

	if (`num_clusters'>0) {
		mata: st_local("temp_clustervars", invtokens(clustervars))
		local vceoption : subinstr local vceoption "<CLUSTERVARS>" "`temp_clustervars'"
	}

}

* 12) Save untransformed data.
*	This allows us to:
*	i) do nested ftests for the FEs,
*	ii) recover the FEs, compute their correlations with xb, check that FE==1

	* We can avoid this if i) nested=check=0 ii) targets={} iii) fast=1
	mata: st_local("any_target_avge", strofreal(any(avge_target :!= "")) ) // saving avge?
	local any_target_hdfe 0 // saving hdfe?
	forv g=1/`N_hdfe' {
		mata: fe2local(`g')
		if (!`is_bivariate' | `is_mock') local hdfe_cvar`g' `cvars'
		// If it's the intercept part of the bivariate absorbed effect, don't add the cvar!
		local hdfe_target`g' `target'
		if ("`target'"!="") local any_target_hdfe 1
	}

	if (`fast') {
		if (`nested' | `check' | `any_target_hdfe' | `any_target_avge' | "`group'"!="") {
			Debug, msg(as text "(option {it:fast} not compatible with other options; disabled)") level(0)
			local fast 0
		}
		else {
			Debug, msg("(option {opt fast} specified; will not save e(sample) or compute correlations)")
		}
	}

	if (!`fast' | `cores'>1) {
		sort `uid'
		tempfile original_vars
		qui save "`original_vars'"
		if (`cores'>1) local parallel_opt `" filename("`original_vars'") uid(`uid') cores(`cores') "'
		Debug, msg("(untransformed dataset saved)") level(2)
	}

* 13) (optional) Compute R2/RSS to run nested Ftests on the FEs
	* a) Compute R2 of regression without FE, to build the joint FTest for all the FEs
	* b) Also, compute RSS of regressions with less FEs so we can run nested FTests on the FEs
	if ("`model'"=="ols" & !`savingcache') {
		qui _regress `vars' `weightexp', noheader notable
		local r2c = e(r2)

		if (`nested') {
			local rss0 = e(rss)
			local subZs
			forv g=1/`=`N_hdfe'-1' {
				Debug, msg("(computing nested model w/`g' FEs)")
				if (`cores'>1) {
					DemeanParallel, varlist(`vars') `maximize_options' num_fe(`g') self(reghdfe) `parallel_opt'
				}
				else {
					Demean, varlist(`vars') `maximize_options' num_fe(`g')	
				}

				qui _regress `vars' `weightexp', noheader notable
				local rss`g' = e(rss)
				qui use "`original_vars'", clear // Back to untransformed dataset
			}
		}
	}

	* Get normalized string of the absvars (i.e. turn -> i.turn)
	local original_absvars
	forv g=1/`N_hdfe' {
		mata: fe2local(`g')
		local original_absvars `original_absvars'  `varlabel'
	}

* Compute summary statistics for the all the regression variables
	if ("`stats'"!="") {
		local tabstat_weight : subinstr local weightexp "[pweight" "[aweight"
		qui tabstat `vars' `tabstat_weight' , stat(`stats') col(stat) save
		tempname statsmatrix
		matrix `statsmatrix' = r(StatTotal)
	}

* 14) Compute residuals for all variables including the AvgEs (overwrites vars!)
	qui ds `vars'
	local NUM_VARS : word count `r(varlist)'
	Debug, msg("(computing residuals for `NUM_VARS' variables)")
	Debug, msg(" - tolerance = `tolerance'")
	Debug, msg(" - max. iter = `maxiterations'")
	if ("`usecache'"=="") {
		if (`cores'>1) {
			DemeanParallel, varlist(`vars') `maximize_options' self(reghdfe) `parallel_opt'
		}
		else {
			Demean, varlist(`vars') `maximize_options'	
		}
	}
	else {
		Debug, msg("(using cache data)")
		drop `vars'
		local handshake_master : char __uid__[handshake]
		char __uid__[handshake]
		// An error in the merge most likely means different # of obs due to missing values in a group but not in other
		// try with if !missing(__uid__) // TODO: Auto-add this by default?
		// TODO: Make this fool-proof when using -over-
		if ("`over'"!="") local using using // This is dangerous
		sort __uid__ // The user may have changed the sort order of the master data
		qui merge 1:1 __uid__ using "`usecache'", keepusing(`vars') assert(match master `using') keep(master match) nolabel sorted
		qui cou if _merge!=3
		Assert r(N)==0, msg(as error "Error: the cache has `r(N)' less observations than the master data" _n ///
			as text " - This is possibly because, when created, it included variables that were missing in cases where the current ones are not.")
		qui drop if _merge!=3
		drop _merge

		local handshake_using : char __uid__[handshake]
		local tolerance_using : char __uid__[tolerance]
		local maxiterations_using : char __uid__[maxiterations]
		Assert (`handshake_master'==`handshake_using'), msg("using dataset does not have the same __uid__")
		Assert abs(`tolerance'-`tolerance_using')<epsdouble(), msg("using dataset not computed with the same tolerance (`tolerance_using')")
		Assert (`maxiterations'==`maxiterations_using'), msg("using dataset not computed with the same maxiterations (`maxiterations_using')")

		local absvar_master `original_absvars'
		local absvar_using : char __uid__[absvars_key]
		Assert ("`absvar_master'"=="`absvar_using'"), msg("using dataset not created with the same absvars")
		char __uid__[absvars_key]
	}

if (`savingcache') {
	Debug, msg("(saving cache and exiting)")
	char __uid__[absvars_key] `original_absvars'
	sort __uid__
	save "`savecache'", replace
	return clear
	ereturn clear
	ereturn local cmdline `"`cmdline'"'
	if ("`over_levels'"!="") ereturn local over_levels = "`over_levels'"
	exit
}

// PART II - REGRESSION

**** <<< START OF UGLY -stages- CODE
assert "`stages'"!=""
if ("`stages'"!="none") {
	Debug, level(2) msg(_n " {title:Stages to run}: " as result "`stages'" _n)
	local backup_fast `fast'
	local num_stages : word count `stages'
	local last_stage : word `num_stages' of `stages'
	assert "`last_stage'"=="iv"
	foreach vargroup in depvar indepvars endogvars instruments {
		local backup_`vargroup' ``vargroup''
		local backup_original_`vargroup' `original_`vargroup''
	}
	local backup_tss = `tss'
}

foreach stage of local stages {
local lhs_endogvars = cond("`stage'"=="first", "`backup_endogvars'", "<none>")

if ("`stage'"=="first") {
	local i_endogvar 0
}
else {
	local i_endogvar
}

foreach lhs_endogvar of local lhs_endogvars {
Assert inlist("`stage'", "none", "iv", "first", "ols", "reduced", "acid")

if ("`stage'"=="iv") {
	local tss = `backup_tss'
	local fast `backup_fast'
	local depvar `backup_depvar'
	local indepvars `backup_indepvars'
	local endogvars `backup_endogvars'
	local instruments `backup_instruments'
	local original_depvar `backup_original_depvar'
	local original_indepvars `backup_original_indepvars'
	local original_endogvars `backup_original_endogvars'
	local original_instruments `backup_original_instruments'
}
else if ("`stage'"=="ols") {
	local tss = `backup_tss'
	local fast 1
	local depvar `backup_depvar'
	local endogvars
	local indepvars `backup_indepvars' `backup_endogvars'
	local instruments
	local original_depvar `backup_original_depvar'
	local original_endogvars
	local original_indepvars `backup_original_indepvars' `backup_original_endogvars'
	local original_instruments
	local vcesuite avar
}
else if ("`stage'"=="reduced") {
	local tss = `backup_tss'
	local fast 1
	local depvar `backup_depvar'
	local indepvars `backup_indepvars' `backup_instruments'
	local endogvars
	local instruments
	local original_depvar `backup_original_depvar'
	local original_indepvars `backup_original_indepvars' `backup_original_instruments'
	local original_endogvars
	local original_instruments
	local vcesuite avar
}
else if ("`stage'"=="acid") {
	local tss = `backup_tss'
	local fast 1
	local depvar `backup_depvar'
	local indepvars `backup_indepvars' `backup_endogvars' `backup_instruments'
	local endogvars
	local instruments
	local original_depvar `backup_original_depvar'
	local original_indepvars `backup_original_indepvars' `backup_original_endogvars' `backup_original_instruments'
	local original_endogvars
	local original_instruments
	local vcesuite avar
}
else if ("`stage'"=="first") {
	local ++ i_endogvar
	local tss = `tss_`lhs_endogvar''
	local fast 1
	local depvar `lhs_endogvar'
	local indepvars `backup_indepvars' `backup_instruments'
	local endogvars
	local instruments
	local original_depvar : word `i_endogvar' of `backup_original_endogvars'
	local original_indepvars `backup_original_indepvars' `backup_original_endogvars' `backup_original_instruments'
	local original_endogvars
	local original_instruments
	local vcesuite avar
}
**** END OF UGLY -stages- CODE >>>> 

* Cleanup
	ereturn clear

* Add back constant
	if (`addconstant') {
		Debug, level(3) msg(_n "adding back constant to regression")
		AddConstant `depvar' `indepvars' `avgevars' `endogvars' `instruments'
	}

* Regress
	if ("`stage'"=="none") Debug, level(2) msg("(running regresion: `model'.`ivsuite')")
	local avge = cond(`N_avge'>0, "__W*__", "")
	local options
	local option_list ///
		depvar indepvars endogvars instruments avgevars ///
		original_depvar original_indepvars original_endogvars ///
		original_instruments original_absvars avge_targets ///
		vceoption vcetype vcesuite ///
		kk suboptions showraw vceunadjusted first weightexp ///
		addconstant /// tells -regress- to hide _cons
		estimator twicerobust // Whether to run or not two-step gmm
	foreach opt of local option_list {
		if ("``opt''"!="") local options `options' `opt'(``opt'')
	}

	* Five wrappers in total, two for iv (ivreg2, ivregress), three for ols (regress, avar, mwc)
	local wrapper "Wrapper_`subcmd'" // regress ivreg2 ivregress
	if ("`subcmd'"=="regress" & "`vcesuite'"=="avar") local wrapper "Wrapper_avar"
	if ("`subcmd'"=="regress" & "`vcesuite'"=="mwc") local wrapper "Wrapper_mwc"

	if (!inlist("`stage'","none", "iv")) local wrapper "Wrapper_avar" // Compatible with ivreg2
	Debug, level(3) msg(_n "call to wrapper:" _n as result "`wrapper', `options'")
	`wrapper', `options'
	
	Assert e(tss)<., msg("within tss is missing (wrapper=`wrapper')")
	
	local subpredict = e(predict) // used to recover the FEs

	if ("`weightvar'"!="") {
		qui su `weightvar', mean
		local sumweights = r(sum)
	}

// PART III - RECOVER FEs AND SAVE RESULTS 

if (`fast') {
	* Copy pasted from below
	Debug, level(3) msg("(avoiding -use- of temporary dataset)")
	tempname b
	matrix `b' = e(b)
	local backup_colnames : colnames `b'
	FixVarnames `backup_colnames'
	local newnames "`r(newnames)'"
	local prettynames "`r(prettynames)'"
	matrix colnames `b' = `newnames'
}
else {
	assert inlist("`stage'", "iv", "none")
* 1) Restore untransformed dataset
	qui use "`original_vars'", clear

* 2) Recover the FEs

	* Predict will get (e+d) from the equation y=xb+d+e
	tempvar resid_d
	if e(df_m)>0 {
		local score = cond("`model'"=="ols", "score", "resid")
		`subpredict' double `resid_d', `score' // Auto-selects the program based on the estimation method		
	}
	else {
		gen double `resid_d' = `depvar'
	}

	* If the eqn doesn't have a constant, we need to save the mean of the resid in order to add it when predicting xb
	if (!`addconstant') {
		su `resid_d', mean
		ereturn `hidden' scalar _cons = r(mean)
	}

	Debug, level(2) msg("(loaded untransformed variables, predicted residuals)")

	* Absorb the residuals to obtain the FEs (i.e. run a regression on just the resids)
	Debug, level(2) tic(31)
	Demean, varlist(`resid_d') `maximize_options' save_fe(1)
	Debug, level(2) toc(31) msg("mata:make_residual on final model took")
	drop `resid_d'

* 3) Compute corr(FE,xb) (do before rescaling by cvar or deleting)
	if ("`model'"=="ols") {
		tempvar xb
		_predict double `xb', xb // -predict- overwrites sreturn, use _predict if needed
		forv g=1/`N_hdfe' { 
			qui corr `xb' __Z`g'__
			local corr`g' = r(rho)
		}
		drop `xb'
	}

* 4) Replace tempnames in the coefs table
	* (e.g. __00001 -> L.somevar)
	* (this needs to be AFTER predict but before deleting FEs and AvgEs)
	tempname b
	matrix `b' = e(b)
	local backup_colnames : colnames `b'
	FixVarnames `backup_colnames'
	local newnames "`r(newnames)'"
	local prettynames "`r(prettynames)'"
	matrix colnames `b' = `newnames'

* 5) Save FEs w/proper name, format
	Save, original_depvar(`original_depvar')
	local keepvars `r(keepvars)'
	if ("`keepvars'"!="") format `fe_format' `keepvars'
* 6) Save AvgEs
	forv g=1/`N_avge' {
		local var __W`g'__
		local target : char `var'[target]
		if ("`target'"!="") {
			rename `var' `target'
			local avge_target`g' `target' // Used by -predict-
			local keepvars `keepvars' `target'
		}
	}

	if ("`keepvars'"!="") format `fe_format' `keepvars' // The format of depvar, saved by -Parse-

* 7) Save dataset with FEs and e(sample)
	keep `uid' `keepvars'
	tempfile output
	qui save "`output'"
} // fast

* 8) Restore original dataset and merge
	if (inlist("`stage'","none", "iv")) restore // Restore user-provided dataset (since -iv- comes at the end, that is done at that stage!)
	if (!`fast') {
		// `saved_group' was created by EstimateDoF.ado
		if (!`saved_group')  local groupdta
		SafeMerge, uid(`uid') file("`output'") groupdta("`groupdta'")
		*cap tsset, noquery // we changed -sortby- when we merged (even if we didn't really resort)
	}

// PART IV - ERETURN OUTPUT

	if (`c(version)'>=12) local hidden hidden // ereturn hidden requires v12+

* Ereturns common to all commands
	ereturn local cmd = "reghdfe"
	ereturn local subcmd = cond(inlist("`stage'", "none", "iv"), "`subcmd'", "regress")
	ereturn local cmdline `"`cmdline'"'
	ereturn local model = cond("`gmm2s'"=="", "`model'", "gmm2s")
	ereturn local model = cond("`cue'"=="", "`model'", "cue")
	ereturn local model = cond("`liml'"=="", "`model'", "liml")
	ereturn local dofadjustments = "`dofadjustments'"
	ereturn local title = "HDFE " + e(title)
	ereturn local title2 =  "Absorbing `N_hdfe' HDFE " + plural(`N_hdfe', "indicator")
	ereturn local predict = "reghdfe_p"
	ereturn local estat_cmd = "reghdfe_estat"
	ereturn local footnote = "reghdfe_footnote"
	ereturn local absvars = "`original_absvars'"
	ereturn local vcesuite = "`vcesuite'"
	ereturn local maximize_options = "`maximize_options'" // In option format; tolerance(..) etc.
	if ("`stage'"!="none") ereturn local iv_depvar = "`backup_original_depvar'"
	ereturn `hidden' local diopts = "`diopts'"
	if ("`over'"!="") {
		ereturn local over = "`over'"
		if ("`over_value'"!="") ereturn local over_value = "`over_value'"
		if ("`over_label'"!="") ereturn local over_label = "`over_label'"
		local fixed_absvars = e(absvars)
		local fixed_absvars : subinstr local fixed_absvars "i.`over'#" "", all
		local fixed_absvars : subinstr local fixed_absvars "i.`over'" "", all
		local fixed_absvars `fixed_absvars' // Trim
		ereturn local absvars = "`fixed_absvars'"
	}

	if ("`e(clustvar)'"!="") {
		mata: st_local("clustvar", invtokens(clustervars_original))
		* With kiefer/dkraay we add a time clustervar
		if ("`clustvar'"!="") ereturn local clustvar "`clustvar'"
		ereturn scalar N_clustervars = `num_clusters'
	}

	* Besides each cmd's naming style (e.g. exogr, exexog, etc.) keep one common one
	foreach cat in depvar indepvars endogvars instruments {
		local vars ``cat''
		if ("`vars'"=="") continue
		ereturn local `cat' "`original_`cat''"
	}
	ereturn local avgevars "`avge'" // bugbug?

	ereturn `hidden' local subpredict = "`subpredict'"
	ereturn `hidden' local prettynames "`prettynames'"
	forv g=1/`N_avge' {
		ereturn `hidden' local avge_target`g' "`avge_target`g''" // Used by -predict-
	}

	* Stata uses e(vcetype) for the SE column headers
	* In the default option, leave it empty.
	* In the cluster and robust options, set it as "Robust"
	ereturn local vcetype = proper("`vcetype'") //
	if (e(vcetype)=="Cluster") ereturn local vcetype = "Robust"
	if (e(vcetype)=="Unadjusted") ereturn local vcetype
	if ("`e(vce)'"=="." | "`e(vce)'"=="") ereturn local vce = "`vcetype'" // +-+-
	Assert inlist("`e(vcetype)'", "", "Robust", "Jackknife", "Bootstrap")
	
	* Clear results that are wrong
	ereturn local ll
	ereturn local ll_0

	ereturn scalar N_hdfe = `N_hdfe'
	if ("`N_avge'"!="") ereturn scalar N_avge = `N_avge'

* Absorbed-specific returns
	ereturn scalar mobility = `M'
	ereturn scalar df_a = `kk'
	forv g=1/`N_hdfe' {
		ereturn scalar M`g' = `M`g''
		ereturn scalar K`g' = `K`g''
		ereturn `hidden' scalar M`g'_exact = `M`g'_exact' // 1 or 0 whether M`g' was calculated exactly or not
		ereturn `hidden' local corr`g' = "`corr`g''" //  cond("`corr`g''"=="", ., "`corr`g''")
		ereturn `hidden' local hdfe_target`g' = "`hdfe_target`g''"
		ereturn `hidden' local hdfe_cvar`g' = "`hdfe_cvar`g''"
		ereturn `hidden' scalar M`g'_nested = `M`g'_nested'
	}

	Assert e(df_r)<. , msg("e(df_r) is missing")
	ereturn scalar r2_within = 1 - e(rss) / e(tss)
	ereturn scalar tss = `tss'
	ereturn scalar mss = e(tss) - e(rss)
	ereturn scalar r2 = 1 - e(rss) / `tss'

	* ivreg2 uses e(r2c) and e(r2u) for centered/uncetered R2; overwrite first and discard second
	if (e(r2c)!=.) {
		ereturn scalar r2c = e(r2)
		ereturn scalar r2u = .
	}

	* Computing Adj R2 with custered SEs is tricky because it doesn't use the adjusted inputs:
	* 1) It uses N instead of N_clust
	* 2) For the DoFs, it uses N - Parameters instead of N_clust-1
	* 3) Further, to compute the parameters, it includes those nested within clusters
	
	* Note that this adjustment is NOT PERFECT because we won't compute the mobility groups just for improving the r2a
	* (when a FE is nested within a cluster, we don't need to compute mobilty groups; but to get the same R2a as other estimators we may want to do it)
	* Instead, you can set by hand the dof() argument and remove -cluster- from the list

	if ("`model'"=="ols" & `num_clusters'>0) Assert e(unclustered_df_r)<., msg("wtf-`vcesuite'")
	local used_df_r = cond(e(unclustered_df_r)<., e(unclustered_df_r), e(df_r)) - `M_due_to_nested'
	ereturn scalar r2_a = 1 - (e(rss)/`used_df_r') / (`tss' / (e(N)-1) )

	ereturn scalar rmse = sqrt( e(rss) / `used_df_r' )

	if (e(N_clust)<.) Assert e(df_r) == e(N_clust) - 1, msg("Error, `wrapper' should have made sure that N_clust-1==df_r")
	*if (e(N_clust)<.) ereturn scalar df_r = e(N_clust) - 1

	if ("`weightvar'"!="") ereturn scalar sumweights = `sumweights'

	if ("`model'"=="ols" & inlist("`vcetype'", "unadjusted", "ols")) {
		ereturn scalar F_absorb = (e(r2)-`r2c') / (1-e(r2)) * e(df_r) / `kk'
		if (`nested') {
			local rss`N_hdfe' = e(rss)
			local temp_dof = e(N) - 1 - e(df_m) // What if there are absorbed collinear with the other RHS vars?
			local j 0
			ereturn `hidden' scalar rss0 = `rss0'
			forv g=1/`N_hdfe' {
				local temp_dof = `temp_dof' - e(K`g') + e(M`g')
				*di in red "g=`g' RSS=`rss`g'' and was `rss`j''.  dof=`temp_dof'"
				ereturn `hidden' scalar rss`g' = `rss`g''
				ereturn `hidden' scalar df_a`g' = e(K`g') - e(M`g')
				ereturn scalar F_absorb`g' = (`rss`j''-`rss`g'') / `rss`g'' * `temp_dof' / e(df_a`g')
				ereturn `hidden' scalar df_r`g' = `temp_dof'
				local j `g'
			}   
		}
	}

	// There is a big assumption here, that the number of other parameters does not increase asymptotically
	// BUGBUG: We should allow the option to indicate what parameters do increase asympt.
	// BUGBUG; xtreg does this: est scalar df_r = min(`df_r':=N-1-K, `df_cl') why was that?

	if ("`savefirst'"!="") ereturn `hidden' scalar savefirst = `savefirst'

	* We have to replace -unadjusted- or else subsequent calls to -suest- will fail
	Subtitle `vceoption' // will set title2, etc. Run after e(bw) and all the others are set!
	if (e(vce)=="unadjusted") ereturn local vce = "ols"

	if ("`stages'"!="none") {
		ereturn local stage = "`stage'"
		ereturn `hidden' local stages = "`stages'"
	}

* Show table and clean up
	ereturn repost b=`b', rename // why here???

	if ("`stage'"!="none") Debug, level(0) msg(_n "{title:Stage: `stage'}" _n)
	if ("`lhs_endogvar'"!="<none>") Debug, level(0) msg("{title:Endogvar: `lhs_endogvar'}")
	Replay
	Attach, notes(`notes') statsmatrix(`statsmatrix') summarize_quietly(`summarize_quietly')

*** <<<< LAST PART OF UGLY STAGE <<<<	
if (!inlist("`stage'","none", "iv")) {
	local estimate_name reghdfe_`stage'`i_endogvar'
	local stored_estimates `stored_estimates' `estimate_name'
	local cmd estimates store `estimate_name', nocopy
	Debug, level(2) msg(" - Storing estimate: `cmd'")
	`cmd'
}
else if ("`stage'"=="iv") {
	* On the last stage, save list of all stored estimates
	assert "`stored_estimates'"!=""
	ereturn `hidden' local stored_estimates = "`stored_estimates'"
}

} // lhs_endogvar
} // stage
*** >>>> LAST PART OF UGLY STAGE >>>>

	Stop

end

// -------------------------------------------------------------------------------------------------

* The idea of this program is to keep the sort order when doing the merges

program define SafeMerge, eclass sortpreserve
syntax, uid(varname numeric) file(string) [groupdta(string)]
	* Merging gives us e(sample) and the FEs / AvgEs
	tempvar merge
	merge 1:1 `uid' using "`file'", assert(master match) nolabel nonotes noreport gen(`merge')
	
	* Add e(sample) from _merge
	tempvar sample
	gen byte `sample' = (`merge'==3)
	la var `sample' "[HDFE Sample]"
	ereturn repost , esample(`sample')
	drop `merge'

	* Add mobility group
	if ("`groupdta'"!="") merge 1:1 `uid' using "`groupdta'", assert(master match) nogen nolabel nonotes noreport sorted
end

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


// -------------------------------------------------------------
// Parsing and basic sanity checks for REGHDFE.ado
// -------------------------------------------------------------
// depvar: dependent variable
// indepvars: included exogenous regressors
// endogvars: included endogenous regressors
// instruments: excluded exogenous regressors

program define Parse

* Remove extra spacing from cmdline (just for aesthetics, run before syntax)
	cap syntax anything(name=indepvars) [if] [in] [fweight aweight pweight/] , SAVEcache(string) [*]
	local savingcache = (`=_rc'==0)

if (`savingcache') {

	* Disable these options
	local fast
	local nested

	syntax anything(name=indepvars) [if] [in] [fweight aweight pweight/] , ///
		Absorb(string) SAVEcache(string) ///
		[Verbose(integer 0) CHECK TOLerance(real 1e-7) MAXITerations(integer 10000) noACCELerate ///
		bad_loop_threshold(integer 1) stuck_threshold(real 5e-3) pause_length(integer 20) ///
		accel_freq(integer 3) accel_start(integer 6) /// Advanced optimization options
		CORES(integer 1) OVER(varname numeric) ///
		DROPSIngletons]

	cap conf file "`savecache'.dta"
	if (`=_rc'!=0) {
		cap conf new file "`savecache'.dta"
		Assert (`=_rc'==0), msg("reghdfe will not be able to save `savecache'.dta")
	}

}
else {
	mata: st_local("cmdline", stritrim(`"reghdfe `0'"') )
	ereturn clear // Clear previous results and drops e(sample)
	syntax anything(id="varlist" name=0 equalok) [if] [in] ///
		[fweight aweight pweight/] , ///
		Absorb(string) ///
		[VCE(string)] ///
		[DOFadjustments(string) GROUP(name)] ///
		[avge(string) EXCLUDESELF] ///
		[Verbose(integer 0) CHECK NESTED FAST] ///
		[TOLerance(real 1e-7) MAXITerations(integer 10000) noACCELerate] ///
		[IVsuite(string) SAVEFIRST FIRST SHOWRAW] /// ESTimator(string)
		[VCEUNADJUSTED] /// Option when running gmm2s with ivregress
		[SMALL Hascons TSSCONS] /// ignored options
		[kiefer] /// excluded
		[SUBOPTions(string)] /// Options to be passed to the estimation command (e.g . to regress)
		[bad_loop_threshold(integer 1) stuck_threshold(real 5e-3) pause_length(integer 20) accel_freq(integer 3) accel_start(integer 6)] /// Advanced optimization options
		[CORES(integer 1)] [USEcache(string)] [OVER(varname numeric)] ///
		[NOTES(string)] /// NOTES(key=value ..)
		[STAGEs(string)] ///
		[noCONstant] /// Disable adding back the intercept (mandatory with -ivreg2-)
		[DROPSIngletons] ///
		[ESTimator(string)] /// GMM2s CUE LIML
		[*] // For display options ; and SUmmarize(stats)
}

* Weight
* We'll have -weight- (fweight|aweight|pweight), -weightvar-, -exp-, and -weightexp-
	if ("`weight'"!="") {
		local weightvar `exp'
		conf var `weightvar' // just allow simple weights
		local weightexp [`weight'=`weightvar']
		local backupweight `weight'
	}

* Cache options
	if ("`usecache'"!="") {
		conf file "`usecache'.dta"
		conf var __uid__
		Assert ("`avge'"==""), msg("option -avge- not allowed with -usecache-")
		Assert ("`avge'"==""), msg("option -nested- not allowed with -usecache-")
	}

* Save locals that will be overwritten by later calls to -syntax-
	local ifopt `if'
	local inopt `in'

* Parse summarize(..)
	local default_stats mean min max
	ParseImplicit, opt(SUmmarize) default(`default_stats') input(`options') syntax([namelist(name=stats)] , [QUIetly]) inject(stats quietly)
	local summarize_quietly = ("`quietly'"!="")
	if ("`stats'"=="" & "`quietly'"!="") local stats `default_stats'

* Coef Table Options
if (!`savingcache') {
	_get_diopts diopts options, `options'
	Assert `"`options'"'=="", msg(`"invalid options: `options'"')
	if ("`hascons'`tsscons'"!="") di in ye "(option `hascons'`tsscons' ignored)"
}

* Over
	if ("`over'"!="") {
		unab over : `over', max(1)
		Assert ("`usecache'"!="" | "`savecache'"!=""), msg("-over- needs to be used together with either -usecache- or -savecache-")
	}

* Verbose
	assert inlist(`verbose', 0, 1, 2, 3, 4) // 3 and 4 are developer options
	mata: VERBOSE = `verbose' // Ugly hack to avoid using a -global-

* Show raw output of called subcommand (e.g. ivreg2)
	local showraw = ("`showraw'"!="")
	
* If true, will use wmatrix(...) vce(unadjusted) instead of the default of setting vce contents equal to wmatrix
* This basically undoes the extra adjustment that ivregress does, so it's comparable with ivreg2
*
* Note: Cannot match exactly the -ivregress- results without vceunadjusted (see test-gmm.do)
* Thus, I will set this to true ALWAYS
	local vceunadjusted = 1 // ("`vceunadjusted'"!="")

* tsset variables, if any
	cap conf var `_dta[_TStvar]'
	if (!_rc) local timevar `_dta[_TStvar]'
	cap conf var `_dta[_TSpanel]'
	if (!_rc) local panelvar `_dta[_TSpanel]'

* Model settings
if (!`savingcache') {

	// local model = cond(strpos(`"`0'"', " ("), "iv", "ols") // Fails with long strs in stata 12<
	local model ols
	foreach _ of local 0 {
		if (substr(`"`_'"', 1, 1)=="(") {
			local model iv
			continue, break
		}
	}
	
	* Estimator
	if ("`estimator'"!="") {
		Assert "`model'"=="iv", msg("reghdfe error: estimator() requires an instrumental-variable regression")
	
		if (substr("`estimator'", 1, 3)=="gmm") local estimator gmm2s

		Assert inlist("`estimator'", "2sls", "gmm2s", "liml", "cue"), msg("reghdfe error: invalid estimator `estimator'")
		if (inlist("`estimator'", "cue")) Assert "`ivsuite'"!="ivregress", msg("reghdfe error: estimator `estimator' only available with the ivreg2 command, you selected ivregress")
	}
	else {
		local estimator 2sls
	}

	* For this, _iv_parse would have been useful, although I don't want to do factor expansions when parsing
	if ("`model'"=="iv") {

		* get part before parentheses
		local wrongparens 1
		while (`wrongparens') {
			gettoken tmp 0 : 0 ,p("(")
			local left `left'`tmp'
			* Avoid matching the parens of e.g. L(-1/2) and L.(var1 var2)
			* Using Mata to avoid regexm() and trim() space limitations
			mata: st_local("tmp1", subinstr("`0'", " ", "") ) // wrong parens if ( and then a number
			mata: st_local("tmp2", substr(strtrim("`left'"), -1) ) // wrong parens if dot
			local wrongparens = regexm("`tmp1'", "^\([0-9-]") | ("`tmp2'"==".")
			if (`wrongparens') {
				gettoken tmp 0 : 0 ,p(")")
				local left `left'`tmp'
			}
		}

		* get part in parentheses
		gettoken right 0 : 0 ,bind match(parens)
		Assert trim(`"`0'"')=="" , msg("error: remaining argument: `0'")

		* now parse part in parentheses
		gettoken endogvars instruments : right ,p("=")
		gettoken equalsign instruments : instruments ,p("=")

		Assert "`endogvars'"!="", msg("iv: endogvars required")
		local 0 `endogvars'
		syntax varlist(fv ts numeric)
		local endogvars `varlist'

		Assert "`instruments'"!="", msg("iv: instruments required")
		local 0 `instruments'
		syntax varlist(fv ts numeric)
		local instruments `varlist'
		
		local 0 `left' // So OLS part can handle it
		Assert "`endogvars'`instruments'"!=""
		
		if ("`ivsuite'"=="") local ivsuite ivreg2
		Assert inlist("`ivsuite'","ivreg2","ivregress") , msg("error: wrong IV routine (`ivsuite'), valid options are -ivreg2- and -ivregress-")
		cap findfile `ivsuite'.ado
		Assert !_rc , msg("error: -`ivsuite'- not installed, please run {stata ssc install `ivsuite'} or change the option -ivsuite-")
	}

* OLS varlist
	syntax varlist(fv ts numeric)
	gettoken depvar indepvars : varlist
	_fv_check_depvar `depvar'

* Extract format of depvar so we can format FEs like this
	fvrevar `depvar', list
	local fe_format : format `r(varlist)' // The format of the FEs and AvgEs that will be saved

* Variables shouldn't be repeated
* This is not perfect (e.g. doesn't deal with "x1-x10") but still helpful
	local allvars `depvar' `indepvars' `endogvars' `instruments'
	local dupvars : list dups allvars
	Assert "`dupvars'"=="", msg("error: there are repeated variables: <`dupvars'>")

	Debug, msg(_n " {title:REGHDFE} Verbose level = `verbose'")
	*Debug, msg("{hline 64}")

* Stages
	assert "`model'"!="" // just to be sure this goes after `model' is set
	local iv_stage iv
	local stages : list stages - iv_stage
	local valid_stages ols first acid reduced
	local wrong_stages : list stages - valid_stages
	Assert "`wrong_stages'"=="", msg("Error, invalid stages(): `wrong_stages'")
	if ("`stages'"!="") {
		Assert "`model'"=="iv", msg("Error, stages() only valid with an IV regression")
		local stages `stages' `iv_stage' // Put -iv- *last* (so it does the -restore-; note that we don't need it first to trim MVs b/c that's done earlier)
		Assert "`avge'"=="", msg("Error, avge not allowed with stages()")
	}
	else {
		local stages none // So we can loop over stages
	}

* Add back constants (place this *after* we define `model')
	local addconstant = ("`constant'"!="noconstant") & !("`model'"=="iv" & "`ivsuite'"=="ivreg2") // also see below
	if (`addconstant' & "`over'"!="") {
		local addconstant 0
		Debug, level(0) msg("Constant will not be reported due to over(); use option -summarize()- or run the command -estat summ- to obtain the summary stats")
	}

* Parse VCE options:
	
	* Note: bw=1 *usually* means just do HC instead of HAC
	* BUGBUG: It is not correct to ignore the case with "bw(1) kernel(Truncated)"
	* but it's too messy to add -if-s everywhere just for this rare case (see also Mark Schaffer's email)

	local 0 `vce'
	syntax [anything(id="VCE type")] , [bw(integer 1)] [KERnel(string)] [dkraay(integer 1)] [kiefer] ///
		[suite(string) TWICErobust]
	if ("`anything'"=="") local anything unadjusted
	Assert `bw'>0, msg("VCE bandwidth must be a positive integer")
	gettoken vcetype clustervars : anything
	* Expand variable abbreviations; but this adds unwanted i. prefixes
	if ("`clustervars'"!="") {
		fvunab clustervars : `clustervars'
		local clustervars : subinstr local clustervars "i." "", all
	}

	* vcetype abbreviations:
	if (substr("`vcetype'",1,3)=="ols") local vcetype unadjusted
	if (substr("`vcetype'",1,2)=="un") local vcetype unadjusted
	if (substr("`vcetype'",1,1)=="r") local vcetype robust
	if (substr("`vcetype'",1,2)=="cl") local vcetype cluster
	if ("`vcetype'"=="conventional") local vcetype unadjusted // Conventional is the name given in e.g. xtreg
	Assert strpos("`vcetype'",",")==0, msg("Unexpected contents of VCE: <`vcetype'> has a comma")

	* Sanity checks on vcetype
	if ("`vcetype'"=="" & "`backupweight'"=="pweight") local vcetype robust
	Assert !("`vcetype'"=="unadjusted" & "`backupweight'"=="pweight"), msg("pweights do not work with unadjusted errors, use a different vce()")
	if ("`vcetype'"=="") local vcetype unadjusted
	Assert inlist("`vcetype'", "unadjusted", "robust", "cluster"), msg("VCE type not supported: `vcetype'")

	* Cluster vars
	local num_clusters : word count `clustervars'
	Assert inlist( (`num_clusters'>0) + ("`vcetype'"=="cluster") , 0 , 2), msg("Can't specify cluster without clustervars and viceversa") // XOR

	* VCE Suite
	local vcesuite `suite'
	if ("`vcesuite'"=="") local vcesuite default
	if ("`vcesuite'"=="default") {
		if (`bw'>1 | `dkraay'>1 | "`kiefer'"!="" | "`kernel'"!="") {
			local vcesuite avar
		}
		else if (`num_clusters'>1) {
			local vcesuite mwc
		}
	}

	Assert inlist("`vcesuite'", "default", "mwc", "avar"), msg("Wrong vce suite: `vcesuite'")
	if (inlist("`vcesuite'", "avar", "mwc")) local addconstant 0 // The constant messes up the VCV

	if ("`vcesuite'"=="mwc") {
		cap findfile tuples.ado
		Assert !_rc , msg("error: -tuples- not installed, please run {stata ssc install tuples} to estimate multi-way clusters.")
	}
	
	if ("`vcesuite'"=="avar" | "`stages'"!="none") {
		cap findfile avar.ado // We use -avar- as default with stages (on the non-iv stages)
		Assert !_rc , msg("error: -avar- not installed, please run {stata ssc install avar} or change the option -vcesuite-")
	}

	* Some combinations are not coded
	Assert !("`ivsuite'"=="ivregress" & (`num_clusters'>1 | `bw'>1 | `dkraay'>1 | "`kiefer'"!="" | "`kernel'"!="") ), msg("option vce(`vce') incompatible with ivregress")
	Assert !("`ivsuite'"=="ivreg2" & (`num_clusters'>2) ), msg("ivreg2 doesn't allow more than two cluster variables")
	Assert !("`model'"=="ols" & "`vcesuite'"=="avar" & (`num_clusters'>2) ), msg("avar doesn't allow more than two cluster variables")
	Assert !("`model'"=="ols" & "`vcesuite'"=="default" & (`bw'>1 | `dkraay'>1 | "`kiefer'"!="" | "`kernel'"!="") ), msg("to use those vce options you need to use -avar- as the vce suite")
	if (`num_clusters'>0) local temp_clustervars " <CLUSTERVARS>"
	if (`bw'==1 & `dkraay'==1 & "`kernel'"!="") local kernel // No point in setting kernel here 
	if (`bw'>1 | "`kernel'"!="") local vceextra `vceextra' bw(`bw') 
	if (`dkraay'>1) local vceextra `vceextra' dkraay(`dkraay') 
	if ("`kiefer'"!="") local vceextra `vceextra' kiefer 
	if ("`kernel'"!="") local vceextra `vceextra' kernel(`kernel')
	if ("`vceextra'"!="") local vceextra , `vceextra'
	local vceoption "`vcetype'`temp_clustervars'`vceextra'" // this excludes "vce(", only has the contents

* DoF Adjustments
	if ("`dofadjustments'"=="") local dofadjustments all
	local 0 , `dofadjustments'
	syntax, [ALL NONE] [PAIRwise FIRSTpair] [CLusters] [CONTinuous]
	opts_exclusive "`all' `none'" dofadjustments
	opts_exclusive "`pairwise' `firstpair'" dofadjustments
	if ("`none'"!="") {
		Assert "`pairwise'`firstpair'`clusters'`continuous'"=="", msg("option {bf:dofadjustments()} invalid; {bf:none} not allowed with other alternatives")
		local dofadjustments
	}
	if ("`all'"!="") {
		Assert "`pairwise'`firstpair'`clusters'`continuous'"=="", msg("option {bf:dofadjustments()} invalid; {bf:all} not allowed with other alternatives")
		local dofadjustments pairwise clusters continuous
	}
	else {
		local dofadjustments `pairwise' `firstpair' `clusters' `continuous'
	}

* Mobility groups
	if ("`group'"!="") conf new var `group'

* IV options
	if ("`small'"!="") di in ye "(note: reghdfe will always use the option -small-, no need to specify it)"

	if ("`model'"=="iv") {
		local savefirst = ("`savefirst'"!="")
		local first = ("`first'"!="")
		if (`savefirst') Assert `first', msg("Option -savefirst- requires -first-")
	}

} // End of !`savingcache'

* Add constant if -dropsingletons-; this applies for both savingcache and normal estimation
	if ("`addconstant'"=="1" & "`dropsingletons'"!="") {
		local addconstant 0
		Debug, level(0) msg("(constant will not be reported because -dropsingletons- changes the reported constant)")
	} 

* Optimization
	if (`maxiterations'==0) local maxiterations 1e8
	Assert (`maxiterations'>0)
	local accelerate = cond("`accelerate'"!="", 0, 1) // 1=Yes
	local check = cond("`check'"!="", 1, 0) // 1=Yes
	local fast = cond("`fast'"!="", 1, 0) // 1=Yes
	local tolerance = strofreal(`tolerance', "%9.1e") // Purely esthetic
	Assert `cores'<=32 & `cores'>0 , msg("At most 32 cores supported")
	if (`cores'>1) {
		cap findfile parallel.ado
		Assert !_rc , msg("error: -parallel- not installed, please run {stata ssc install parallel}")
	}
	local opt_list tolerance maxiterations check accelerate ///
		bad_loop_threshold stuck_threshold pause_length accel_freq accel_start
	foreach opt of local opt_list {
		if ("``opt''"!="") local maximize_options `maximize_options' `opt'(``opt'')
	}

* Varnames underlying tsvars and fvvars (e.g. i.foo L(1/3).bar -> foo bar)
	foreach vars in depvar indepvars endogvars instruments {
		if ("``vars''"!="") {
			fvrevar ``vars'' , list
			local basevars `basevars' `r(varlist)'
		}
	}

if (!`savingcache') {
* Nested
	local nested = cond("`nested'"!="", 1, 0) // 1=Yes
	if (`nested' & !("`model'"=="ols" & "`vcetype'"=="unadjusted") ) {
		Debug, level(0) msg("(option nested ignored, only works with OLS and conventional/unadjusted VCE)") color("error")
	}

* How can we do the same regression from a standard stata command?
* (useful for benchmarking and testing correctness of results)
	local subcmd = cond("`model'"=="ols" ,"regress", "`ivsuite'")

* _fv_check_depvar overwrites the local -weight-
	local weight `backupweight'
	Assert inlist( ("`weight'"!="") + ("`weightvar'"!="") + ("`weightexp'"!="") , 0 , 3 ) , msg("not all 3 weight locals are set")

* Return values
	local names cmdline diopts model ///
		ivsuite showraw vceunadjusted ///
		depvar indepvars endogvars instruments savefirst first ///
		vceoption vcetype vcesuite vceextra num_clusters clustervars bw kernel dkraay /// vceextra
		dofadjustments ///
		if in group check fast nested fe_format ///
		tolerance maxiterations accelerate maximize_options ///
		subcmd suboptions ///
		absorb avge excludeself ///
		timevar panelvar basevars ///
		addconstant ///
		weight weightvar exp weightexp /// type of weight (fw,aw,pw), weight var., and full expr. ([fw=n])
		cores savingcache usecache over ///
		stats summarize_quietly notes stages ///
		dropsingletons estimator twicerobust
}

if (`savingcache') {
	local names maximize_options cores if in timevar panelvar indepvars basevars ///
		absorb savecache savingcache fast nested check over ///
		weight weightvar exp weightexp /// type of weight (fw,aw), weight var., and full expr. ([fw=n])
		tolerance maxiterations /// Here just used for -verbose- and cache handshake purposes
		dropsingletons
}

	local if `ifopt'
	local in `inopt'

	Debug, level(3) newline
	Debug, level(3) msg("Parsed options:")
	foreach name of local names {
		if (`"``name''"'!="") Debug, level(3) msg("  `name' = " as result `"``name''"')
		c_local `name' `"``name''"' // Inject values into caller (reghdfe.ado)
	}
	// Debug, level(3) newline
end

program define ParseImplicit
* Parse options in the form NAME|NAME(arguments)
* Inject: what locals to inject (depend on -syntax)
* Default: default value for implicit form
* XOR: opt is mandatory (one of the two versions)
	syntax, opt(name local) default(string) syntax(string asis) [input(string asis)] inject(namelist local) [XOR]

	* First see if the implicit version is possible
	local lower_opt = lower("`opt'")
	local 0 , `input'
	cap syntax, `opt' [*]
	if ("`xor'"=="") local capture capture
	local rc = _rc
	if (`rc') {
		`capture' syntax, `opt'(string asis) [*]
		if ("`capture'"!="" & _rc) exit
	}
	else {
		local `lower_opt' `default'
	}
	local 0 ``lower_opt''
	syntax `syntax'
	foreach loc of local inject {
		c_local `loc' ``loc''
	}
	c_local options `options'
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


//------------------------------------------------------------------------------
// Expand Factor Variables, interactions, and time-series vars
//------------------------------------------------------------------------------
// This basically wraps -fvrevar-, adds labels, and drops omitted/base

program define ExpandFactorVariables, rclass
syntax varlist(min=1 numeric fv ts) [if] [,setname(string)] [CACHE]

	local expanded_msg `"" - variable expansion for `setname': " as result "`varlist'" as text " ->""'

	* It's (usually) a waste to add base and omitted categories
	* EG: if we use i.foreign#i.rep78 , several categories will be redundant, seen as e.g. "0b.foreign" in -char list-
	* We'll also exclude base categories that don't have the "bn" option (to have no base)

	* However, if there is a cont. interaction, then we DO want the base categories!

	* Loop for each var and then expand them into i.var -> 1.var.. and loop
	* Why two loops? B/c I want to save each var expansion to allow for a cache

	if ("`cache'"!="") mata: varlist_cache = asarray_create()

	local newvarlist
	* I can't do a simple foreach!
	* Because a factor expression could be io(3 4).rep78
	* and foreach would split the parens in half
	while (1) {
	gettoken fvvar varlist : varlist, bind
	if ("`fvvar'"=="") continue, break

		fvrevar `fvvar' `if' // , stub(__V__) // stub doesn't work in 11.2
		local contents

		foreach var of varlist `r(varlist)' {
			
			* Get readable varname
			local fvchar : char `var'[fvrevar]
			local tschar : char `var'[tsrevar]
			local name `fvchar'`tschar'
			local color input
			if ("`name'"=="") {
				local name `var'
				local color result
			}
			char `var'[name] `name'
			la var `var' "[Tempvar] `name'"

			* See if the factor can be dropped safely
			if (substr("`var'", 1, 2)=="__") {
				local color result
				local parts : subinstr local fvchar "#" " ", all
				local continteraction = strpos("`fvchar'", "c.")>0
				foreach part of local parts {
					*di as error "part=<`part'> cont=`continteraction' all=<`fvchar'>"
					* "^[0-9]+b\." -> "b.*\."
					if (regexm("`part'", "b.*\.") & !`continteraction') | regexm("`part'", "o.*\.") {
						local color error	
					}
				}
				if ("`color'"=="error") {
					local color result
				}


				* Need to rename it, or else it gets dropped since its a tempvar
				if ("`color'"!="error") {
					local newvarbase : subinstr local name "." "__", all // pray that no variable has three _
					local newvarbase : subinstr local newvarbase "#" "_X_", all // idem
					local newvarbase : permname __`newvarbase', length(30)

					* In what cases will just using newvarbase fail???
					local i // Empty
					while (1) {
						local newvar "`newvarbase'`i'"
					
						if ("`i'"=="") {
							local i 1
						}
						else {
							local ++i
						}

						Assert `i'<1000, msg("Couldn't create tempvar for `var' (`name')")
						cap conf new var `newvar', exact
						if _rc==0 {
							continue, break
						}
					}

					rename `var' `newvar'
					local var `newvar'
				}
			}

			* Save contents of the expansion for optional -cache-			
			if ("`color'"!="error") {
				local contents `contents' `var'
			}
			
			* Set debug message
			local expanded_msg `"`expanded_msg' as `color' " `name'" as text " (`var')""'
		}

		if ("`cache'"!="") mata: asarray(varlist_cache, "`fvvar'", "`contents'")
		Assert "`contents'"!="", msg("error: variable -`fvvar'- in varlist -`varlist'- in category -`setname'- is  empty after factor expansion")
		local newvarlist `newvarlist' `contents'
	}

	* Yellow=Already existed, White=Created, Red=NotCreated (omitted or base)
	Debug, level(3) msg(`expanded_msg')
	return local varlist "`newvarlist'"
end


//------------------------------------------------------------------------------
// Name tempvars into e.g. L.x i1.y i2.y AvgE:z , etc.
//------------------------------------------------------------------------------

program define FixVarnames, rclass
local vars `0'

	foreach var of local vars {
		local newname
		local pretyname

		* -var- can be <o.__W1__>
		if ("`var'"=="_cons") {
			local newname `var'
			local prettyname `var'
		}
		else {
			fvrevar `var', list
			local basevar "`r(varlist)'"
			local label : var label `basevar'
			local is_avge = regexm("`basevar'", "^__W[0-9]+__$")
			local is_temp = substr("`basevar'",1,2)=="__"
			local is_omitted = strpos("`var'", "o.")
			local prefix = cond(`is_omitted'>0, "o.", "")
			local name : char `basevar'[name]

			if (`is_avge') {
				local avge_str : char `basevar'[avge_equation]
				local name : char `basevar'[name]
				local prettyname `avge_str':`prefix'`name'

				local newname : char `basevar'[target]
				if ("`newname'"=="") local newname `var'
			}
			else if (`is_temp' & "`name'"!="") {
				local newname `prefix'`name'
				
				* Fix bug when the var is omitted:
				local bugmatch = regexm("`newname'", "^o\.([0-9]+)b?\.(.+)$")
				if (`bugmatch') {
					local newname = regexs(1) + "o." + regexs(2) // EG: 1o.var
				}

				local prettyname `newname'
			}
			else {
				local newname `var'
				local prettyname `newname'
			}
		}
		
		*di in red " var=<`var'> --> new=<`newname'> pretty=<`prettyname'>"
		Assert ("`newname'"!="" & "`prettyname'"!=""), ///
			msg("var=<`var'> --> new=<`newname'> pretty=<`prettyname'>")
		local newnames `newnames' `newname'
		local prettynames `prettynames' `prettyname'
	}

	local A : word count `vars'
	local B : word count `newnames'
	local C : word count `prettynames'
	Assert `A'==`B', msg("`A' vars but `B' newnames")
	Assert `A'==`C', msg("`A' vars but `C' newnames")
	
	***di as error "newnames=`newnames'"
	***di as error "prettynames=`prettynames'"

	return local newnames "`newnames'"
	return local prettynames "`prettynames'"
end

program define Wrapper_regress, eclass
	syntax , depvar(varname) [indepvars(varlist) avgevars(varlist)] ///
		original_absvars(string) original_depvar(string) [original_indepvars(string) avge_targets(string)] ///
		vceoption(string asis) ///
		kk(integer) ///
		[weightexp(string)] ///
		addconstant(integer) ///
		[SUBOPTions(string)] [*] // [*] are ignored!
	
	if ("`options'"!="") Debug, level(3) msg("(ignored options: `options')")
	mata: st_local("vars", strtrim(stritrim( "`depvar' `indepvars' `avgevars'" )) ) // Just for esthetic purposes
	if (`c(version)'>=12) local hidden hidden

* Convert -vceoption- to what -regress- expects
	gettoken vcetype clustervars : vceoption
	local clustervars `clustervars' // Trim
	local vceoption : subinstr local vceoption "unadjusted" "ols"
	local vceoption "vce(`vceoption')"

* Hide constant
	if (!`addconstant') {
		local nocons noconstant
		local kk = `kk' + 1
	}

* Note: the dof() option of regress is *useless* with robust errors,
* and overriding e(df_r) is also useless because -test- ignores it,
* so we have to go all the way and do a -post- from scratch

* Obtain K so we can obtain DoF = N - K - kk
* This is already done by regress EXCEPT when clustering
* (but we still need the unclustered version for r2_a, etc.)
	_rmcoll `indepvars' `avgevars' `weightexp', forcedrop
	local varlist = r(varlist)
	if ("`varlist'"==".") local varlist
	local K : list sizeof varlist

* Run -regress-
	local subcmd regress `vars' `weightexp', `vceoption' `suboptions' `nocons' noheader notable
	Debug, level(3) msg("Subcommand: " in ye "`subcmd'")
	qui `subcmd'
	
	local N = e(N) // We couldn't just use c(N) due to possible frequency weights
	local WrongDoF = `N' - `addconstant' - `K'
	local CorrectDoF = `WrongDoF' - `kk' // kk = Absorbed DoF
	if ("`vcetype'"!="cluster") Assert e(df_r)==`WrongDoF', msg("e(df_r) doesn't match: `e(df_r)'!=`WrongDoF'")


* Store results for the -ereturn post-
	tempname b V
	matrix `b' = e(b)
	matrix `V' = e(V)
	local N = e(N)
	local marginsok = e(marginsok)
	local rmse = e(rmse)
	local rss = e(rss)
	local tss = e(mss) + e(rss) // Regress doesn't report e(tss)
	local N_clust = e(N_clust)

	local predict = e(predict)
	local cmd = e(cmd)
	local cmdline = e(cmdline)
	local title = e(title)

	* Fix V
	if (`K'>0) matrix `V' = `V' * (`WrongDoF' / `CorrectDoF')

	* DoF
	if ("`vcetype'"=="cluster") Assert e(df_r) == e(N_clust) - 1
	local df_r = cond( "`vcetype'"=="cluster" , e(df_r) , max( `CorrectDoF' , 0 ) )

	capture ereturn post `b' `V' `weightexp', dep(`depvar') obs(`N') dof(`df_r') properties(b V)
	local rc = _rc
	Assert inlist(_rc,0,504), msg("error `=_rc' when adjusting the VCV") // 504 = Matrix has MVs
	Assert `rc'==0, msg("Error: estimated variance-covariance matrix has missing values")
	ereturn local marginsok = "`marginsok'"
	ereturn local predict = "`predict'"
	ereturn local cmd = "`cmd'"
	ereturn local cmdline = "`cmdline'"
	ereturn local title = "`title'"
	ereturn local clustvar = "`clustervars'"
	ereturn scalar rmse = `rmse'
	ereturn scalar rss = `rss'
	ereturn scalar tss = `tss'
	if (`N_clust'<.) ereturn scalar N_clust = `N_clust'
	if (`N_clust'<.) ereturn scalar N_clust1 = `N_clust'
	ereturn `hidden' scalar unclustered_df_r = `CorrectDoF' // Used later in R2 adj

* Compute model F-test
	if (`K'>0) {
		qui test `indepvars' `avge' // Wald test
		ereturn scalar F = r(F)
		ereturn scalar df_m = r(df)
		ereturn scalar rank = r(df)+1 // Add constant
		if missing(e(F)) di as error "WARNING! Missing FStat"
	}
	else {
		ereturn scalar F = 0
		ereturn scalar df_m = 0
		ereturn scalar rank = 1
	}

	mata: st_local("original_vars", strtrim(stritrim( "`original_depvar' `original_indepvars' `avge_targets' `original_absvars'" )) )
end

program define Wrapper_mwc, eclass
* This will compute an ols regression with 2+ clusters
syntax , depvar(varname) [indepvars(varlist) avgevars(varlist)] ///
	original_absvars(string) original_depvar(string) [original_indepvars(string) avge_targets(string)] ///
	vceoption(string asis) ///
	kk(integer) ///
	[weightexp(string)] ///
	addconstant(integer) ///
	[SUBOPTions(string)] [*] // [*] are ignored!

	if ("`options'"!="") Debug, level(3) msg("(ignored options: `options')")
	mata: st_local("vars", strtrim(stritrim( "`depvar' `indepvars' `avgevars'" )) ) // Just for esthetic purposes
	if (`c(version)'>=12) local hidden hidden

* Parse contents of VCE()
	local 0 `vceoption'
	syntax namelist(max=11) // Of course clustering by anything beyond 2-3 is insane
	gettoken vcetype clustervars : namelist
	assert "`vcetype'"=="cluster"
	local clustervars `clustervars' // Trim

* Hide constant
	if (!`addconstant') {
		local nocons noconstant
		local kk = `kk' + 1
	}

* Obtain e(b), e(df_m), and resids
	local subcmd regress `depvar' `indepvars' `avgevars' `weightexp', `nocons'
	Debug, level(3) msg("Subcommand: " in ye "`subcmd'")
	qui `subcmd'

	local K = e(df_m)
	local WrongDoF = e(df_r)

	* Store some results for the -ereturn post-
	tempname b
	matrix `b' = e(b)
	local N = e(N)
	local marginsok = e(marginsok)
	local rmse = e(rmse)
	local rss = e(rss)
	local tss = e(mss) + e(rss) // Regress doesn't report e(tss)

	local predict = e(predict)
	local cmd = e(cmd)
	local cmdline = e(cmdline)
	local title = e(title)

	* Compute the bread of the sandwich D := inv(X'X/N)
	tempname XX invSxx
	qui mat accum `XX' = `indepvars' `avgevars' `weightexp', `nocons'
	mat `invSxx' = syminv(`XX') // This line is different from <Wrapper_avar>

	* Resids
	tempvar resid
	predict double `resid', resid

	* DoF
	local df_r = max( `WrongDoF' - `kk' , 0 )

* Use MWC to get meat of sandwich "M" (notation: V = DMD)
	local size = rowsof(`invSxx')
	tempname M V // Will store the Meat and the final Variance
	matrix `V' = J(`size', `size', 0)

* This gives all the required combinations of clustervars (ssc install tuples)
	tuples `clustervars' // create locals i) ntuples, ii) tuple1 .. tuple#
	tempvar group
	local N_clust = .
	local j 0

	forval i = 1/`ntuples' {
		matrix `M' =  `invSxx'
		local vars `tuple`i''
		local numvars : word count `vars'
		local sign = cond(mod(`numvars', 2), "+", "-") // + with odd number of variables, - with even

		GenerateID `vars', gen(`group')
		
		if (`numvars'==1) {
			su `group', mean
			local ++j
			local h : list posof "`vars'" in clustervars
			local N_clust`h' = r(max)

			local N_clust = min(`N_clust', r(max))
			Debug, level(2) msg(" - multi-way-clustering: `vars' has `r(max)' groups")
		}
		
		* Compute the full sandwich (will be saved in `M')

		_robust `resid' `weightexp', variance(`M') minus(0) cluster(`group') // Use minus==1 b/c we adjust the q later
		Debug, level(3) msg(as result "`sign' `vars'")
		* Add it to the other sandwiches
		matrix `V' = `V' `sign' `M'
		drop `group'
	}

	local N_clustervars = `j'

* If the VCV matrix is not positive-semidefinite, use the fix from
* Cameron, Gelbach & Miller - Robust Inference with Multi-way Clustering (JBES 2011)
* 1) Use eigendecomposition V = U Lambda U' where U are the eigenvectors and Lambda = diag(eigenvalues)
* 2) Replace negative eigenvalues into zero and obtain FixedLambda
* 3) Recover FixedV = U * FixedLambda * U'
* This will fail if V is not symmetric (we could use -mata makesymmetric- to deal with numerical precision errors)

	mata: fix_psd("`V'") // This will update `V' making it PSD
	assert inlist(`eigenfix', 0, 1)
	if (`eigenfix') Debug, level(0) msg("Warning: VCV matrix was non-positive semi-definite; adjustment from Cameron, Gelbach & Miller applied.")

	local M = `N_clust' // cond( `N_clust' < . , `N_clust' , `N' )
	local q = ( `N' - 1 ) / `df_r' * `M' / (`M' - 1) // General formula, from Stata PDF
	matrix `V' = `V' * `q'

	* At this point, we have the true V and just need to add it to e()

	local unclustered_df_r = `df_r' // Used later in R2 adj
	local df_r = `M' - 1 // Cluster adjustment

	capture ereturn post `b' `V' `weightexp', dep(`depvar') obs(`N') dof(`df_r') properties(b V)

	local rc = _rc
	Assert inlist(_rc,0,504), msg("error `=_rc' when adjusting the VCV") // 504 = Matrix has MVs
	Assert `rc'==0, msg("Error: estimated variance-covariance matrix has missing values")
	ereturn local marginsok = "`marginsok'"
	ereturn local predict = "`predict'"
	ereturn local cmd = "`cmd'"
	ereturn local cmdline = "`cmdline'"
	ereturn local title = "`title'"
	ereturn scalar rmse = `rmse'
	ereturn scalar rss = `rss'
	ereturn scalar tss = `tss'
	ereturn `hidden' scalar unclustered_df_r = `unclustered_df_r'

	ereturn local clustvar = "`clustervars'"
	assert `N_clust'<.
	ereturn scalar N_clust = `N_clust'
	forval i = 1/`N_clustervars' {
		ereturn scalar N_clust`i' = `N_clust`i''
	}

* Compute model F-test
	if (`K'>0) {
		qui test `indepvars' `avge' // Wald test
		if (r(drop)==1) Debug, level(0) msg("Warning: Some variables were dropped by the F test due to collinearity (or insufficient number of clusters).")
		ereturn scalar F = r(F)
		ereturn scalar df_m = r(df)
		ereturn scalar rank = r(df)+1 // Add constant
		if missing(e(F)) di as error "WARNING! Missing FStat"
	}
	else {
		ereturn scalar F = 0
		ereturn df_m = 0
		ereturn scalar rank = 1
	}

* ereturns specific to this command
	mata: st_local("original_vars", strtrim(stritrim( "`original_depvar' `original_indepvars' `avge_targets' `original_absvars'" )) )

end

program define Wrapper_avar, eclass
	syntax , depvar(varname) [indepvars(varlist) avgevars(varlist)] ///
		original_absvars(string) original_depvar(string) [original_indepvars(string) avge_targets(string)] ///
		vceoption(string asis) ///
		kk(integer) ///
		[weightexp(string)] ///
		addconstant(integer) ///
		[SUBOPTions(string)] [*] // [*] are ignored!

	if ("`options'"!="") Debug, level(3) msg("(ignored options: `options')")
	mata: st_local("vars", strtrim(stritrim( "`depvar' `indepvars' `avgevars'" )) ) // Just for esthetic purposes
	if (`c(version)'>=12) local hidden hidden

	local tmpweightexp = subinstr("`weightexp'", "[pweight=", "[aweight=", 1)

* Convert -vceoption- to what -avar- expects
	local 0 `vceoption'
	syntax namelist(max=3) , [bw(integer 1) dkraay(integer 1) kernel(string) kiefer]
	gettoken vcetype clustervars : namelist
	local clustervars `clustervars' // Trim
	Assert inlist("`vcetype'", "unadjusted", "robust", "cluster")
	local vceoption = cond("`vcetype'"=="unadjusted", "", "`vcetype'")
	if ("`clustervars'"!="") local vceoption `vceoption'(`clustervars')
	if (`bw'>1) local vceoption `vceoption' bw(`bw')
	if (`dkraay'>1) local vceoption `vceoption' dkraay(`dkraay')
	if ("`kernel'"!="") local vceoption `vceoption' kernel(`kernel')
	if ("`kiefer'"!="") local vceoption `vceoption' kiefer

* Hide constant
	if (!`addconstant') {
		local nocons noconstant
		local kk = `kk' + 1
	}

* Before -avar- we need:
*	i) inv(X'X)
*	ii) DoF lost due to included indepvars
*	iii) resids
* Note: It would be shorter to use -mse1- (b/c then invSxx==e(V)*e(N)) but then I don't know e(df_r)
	local subcmd regress `depvar' `indepvars' `avgevars' `weightexp', `nocons'
	Debug, level(3) msg("Subcommand: " in ye "`subcmd'")
	qui `subcmd'
	qui cou if !e(sample)
	assert r(N)==0

	local K = e(df_m) // Should also be equal to e(rank)+1
	local WrongDoF = e(df_r)

	* Store some results for the -ereturn post-
	tempname b
	matrix `b' = e(b)
	local N = e(N)
	local marginsok = e(marginsok)
	local rmse = e(rmse)
	local rss = e(rss)
	local tss = e(mss) + e(rss) // Regress doesn't report e(tss)

	local predict = e(predict)
	local cmd = e(cmd)
	local cmdline = e(cmdline)
	local title = e(title)

	* Compute the bread of the sandwich inv(X'X/N)
	tempname XX invSxx
	qui mat accum `XX' = `indepvars' `avgevars' `tmpweightexp', `nocons'
	* WHY DO I NEED TO REPLACE PWEIGHT WITH AWEIGHT HERE?!?
	
	* (Is this precise enough? i.e. using -matrix- commands instead of mata?)
	mat `invSxx' = syminv(`XX' * 1/`N')
	
	* Resids
	tempvar resid
	predict double `resid', resid

	* DoF
	local df_r = max( `WrongDoF' - `kk' , 0 )

* Use -avar- to get meat of sandwich
	local subcmd avar `resid' (`indepvars' `avgevars') `weightexp', `vceoption' `nocons' // dofminus(0)
	Debug, level(3) msg("Subcommand: " in ye "`subcmd'")
	cap `subcmd'
	local rc = _rc
	if (`rc') {
		di as error "Error in -avar- module:"
		noi `subcmd'
		exit 198
	}

	local N_clust = r(N_clust)
	local N_clust1 = cond(r(N_clust1)<., r(N_clust1), r(N_clust))
	local N_clust2 = r(N_clust2)

* Get the entire sandwich
	* Without clusters it's as if every obs. is is own cluster
	local M = cond( r(N_clust) < . , r(N_clust) , r(N) )
	local q = ( `N' - 1 ) / `df_r' * `M' / (`M' - 1) // General formula, from Stata PDF
	tempname V

	* A little worried about numerical precision
	matrix `V' = `invSxx' * r(S) * `invSxx' / r(N) // Large-sample version
	matrix `V' = `V' * `q' // Small-sample adjustments
	* At this point, we have the true V and just need to add it to e()

* Avoid corner case error when all the RHS vars are collinear with the FEs
	local unclustered_df_r = `df_r' // Used later in R2 adj
	if (`dkraay'>1) local clustervars "`_dta[_TStvar]'"
	if ("`clustervars'"!="") local df_r = `M' - 1

	capture ereturn post `b' `V' `weightexp', dep(`depvar') obs(`N') dof(`df_r') properties(b V)

	local rc = _rc
	Assert inlist(_rc,0,504), msg("error `=_rc' when adjusting the VCV") // 504 = Matrix has MVs
	Assert `rc'==0, msg("Error: estimated variance-covariance matrix has missing values")
	ereturn local marginsok = "`marginsok'"
	ereturn local predict = "`predict'"
	ereturn local cmd = "`cmd'"
	ereturn local cmdline = "`cmdline'"
	ereturn local title = "`title'"
	ereturn local clustvar = "`clustervars'"

	ereturn scalar rmse = `rmse'
	ereturn scalar rss = `rss'
	ereturn scalar tss = `tss'
	if ("`N_clust'"!="") ereturn scalar N_clust = `N_clust'
	if ("`N_clust1'"!="") ereturn scalar N_clust1 = `N_clust1'
	if ("`N_clust2'"!="") ereturn scalar N_clust2 = `N_clust2'
	ereturn `hidden' scalar unclustered_df_r = `unclustered_df_r'

	if (`bw'>1) {
		ereturn scalar bw = `bw'
		if ("`kernel'"=="") local kernel Bartlett // Default
	}
	if ("`kernel'"!="") ereturn local kernel = "`kernel'"
	if ("`kiefer'"!="") ereturn local kiefer = "`kiefer'"
	if (`dkraay'>1) ereturn scalar dkraay = `dkraay'

* Compute model F-test
	if (`K'>0) {
		qui test `indepvars' `avge' // Wald test
		if (r(drop)==1) Debug, level(0) msg("Warning: Some variables were dropped by the F test due to collinearity (or insufficient number of clusters).")
		ereturn scalar F = r(F)
		ereturn scalar df_m = r(df)
		ereturn scalar rank = r(df)+1 // Add constant
		if missing(e(F)) di as error "WARNING! Missing FStat"
	}
	else {
		ereturn scalar F = 0
		ereturn df_m = 0
		ereturn scalar rank = 1
	}
	
* ereturns specific to this command
	mata: st_local("original_vars", strtrim(stritrim( "`original_depvar' `original_indepvars' `avge_targets' `original_absvars'" )) )
end

program define Wrapper_ivregress, eclass
	syntax , depvar(varname) endogvars(varlist) instruments(varlist) ///
		[indepvars(varlist) avgevars(varlist)] ///
		original_depvar(string) original_endogvars(string) original_instruments(string) ///
		[original_indepvars(string) avge_targets(string)] ///
		vceoption(string asis) ///
		KK(integer) ///
		[weightexp(string)] ///
		addconstant(integer) ///
		SHOWRAW(integer) first(integer) vceunadjusted(integer) ///
		[ESTimator(string) TWICErobust(string)] ///
		[SUBOPTions(string)] [*] // [*] are ignored!

	if ("`options'"!="") Debug, level(3) msg("(ignored options: `options')")
	mata: st_local("vars", strtrim(stritrim( "`depvar' `indepvars' `avgevars' (`endogvars'=`instruments')" )) )
	if (`c(version)'>=12) local hidden hidden

	local opt_estimator = cond("`estimator'"=="gmm2s", "gmm", "`estimator'")

	* Convert -vceoption- to what -ivreg2- expects
	local 0 `vceoption'
	syntax namelist(max=2)
	gettoken vceoption clustervars : namelist
	local clustervars `clustervars' // Trim
	Assert inlist("`vceoption'", "unadjusted", "robust", "cluster")
	if ("`clustervars'"!="") local vceoption `vceoption' `clustervars'
	local vceoption "vce(`vceoption')"

	if ("`estimator'"=="gmm2s") {
		local wmatrix : subinstr local vceoption "vce(" "wmatrix("
		if ("`twicerobust'"=="") {
			local vceoption = cond(`vceunadjusted', "vce(unadjusted)", "")			
		}
	}
	
	* Note: the call to -ivregress- could be optimized.
	* EG: -ivregress- calls ereturn post .. ESAMPLE(..) but we overwrite the esample and its SLOW
	* But it's a 1700 line program so let's not worry about it

* Hide constant
	if (!`addconstant') {
		local nocons noconstant
		local kk = `kk' + 1
	}

* Show first stage
	if (`first') {
		local firstoption "first"
	}

* Subcmd
	local subcmd ivregress `opt_estimator' `vars' `weightexp', `wmatrix' `vceoption' small `nocons' `firstoption' `suboptions'
	Debug, level(3) msg("Subcommand: " in ye "`subcmd'")
	local noise = cond(`showraw', "noi", "qui")
	`noise' `subcmd'
	
	* Fix DoF if needed
	local N = e(N)
	local K = e(df_m)
	local WrongDoF = `N' - `addconstant' - `K'
	local CorrectDoF = `WrongDoF' - `kk'
	Assert !missing(`CorrectDoF')

	* We should have used M/M-1 instead of N/N-1, but we are making ivregress to do the wrong thing by using vce(unadjusted) (which makes it fit with ivreg2)
	local q 1
	if ("`estimator'"=="gmm2s" & "`clustervars'"!="") {
		local N = e(N)
		tempvar group
		GenerateID `clustervars', gen(`group')
		su `group', mean
		drop `group'
		local M = r(max) // N_clust
		local q = ( `M' / (`M' - 1)) / ( `N' / (`N' - 1) ) // multiply correct, divide prev wrong one
		ereturn scalar df_r = `M' - 1
	}

	tempname V
	matrix `V' = e(V) * (`WrongDoF' / `CorrectDoF') * `q'
	ereturn repost V=`V'
	
	if ("`clustervars'"=="") ereturn scalar df_r = `CorrectDoF'

	* ereturns specific to this command
	mata: st_local("original_vars", strtrim(stritrim( "`original_depvar' `original_indepvars' `avge_targets' `original_absvars' (`original_endogvars'=`original_instruments')" )) )
	ereturn scalar F = e(F) * `CorrectDoF' / `WrongDoF'

	ereturn scalar tss = e(mss) + e(rss) // ivreg2 doesn't report e(tss)
	ereturn `hidden' scalar unclustered_df_r = `CorrectDoF' // Used later in R2 adj
end

program define Wrapper_ivreg2, eclass
	syntax , depvar(varname) endogvars(varlist) instruments(varlist) ///
		[indepvars(varlist) avgevars(varlist)] ///
		original_depvar(string) original_endogvars(string) original_instruments(string) ///
		[original_indepvars(string) avge_targets(string)] ///
		[original_absvars(string) avge_targets] ///
		vceoption(string asis) ///
		KK(integer) ///
		[SHOWRAW(integer 0)] first(integer) [weightexp(string)] ///
		addconstant(integer) ///
		[ESTimator(string)] ///
		[SUBOPTions(string)] [*] // [*] are ignored!
	if ("`options'"!="") Debug, level(3) msg("(ignored options: `options')")
	if (`c(version)'>=12) local hidden hidden
	
	* Disable some options
	local 0 , `suboptions'
	syntax , [SAVEFPrefix(name)] [*] // Will ignore SAVEFPREFIX
	local suboptions `options'
	assert (`addconstant'==0)

	* Convert -vceoption- to what -ivreg2- expects
	local 0 `vceoption'
	syntax namelist(max=3) , [bw(string) dkraay(string) kernel(string) kiefer]
	gettoken vcetype clustervars : namelist
	local clustervars `clustervars' // Trim
	Assert inlist("`vcetype'", "unadjusted", "robust", "cluster")
	local vceoption = cond("`vcetype'"=="unadjusted", "", "`vcetype'")
	if ("`clustervars'"!="") local vceoption `vceoption'(`clustervars')
	if ("`bw'"!="") local vceoption `vceoption' bw(`bw')
	if ("`dkraay'"!="") local vceoption `vceoption' dkraay(`dkraay')
	if ("`kernel'"!="") local vceoption `vceoption' kernel(`kernel')
	if ("`kiefer'"!="") local vceoption `vceoption' kiefer
	
	mata: st_local("vars", strtrim(stritrim( "`depvar' `indepvars' `avgevars' (`endogvars'=`instruments')" )) )
	
	if (`first') {
		local firstoption "first savefirst"
	}

	if ("`estimator'"!="2sls") local opt_estimator `estimator'
	
	* Variables have already been demeaned, so we need to add -nocons- or the matrix of orthog conditions will be singular
	if ("`cue'"=="") {
		local nocons nocons // Exception to get the same results as ivreg2, partial
	}
	else {
		local nocons nocons // partial(cons)
	}

	local subcmd ivreg2 `vars' `weightexp', `vceoption' `firstoption' small sdofminus(`=`kk'+1') `nocons' `opt_estimator' `suboptions'
	Debug, level(3) msg(_n "call to subcommand: " _n as result "`subcmd'")
	local noise = cond(`showraw', "noi", "qui")
	`noise' `subcmd'
	if ("`noise'"=="noi") di in red "{hline 64}" _n "{hline 64}"
	ereturn scalar tss = e(mss) + e(rss) // ivreg2 doesn't report e(tss)
	ereturn scalar unclustered_df_r = e(N) - e(df_m)

	if ("`e(vce)'"=="robust cluster") ereturn local vce = "cluster"

	if !missing(e(ecollin)) {
		di as error "endogenous covariate <`e(ecollin)'> was perfectly predicted by the instruments!"
		error 2000
	}

	if (`first') {
		ereturn `hidden' local first_prefix = "_ivreg2_"
		ereturn `hidden' local ivreg2_firsteqs = e(firsteqs)
		ereturn local firsteqs
	}

	foreach cat in exexog insts instd {
		FixVarnames `e(`cat')'
		ereturn local `cat' = "`r(newnames)'"
	}

	if (`first') {
		* May be a problem if we ran out of space for storing estimates
		local ivreg2_firsteqs "`e(ivreg2_firsteqs)'"
		tempname hold
		estimates store `hold' , nocopy
		foreach fs_eqn in `ivreg2_firsteqs' {
			qui estimates restore `fs_eqn'
			FixVarnames `e(depvar)'
			ereturn local depvar = r(prettynames)
			FixVarnames `e(inexog)'
			ereturn local inexog = r(prettynames)

			tempname b
			matrix `b' = e(b)
			local backup_colnames : colnames `b'
			FixVarnames `backup_colnames'
			matrix colnames `b' = `r(prettynames)' // newnames? prettynames?
			ereturn repost b=`b', rename

			estimates store `fs_eqn', nocopy
		}
		qui estimates restore `hold'
		estimates drop `hold'

	}

	* ereturns specific to this command
	mata: st_local("original_vars", strtrim(stritrim( "`original_depvar' `original_indepvars' `avge_targets' `original_absvars' (`original_endogvars'=`original_instruments')" )) )
end

program define AddConstant
	syntax varlist(numeric)
	foreach var of local varlist {
		local mean : char `var'[mean]
		assert "`mean'"!=""
		assert !missing(`mean')
		qui replace `var' = `var' + `mean'
	}
end

program define Attach, eclass
	syntax, [NOTES(string)] [statsmatrix(string)] summarize_quietly(integer)
	
	* Summarize
	* This needs to happen after all the missing obs. have been dropped and the only obs. are those that *WILL* be in the regression
	if ("`statsmatrix'"!="") {
		* Update beta vector
		* ...

		ereturn matrix summarize = `statsmatrix', copy // If we move instead of copy, stages() will fail
		if (!`summarize_quietly' & "`statsmatrix'"!="") {
			di as text _n "{sf:Regression Summary Statistics}" _c
			matlist e(summarize)', border(top bottom) twidth(18) rowtitle(Variable)
		}
	}

	* Parse key=value options and append to ereturn as hidden
	mata: st_local("notes", strtrim(`"`notes'"')) // trim (supports large strings)
	local keys
	while (`"`notes'"'!="") {
		gettoken key notes : notes, parse(" =")
		Assert !inlist("`key'","sample","time"), msg("Key cannot be -sample- or -time-") // Else -estimates- will fail
		gettoken _ notes : notes, parse("=")
		gettoken value notes : notes, quotes
		local keys `keys' `key'
		ereturn hidden local `key' `value'
	}
	if ("`keys'"!="") ereturn hidden local keys `keys'

end


// -------------------------------------------------------------
// Display Regression Table
// -------------------------------------------------------------

 program define Replay, eclass
	syntax , [*]
	Assert e(cmd)=="reghdfe"
	local subcmd = e(subcmd)
	Assert "`subcmd'"!="" , msg("e(subcmd) is empty")

	* Add pretty names for AvgE variables
	tempname b
	matrix `b' = e(b)
	local backup_colnames : colnames `b'
	matrix colnames `b' = `e(prettynames)'
	local savefirst = e(savefirst)
	local suboptions = e(suboptions)

	local diopts = "`e(diopts)'"
	if ("`options'"!="") { // Override
		_get_diopts diopts /* options */, `options'
	}

	if ("`subcmd'"=="ivregress") {
		* Don't want to display anova table or footnote
		_coef_table_header
		_coef_table, `diopts' bmatrix(`b') vmatrix(e(V)) // plus 
	}
	else if ("`subcmd'"=="ivreg2") {
		* Backup before showing both first and second stage
		tempname hold
		
		if ("`e(ivreg2_firsteqs)'"!="") {
			estimates store `hold'

			local i 0
			foreach fs_eqn in `e(ivreg2_firsteqs)' {
				local instrument  : word `++i' of `e(instd)'
				di as input _n "{title:First stage for `instrument'}"
				estimates replay `fs_eqn' , nohead `diopts'
				if (!`savefirst') estimates drop `fs_eqn'
			}

			ereturn clear
			qui estimates restore `hold'
			di as input _n "{title:Second stage}"
		}

		estimates store `hold'
		ereturn repost b=`b', rename
		ereturn local cmd = "`subcmd'"
		`subcmd' , `diopts'
		ereturn clear // Need this because -estimates restore- behaves oddly
		qui estimates restore `hold'
		assert e(cmd)=="reghdfe"
		estimates drop `hold'


		*ereturn local cmd = "reghdfe"
		*matrix `b' = e(b)
		*matrix colnames `b' = `backup_colnames'
		*ereturn repost b=`b', rename
	}
	else {

		* Regress-specific code, because it doesn't play nice with ereturn
		sreturn clear 

		if "`e(prefix)'" != "" {
			_prefix_display, `diopts'
			exit
		}
		

		*_coef_table_header
		Header

		di
		local plus = cond(e(model)=="ols" & inlist("`e(vce)'", "unadjusted", "ols"), "plus", "")
		_coef_table, `plus' `diopts' bmatrix(`b') vmatrix(e(V))
	}

	reghdfe_footnote
	* Revert AvgE else -predict- and other commands will choke


end


* (Modified from _coef_table_header.ado)

program define Header
	if !c(noisily) exit

	tempname left right
	.`left' = {}
	.`right' = {}

	local width 78
	local colwidths 1 30 51 67
	local i 0
	foreach c of local colwidths {
		local ++i
		local c`i' `c'
		local C`i' _col(`c')
	}

	local c2wfmt 10
	local c4wfmt 10
	local max_len_title = `c3' - 2
	local c4wfmt1 = `c4wfmt' + 1
	local title  `"`e(title)'"'
	local title2 `"`e(title2)'"'
	local title3 `"`e(title3)'"'
	local title4 `"`e(title4)'"'
	local title5 `"`e(title5)'"'

	// Right hand header ************************************************

	*N obs
	.`right'.Arrpush `C3' "Number of obs" `C4' "= " as res %`c4wfmt'.0f e(N)

	* Ftest
	if `"`e(chi2)'"' != "" | "`e(df_r)'" == "" {
		Chi2test `right' `C3' `C4' `c4wfmt'
	}
	else {
		Ftest `right' `C3' `C4' `c4wfmt'
	}

	* display R-squared
	if !missing(e(r2)) {
		.`right'.Arrpush `C3' "R-squared" `C4' "= " as res %`c4wfmt'.4f e(r2)
	}
	*if !missing(e(r2_p)) {
	*	.`right'.Arrpush `C3' "Pseudo R2" `C4' "= " as res %`c4wfmt'.4f e(r2_p)
	*}
	if !missing(e(r2_a)) {
		.`right'.Arrpush `C3' "Adj R-squared" `C4' "= " as res %`c4wfmt'.4f e(r2_a)
	}
	if !missing(e(r2_within)) {
		.`right'.Arrpush `C3' "Within R-sq." `C4' "= " as res %`c4wfmt'.4f e(r2_within)
	}
	if !missing(e(rmse)) {
		.`right'.Arrpush `C3' "Root MSE" `C4' "= " as res %`c4wfmt'.4f e(rmse)
	}

	// Left hand header *************************************************

	* make title line part of the header if it fits
	local len_title : length local title
	forv i=2/5 {
		if (`"`title`i''"'!="") {
			local len_title = max(`len_title',`:length local title`i'')
		}
	}
	
	if `len_title' < `max_len_title' {
		.`left'.Arrpush `"`"`title'"'"'
		local title
		forv i=2/5 {
			if `"`title`i''"' != "" {
					.`left'.Arrpush `"`"`title`i''"'"'
					local title`i'
			}
		}
		.`left'.Arrpush "" // Empty
	}

	* Clusters
	local kr = `.`right'.arrnels' // number of elements in the right header
	local kl = `.`left'.arrnels' // number of elements in the left header
	local N_clustervars = e(N_clustervars)
	if (`N_clustervars'==.) local N_clustervars 0
	local space = `kr' - `kl' - `N_clustervars'
	local clustvar = e(clustvar)
	forv i=1/`space' {
		.`left'.Arrpush ""
	}
	forval i = 1/`N_clustervars' {
		gettoken cluster clustvar : clustvar
		local num = e(N_clust`i')
		.`left'.Arrpush `C1' "Number of clusters (" as res "`cluster'" as text  ") " `C2' as text "= " as res %`c2wfmt'.0f `num'
	}
	
	HeaderDisplay `left' `right' `"`title'"' `"`title2'"' `"`title3'"' `"`title4'"' `"`title5'"'
end

program define HeaderDisplay
		args left right title1 title2 title3 title4 title5

		local nl = `.`left'.arrnels'
		local nr = `.`right'.arrnels'
		local K = max(`nl',`nr')

		di
		if `"`title1'"' != "" {
				di as txt `"`title'"'
				forval i = 2/5 {
					if `"`title`i''"' != "" {
							di as txt `"`title`i''"'
					}
				}
				if `K' {
						di
				}
		}

		local c _c
		forval i = 1/`K' {
				di as txt `.`left'[`i']' as txt `.`right'[`i']'
		}
end

program define Ftest
		args right C3 C4 c4wfmt is_svy

		local df = e(df_r)
		if !missing(e(F)) {
				.`right'.Arrpush                                ///
						 `C3' "F("                              ///
				   as res %4.0f e(df_m)                         ///
				   as txt ","                                   ///
				   as res %7.0f `df' as txt ")" `C4' "= "       ///
				   as res %`c4wfmt'.2f e(F)
				.`right'.Arrpush                                ///
						 `C3' "Prob > F" `C4' "= "              ///
				   as res %`c4wfmt'.4f Ftail(e(df_m),`df',e(F))
		}
		else {
				local dfm_l : di %4.0f e(df_m)
				local dfm_l2: di %7.0f `df'
				local j_robust "{help j_robustsingular##|_new:F(`dfm_l',`dfm_l2')}"
				.`right'.Arrpush                                ///
						  `C3' "`j_robust'"                     ///
				   as txt `C4' "= " as result %`c4wfmt's "."
				.`right'.Arrpush                                ///
						  `C3' "Prob > F" `C4' "= " as res %`c4wfmt's "."
		}
end

program define Chi2test

		args right C3 C4 c4wfmt

		local type `e(chi2type)'
		if `"`type'"' == "" {
				local type Wald
		}
		if !missing(e(chi2)) {
				.`right'.Arrpush                                ///
						  `C3' "`type' chi2("                   ///
				   as res e(df_m)                               ///
				   as txt ")" `C4' "= "                         ///
				   as res %`c4wfmt'.2f e(chi2)
				.`right'.Arrpush                                ///
						  `C3' "Prob > chi2" `C4' "= "          ///
				   as res %`c4wfmt'.4f chi2tail(e(df_m),e(chi2))
		}
		else {
				local j_robust                                  ///
				"{help j_robustsingular##|_new:`type' chi2(`e(df_m)')}"
				.`right'.Arrpush                                ///
						  `C3' "`j_robust'"                     ///
				   as txt `C4' "= " as res %`c4wfmt's "."
				.`right'.Arrpush                                ///
						  `C3' "Prob > chi2" `C4' "= "          ///
				   as res %`c4wfmt's "."
		}
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
		local M`g'_nested = 0
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
				local M`g'_nested = 1
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
		return scalar M`g'_nested = `M`g'_nested'
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

