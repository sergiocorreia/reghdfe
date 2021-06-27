// --------------------------------------------------------------------------
// Helper functions
// --------------------------------------------------------------------------

mata:

`Matrix' inv_vec(`Vector' x, `Integer' num_cols)
{
	// Invert the vec() function
	//return(num_cols == 1 ? x : transposeonly(rowshape(x, num_cols)))
	return(num_cols == 1 ? x : colshape(x, num_cols))
}


`Vector' fast_vec(`Matrix' x)
{
	//return(cols(x) == 1 ? x :  vec(x))
	return(cols(x) == 1 ? x :  colshape(x, 1))
}


`Matrix' givens_rotation(`Real' a, `Real' b, `Real' r)
{
	// See http://www.netlib.org/lapack/lawnspdf/lawn150.pdf for technical details and an improved alternative
	`Real' c, s
	if (b) {
		r = hypot(a, b) // sqrt(a ^ 2 + b ^ 2)
		c = a / r
		s = b / r
	}
	else {
		// If a and b are zero, then r=0 and c=s=.
		// To avoid these issues, we follow the "Stable Calculation" of https://en.wikipedia.org/wiki/Givens_rotation
		r = a
		c = 1
		s = 0
	}
	return(c, -s \ s, c)
}


`Real' hypot(`Real' x, `Real' y)
{
	// Naive algorithm: sqrt(x ^ 2 + y ^ 2)
	
	// This algorithm: corrected unfused from page 11 of https://arxiv.org/pdf/1904.09481.pdf
	
	// See also:
	// - https://en.wikipedia.org/wiki/Hypot
	// - https://arxiv.org/abs/1904.09481
	// - https://walkingrandomly.com/?p=6633

	`Real' ax, ay, scale, h, delta
	ax = abs(x)
	ay = abs(y)

	// Ensure ax >= ay
	if (ax < ay) swap(ax, ay)

	// Operands vary widely (ay is much smaller than ay)
	if (ay <= ax * sqrt(epsilon(0.5))) return(ax)

	// Operands do not vary widely
	scale = epsilon(sqrt(smallestdouble())) // rescaling constant
	if (ax > sqrt(maxdouble()/2)) {
		ax = ax * scale
		ay = ay * scale
		scale = 1 / scale
	}
	else if (ay < sqrt(smallestdouble())) {
		ax = ax / scale
		ay = ay / scale
	}
	else {
		scale = 1
	}
	h = sqrt(ax^2+ay^2)

	// Single branch
	delta = h - ax
	h = h - (delta*(2*(ax-ay)) + (2*delta-ay)*ay + delta^2) / (2*h)
	return(h*scale)
}


`Void' assign(`Vector' x, `Real' x1, `Real' x2)
{
	// "assign(x, a, b)" is equivalent to "a, b = x"
	assert(rows(x) == 2)
	assert_msg(!hasmissing(x), "input x has missing values") // BUGBUG remove once we are done debugging
	x1 = x[1]
	x2 = x[2]
}


`Vector' normalize(`Vector' x, `Vector' weights, `Real' norm, `String' msg)
{
	`Vector' normalized_x
	assert(weights!=.) // If we don't want to pass any weights just set them to "1", not to missing
	norm = vector_norm(x, weights, msg)

	if (norm < epsilon(1)) {
		norm = 0
		return(x)
	}
	else {
		return(x / norm)
	}
}


`Real' vector_norm(`vector' x, `Vector' weights, `String' msg)
{
	// Compute 2-norm (much faster than calling -norm-)

	`Real'		ans
	if (weights != 1) assert_msg(rows(x) == rows(weights), "weights non-conforming size")
	ans = sqrt(quadcross(x, weights, x))
	assert_msg(!missing(ans), msg) // TODO: remove in production
	return(ans)
	
	// Note: matrix_norm() computed as not used; formula for Frobenius norm was:
	// sqrt(sum(x :^ 2)) == sqrt(trace(cross(x, x)))
}


`Variables' panelmean(`Variables' y, `FE_Factor' f)
{
	assert(0) // replaced by f.panelmean() which feels more natural
	assert_boolean(f.has_weights)
	assert_msg(f.panel_is_setup, "F.panel_setup() must be run beforehand")
	if (f.has_weights) {
		return(editmissing(f.panelsum(y) :/ f.weighted_counts, 0))
	}
	else {
		return(f.panelsum(y) :/ f.counts)
	}
}


`Matrix' precompute_inv_xx(`FE_Factor' f)
{
	`Integer'               i, L, K, offset
	`Vector'             	tmp_x, tmp_w
	`Matrix'                inv_xx, tmp_inv_xx, sorted_x
	`RowVector'             tmp_xmeans

	// note: f.x and f.weights must be already sorted by the factor f
	assert_boolean(f.has_weights)
	assert_boolean(f.has_intercept)

	L = f.num_levels
	K = cols(f.unsorted_x)
	inv_xx = J(L * K, K, .)
	sorted_x = f.sort(f.unsorted_x)

	for (i=1; i<=L; i++) {
		tmp_x = panelsubmatrix(sorted_x, i, f.info)
		tmp_w = f.has_weights ? panelsubmatrix(f.weights, i, f.info) : 1
		if (f.has_intercept) {
			tmp_xmeans = K > 1 ? f.x_means[i, .] : f.x_means[i]
			tmp_inv_xx = invsym(quadcrossdev(tmp_x, tmp_xmeans, tmp_w, tmp_x, tmp_xmeans))
		}
		else {
			tmp_inv_xx = invsym(quadcross(tmp_x, tmp_w, tmp_x))
		}
		offset = K * (i - 1)
		inv_xx[|offset + 1, 1 \ offset + K , . |] = tmp_inv_xx
	}
	return(inv_xx)
}


`Variables' panelsolve_invsym(`Variables' y, `FE_Factor' f, | `Matrix' alphas)
{
	`Integer'               i, L, K, offset
	`Boolean'               save_alphas
	`Matrix'             	xbd, tmp_xbd, tmp_x, tmp_y, tmp_xy, tmp_inv_xx, sorted_x
	`Vector'             	tmp_w, b
	`RowVector'             tmp_xmeans, tmp_ymeans

	// note: y, f.x, and f.weights must be already sorted by the factor f
	assert_boolean(f.has_weights)
	assert_boolean(f.has_intercept)
	assert_msg(rows(f.unsorted_x) == rows(y))

	save_alphas = args()>=3 & alphas!=J(0,0,.)
	if (save_alphas) assert(cols(y)==1)

	L = f.num_levels
	K = cols(f.unsorted_x)
	xbd = J(rows(y), cols(y), .)
	sorted_x = f.sort(f.unsorted_x)

	for (i=1; i<=L; i++) {
		tmp_y = panelsubmatrix(y, i, f.info)
		tmp_x = panelsubmatrix(sorted_x, i, f.info)
		tmp_w = f.has_weights ? panelsubmatrix(f.weights, i, f.info) : 1
		offset = K * (i - 1)
		tmp_inv_xx = f.inv_xx[|offset + 1, 1 \ offset + K , . |]

		if (f.has_intercept) {
			tmp_ymeans = mean(tmp_y, tmp_w)
			tmp_xmeans = K > 1 ? f.x_means[i, .] : f.x_means[i]
			tmp_xy = quadcrossdev(tmp_x, tmp_xmeans, tmp_w, tmp_y, tmp_ymeans)
			b = tmp_inv_xx * tmp_xy
			tmp_xbd = (tmp_x :- tmp_xmeans) * b :+ tmp_ymeans
			if (save_alphas) alphas[i, .] = tmp_ymeans - tmp_xmeans * b, b'
		}
		else {
			tmp_xy = quadcross(tmp_x, tmp_w, tmp_y)
			b = tmp_inv_xx * tmp_xy
			tmp_xbd = tmp_x * b
			if (save_alphas) alphas[i, .] = b'
		}
		xbd[|f.info[i,1], 1 \ f.info[i,2], .|] = tmp_xbd
	}
	return(f.invsort(xbd))
}


`RowVector' function compute_stdevs(`Matrix' A)
{
	`RowVector'				stdevs
	`Integer'				K, N
	// Note: each col of A will have stdev of 1 unless stdev is quite close to 0

	// We don't need good accuracy for the stdevs, so we have a few alternatives:
	// Note: cross(1,A) is the same as colsum(A), but faster
	// Note: cross(A, A) is very fast, but we only need the main diagonals
	// [A: 1sec] stdevs = sqrt( (colsum(A:*A) - (cross(1, A) :^ 2 / N)) / (N-1) )
	// [B: .61s] stdevs = sqrt( (diagonal(cross(A, A))' - (cross(1, A) :^ 2 / N)) / (N-1) )
	// [C: .80s] stdevs = diagonal(sqrt(variance(A)))'
	// [D: .67s] means = cross(1, A) / N; stdevs =  sqrt(diagonal(crossdev(A, means, A, means))' / (N-1))

	N = rows(A)
	K = cols(A)
	stdevs = J(1, K, .)

	// (A) Very precise

	// (B) Precise
	// means = cross(1, A) / N
	// stdevs =  sqrt(diagonal(quadcrossdev(A, means, A, means))' / (N-1))

	// (C) 20% faster; don't use it if you care about accuracy
	stdevs = sqrt( (diagonal(cross(A, A))' - (cross(1, A) :^ 2 / N)) / (N-1) )
	assert_msg(!missing(stdevs), "stdevs are missing; is N==1?") // Shouldn't happen as we don't expect N==1
	stdevs = colmax(( stdevs \ J(1, K, 1e-3) ))

	// (D) Equilibrate matrix columns instead of standardize (i.e. just divide by column max)
	// _perhapsequilc(A, stdevs=.)
	// stdevs = 1 :/ stdevs
	// assert_msg(!missing(stdevs), "stdevs are missing; is N==1?")

	// (E) Don't do anything
	// stdevs = J(1, cols(A), 1)

	return(stdevs)
}


`RowVector' function standardize_data(`Matrix' A)
{
	`RowVector'				stdevs
	assert_msg(!isfleeting(A), "input cannot be fleeting")
	stdevs = compute_stdevs(A)
	A = A :/ stdevs
	return(stdevs)
}


`Variables' st_data_pool(`Vector'  sample,		// observations to load
					     `Varlist' vars,		// variables to load
					     `Integer' poolsize,	// how many variables to load each time
					   | `Boolean' compact,		// whether to trim the dataset to save memory or not
					     `Varlist' keepvars,	// what extra vars to keep if we trim the dataset (clustervars, timevar, panelvar)
					     `Boolean' verbose)
{
	`Integer'               i, j, k
	`Varlist'				temp_keepvars
	`Variables'				data

	if (args()<4 | compact == .) compact = 0
	if (args()<6 | verbose == .) verbose = 0

	k = cols(vars)
	assert_msg(poolsize > 0, "poolsize must be a positive integer")

	if (k <= poolsize) {
		data = st_data(sample, vars)
		if (compact) compact_dataset(keepvars)
	}
	else {
		data = J(rows(sample), 0, .)
		for (i=1; i<=k; i=i+poolsize) {
			j = i + poolsize - 1
			if (j>k) j = k
			data = data, st_data(sample, vars[i..j])
			if (compact) {
				temp_keepvars = j == k ? keepvars :  keepvars , vars[j+1..k]
				compact_dataset(temp_keepvars)
			}
		}
	}

	assert_msg(k == cols(data), "could not load data into Mata correctly; please contact author")
	return(data)
}


`Void' compact_dataset(`Varlist' keepvars)
{
	keepvars = tokens(invtokens(keepvars)) // deal with empty strings
	if (cols(keepvars)) {
		stata(sprintf("fvrevar %s, list", invtokens(keepvars)))
		stata(sprintf("novarabbrev keep %s", st_global("r(varlist)")))
	}
	else {
		stata("clear")
	}
}


`Void' add_undocumented_options(`String' object, `String' options, `Integer' verbose)
{
	`StringRowVector'		tokens
	`Integer'				num_tokens, i
	`String'				key, val, cmd, msg

	if (options == "") return
	tokens = tokens(options, " ()")
	
	msg = sprintf("options {bf:%s} not allowed", options)
	assert_msg(!mod(cols(tokens), 4), msg, 198, 0) // 198: option not allowed; 0:traceback off

	num_tokens = trunc(cols(tokens) / 4)
	for (i=0; i<num_tokens; i++) {
		key = tokens[4*i+1]
		assert_msg(tokens[4*i+2] == "(")
		val = tokens[4*i+3]
		assert_msg(tokens[4*i+4] == ")")
		if (verbose > 1) printf("{txt}    %s.%s = %s\n", object, key, val)
		cmd = sprintf("cap mata: %s.%s = %s", object, key, val)
		stata(cmd)
		assert_msg(!st_numscalar("c(rc)"), sprintf("option {bf:%s} not allowed", key), 198, 0)
	}
}


// This short-but-difficult function receives a list of dropped rows (idx) at the group level
// And returns the list of dropped rows (indiv_idx) at the individual level
`Vector' get_indiv_idx(`Vector' sample, `Vector' indiv_sample, `Varname' group_id, `Vector' idx)
{
	`Factor'				F
	`Vector'				mask, dropped_levels, indiv_idx

	// varnames, touse, verbose, method, SORT_LEVELS, count_levels, hash_ratio, save_keys
	F = factor(group_id, indiv_sample, ., "", 1, 1, ., 0) // not completely sure that we need sort_levels==1
	
	// Map observations to group identifiers
	mask = J(st_nobs(), 1, 0)
	mask[indiv_sample] = F.levels

	// List what groups we are dropping
	dropped_levels = mask[sample[idx]]

	// Create mask indicating whether we drop each group or not
	mask = create_mask(F.num_levels, 0, dropped_levels, 1)

	// Select rows within indiv_sample that correspond to the groups we are dropping
	indiv_idx = selectindex(mask[F.levels])

	return(indiv_idx)
}


// --------------------------------------------------------------------------
// Divide two row vectors but adjust the denominator if it's too small
// --------------------------------------------------------------------------
`RowVector' safe_divide(`RowVector' numerator, `RowVector' denominator, | `Real' epsi) {
	 // If the numerator goes below machine precision, we lose accuracy
	 // If the denominator goes below machine precision, the division explodes
	 if (args()<3 | epsi==.) epsi = epsilon(1)
	return( numerator :/ colmax(denominator \ J(1,cols(denominator),epsi)) )
}

end
