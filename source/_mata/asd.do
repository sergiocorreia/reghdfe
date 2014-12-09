clear all

local Varlist 		string scalar
local Integer 		real scalar
local VarByFE 		real colvector // Should be levels*1
local Series		real colvector // Should be N*1
local Matrix		real matrix
local SharedData 	external struct FixedEffect vector

mata:
mata set matastrict on

// -------------------------------------------------------------
// REGRESS_BY_GROUP: Multivariate regression on constant and at least 1 var.
// -------------------------------------------------------------
// Returns block-column matrix with estimates by group; last estimate is the constant
`Matrix' function regress_by_group(`Series' y, `Matrix' x, `Series' index, 
	`VarByFE' offset, `VarByFE' count, `Matrix' invxx, `Series' sorted_weight)
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
		predicted[| j_lower , 1 \ j_upper , . |] = b[K] :+ tmp_x * b[|1 \ K-1|]
		j_lower = j_upper + 1
	}
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


void function testit() {
	real colvector y, group, weight, index, sum_count, count, sorted_weight
	real matrix v, invxx, predicted
	external real matrix betas
	
	y = st_data(., "y")
	group = st_data(., "g")
	v = st_data(., "v1 v2")
	weight = st_data(., "w")
	index = order(group,1) // index[1] has a member of the first group, and so on
	// group, index
	
	count = count_by_group(group, index)
	sum_count = quadrunningsum(count)
	sorted_weight = weight[index, 1]
	count = count_by_group(group, index, 0, sorted_weight)	
	// count, sum_count
	
	invxx = compute_invxx(v ,index, sum_count, count, sorted_weight)
	predicted = regress_by_group(y,v,index, sum_count, count, invxx, sorted_weight)
	
	//predicted
	betas
	
	
}

end

include reghdfe.mata
*FEs[g].indexfrom0 = order(FEs[g].group,1)


cls
set more off
set obs 25
set seed 4352458
gen y = uniform()
gen v1 = uniform()
gen v2 = uniform()
gen g = uniform()<0.25	
replace g= 1 in 1
replace g = g + g[_n-1] in 2/l
gen w = 1 + int(10 * uniform())
gen rand = uniform()
sort rand
drop rand
mata: testit()

tab g [fw=w] // count_by_group

// calculate_invxx: element 1 , 1+g, and their cross
//reg y ibn.g ibn.g#c.v [fw=w], hascons
//matrix A = e(V) * e(df_r) / e(rss)
//matrix list A

reg y ibn.g ibn.g#c.v1 ibn.g#c.v2 [fw=w], hascons

*br


exit
//	outdata = bivariate_by_group(indata, FEs[g_to].indexfrom0,
//  FEs[g_to].sum_count, FEs[g_to].count,
//		FEs[g_to].v, FEs[g_to].v_meanvar, FEs[g_to].sorted_weight)
	
