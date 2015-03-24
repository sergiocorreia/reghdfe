* Do not run this test with the normal ones
* It's purpose is just to test/debug the components
* of reghdfe.mata

// -------------------------------------------------------------------------------------------------
// MATA FUNCTION HEADERS, JUST FOR REFERENCE
// -------------------------------------------------------------------------------------------------
* void function initialize() 
* void function add_fe
* void function fe2local(`Integer' g) 
* void function prepare()
* `VarByFE' function transform(`VarByFE' indata, `Integer' g_from, `Integer' g_to) 
* `VarByFE' function count_by_group(`Series' group, `Series' index, | `Series' v, `Series' sorted_weight)
* `VarByFE' function mean_by_group(`Series' indata, `Series' index, `VarByFE' sum_count, `VarByFE' counts_to, `Series' * sorted_weight)
* `VarByFE' function remean_by_group(`VarByFE' indata, `Series' index, `VarByFE' sum_count, `VarByFE' counts_to, `Series'*  sorted_weight)
* `Matrix' function regress_by_group(`Series' y, `Matrix' x, `Series' index, 
* `Matrix' function compute_invxx(`Matrix' x, `Series' index, `VarByFE' offset, `VarByFE' count, `Series' sorted_weight)
* void function make_residual(...)
// -------------------------------------------------------------------------------------------------

// -------------------------------------------------------------------------------------------------
// Shortcuts
// -------------------------------------------------------------------------------------------------
* local Varlist 		string scalar
* local Integer 		real scalar
* local VarByFE 		real colvector // Should be levels*1
* local Series		real colvector // Should be N*1
* local Matrix		real matrix
* local SharedData 	external struct FixedEffect vector
// -------------------------------------------------------------------------------------------------



clear all
pr drop _all
discard

pwd
clear all
qui include ../source/_mata/reghdfe.mata
qui include ../source/_common/Assert.ado
qui include ../source/_common/Debug.ado
qui include ../source/_hdfe/CheckCorrectOrder.ado
qui include ../source/_hdfe/Start.ado
	qui include ../source/_hdfe/ParseOneAbsvar.ado
qui include ../source/_hdfe/Precompute.ado
	qui include ../source/_hdfe/GenerateID.ado
qui include ../source/_hdfe/Demean.ado


// -------------------------------------------------------------------------------------------------
// Quick checks
// -------------------------------------------------------------------------------------------------
mata
group = 1\1\2\2\1\1\3
index = order(group, 1)
freqs = 1\1\1\1\1\1\30
sorted_freqs = freqs[index, 1]
ans = count_by_group(group, index, 0, sorted_freqs)
assert(all(ans==(4\2\30)))
end

* ssc install moremata
mata
x = 1.5 \ 2.3 \ 6.5 \ 4.2 \ 1.1 \ 4.6 \ 0.0
group = 1\1\2\2\1\1\3
freqs = 1\3\1\2\1\10\30

index = order(group, 1)
sorted_freqs = freqs[index, 1]

counts_to = count_by_group(group, index)
sum_count = quadrunningsum(counts_to)

counts_to = count_by_group(group, index, 0 , sorted_freqs)

// x , group , freqs
ans = mean_by_group(x, index, sum_count, counts_to, sorted_freqs)
//check = 3.915385 \ 5.35 \ 0
check = 3.7 \ 4.96666 \ 0
ans , check

assert(sum(abs(ans-check))<0.01)
end

// -------------------------------------------------------------------------------------------------

sysuse auto, clear
local variables price weight length 
areg `variables' i.trunk, a(turn)

mata: VERBOSE = 3
Start, absorb(turn trunk)
Precompute, keep(`variables')

*Demean, varlist(`variables')
foreach var of local variables {
	mata: make_residual("`var'", "resid", 1e-07, 100) // , 0, 1, -1, 1, .005, 20, 3, 6)
	replace `var' = resid
	drop resid
}

regress `variables', nocons
exit
