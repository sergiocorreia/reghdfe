* TLDR:
* As [e(N)-e(df_a)] goes to infinity, the fraction [N-L]/[N-L-k] goes to 1 and R2a->R2
* Thus, we can't let N grow too much
* If the panel is balanced, then (N-L)/N corrects the xtreg_fe adjusted R2
* If N=LT, then the above is (T-1)/T, so the bias is larger for small T
* The bias seems small if e.g. R2 is 0.5, but is more relevant for smaller R2s
* This is in turn caused by larger e's compared to xb

cap cls
clear all

local N 10
local T 3
local NT = `N' * `T'
set obs `NT'

gen id = int( (_n-1) / `N' )
bys id: gen t = _n
tab id
tsset id t

gen double wrong = .
gen double correct = .

set obs 100

forv i=1/100 {
	cap drop a e x y
	di "." _c
	* y_it = b x_it + a_i + e_it
	qui gen a = rnormal() in 1/`NT'
	qui by id: replace a = a[1]

	qui gen e = 9 * rnormal() in 1/`NT'
	qui gen x = rnormal() in 1/`NT'
	qui gen y = x + a + e in 1/`NT'

	qui xtreg y x in 1/`NT', fe
	*di e(r2)
	*di e(r2_a)
	qui replace wrong = e(r2_a) in `i'

	qui reghdfe y x in 1/`NT', a(id) fast
	*di e(r2_within)
	*di e(r2_a_within)
	qui replace correct = e(r2_a_within) in `i'
}

su wrong correct, d
