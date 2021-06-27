/*

Programmer options

reghdfe, nopartial
reghdfe, noregress
reghdfe, keepmata
reghdfe, fastregress

asd

*/


sysuse auto
reghdfe price weight length gear i.foreign, a(turn trunk) vce(cluster turn) noregress
// validate that subsequent calls have the same options (absorb, vce, etc etc)
// best way to do this is to prohibit any options in the subsequent call, and just allow y x .. , cache(use)
// Not even IF/IN!!!

loc depvar weight
loc indepvars length gear i.foreign

tempvar touse
mata: HDFE.save_touse("`touse'", 1)

ms_expand_varlist `indepvars' if `touse'
return list
loc indepvars		"`r(varlist)'"
loc fullindepvars	"`r(fullvarlist)'"
loc not_omitted		"`r(not_omitted)'"

mata: HDFE.solution.depvar = "`depvar'"
mata: HDFE.solution.indepvars = tokens("`indepvars'")
mata: HDFE.solution.fullindepvars = tokens("`fullindepvars'")
mata: HDFE.solution.indepvar_status = !strtoreal(tokens("1 `not_omitted'")) // indepvar_status[i]=1 for variables omitted due to being a basevar (hack: the first element is the depvar)


*<<CONTINUE HERE
* unsure how the hdfe_data and hdfe_tss were used..

* Regress
if ("`keepmata'" != "") mata: hdfe_data = HDFE.solution.data // uses memory!
if ("`keepmata'" != "") mata: hdfe_tss = HDFE.solution.tss
if (`timeit') timer on 25
mata: reghdfe_solve_ols(HDFE, HDFE.solution, "vce_small")
if (`timeit') timer off 25
mata: HDFE.solution.cmdline = HDFE.solution.cmd + " " + st_local("0")

* Restore
if (`compact') {
	if (`verbose' > 0) di as text "## Restoring dataset"
	restore
}

if (`verbose' > 0) di as text "{title:Posting results to e() and displaying them}" _n

* Post regression results
* 1) Expand 'b' and 'V' to add base/omitted vars; expand fullindepvars/indepvars to add _cons
tempname b V
mata: HDFE.solution.expand_results("`b'", "`V'", HDFE.verbose)

* 2) Run "ereturn post"
loc store_sample = ("`nosample'"=="")
EreturnPost `touse' `b' `V' `store_sample'


*...
