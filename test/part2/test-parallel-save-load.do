* ===========================================================================
* Test that we can save+load our Mata objects without segfaults
* ===========================================================================

	clear all
	*cap ado uninstall reghdfe

	mata: mata set matastrict on

	include "reghdfe.mata", adopath


// --------------------------------------------------------------------------
// Validate that we can save and load Factor objects
// --------------------------------------------------------------------------
	
sysuse auto
cap erase "data0.tmp"
mata:
	F = Factor(2)
	F[1] = factor("turn")
	F[2] = factor("trunk")

	F[1].cleanup_before_saving()
	F[2].cleanup_before_saving()

	fn = "data0.tmp"
	fh = fopen(fn, "w")
	fputmatrix(fh, F)
	fclose(fh)

	fh = fopen(fn, "r")
	F2 = fgetmatrix(fh)
	fclose(fh)

	assert(F2[1].varlist == "turn")
end


// --------------------------------------------------------------------------
// Validate that we can save and load FE_Factor() objects
// --------------------------------------------------------------------------

cap erase "data0.tmp"
mata:
	F = fe_factor("turn")
	F.absvar = F.varlist
	F.cleanup_before_saving()

	fn = "data0.tmp"
	fh = fopen(fn, "w")
	fputmatrix(fh, F)
	fclose(fh)

	fh = fopen(fn, "r")
	F = fgetmatrix(fh)
	fclose(fh)

	assert(F.absvar == "turn")
end


// --------------------------------------------------------------------------
// Validate that we can save and load Bipartite() objects
// --------------------------------------------------------------------------

cap erase "data0.tmp"
use "../../test/part2/toy-patents-long", clear
bys patent_id: gen byte touse = _n == 1
gen byte indiv_touse = 1

mata:
	FG = factor("patent_id")
	FI = factor("inventor_id")
	bg = BipartiteGraph()

	bg.init(&FG, &FI, 1)

	bg.cleanup_before_saving()
	FG.cleanup_before_saving()
	FI.cleanup_before_saving()

	fn = "data0.tmp"
	fh = fopen(fn, "w")
	fputmatrix(fh, bg)
	fclose(fh)

	fh = fopen(fn, "r")
	bg = fgetmatrix(fh)
	fclose(fh)
end


// --------------------------------------------------------------------------
// Validate that we can save and load Individual_Factor() objects
// --------------------------------------------------------------------------

cap erase "data0.tmp"

use "../../test/part2/toy-patents-long", clear
bys patent_id: gen byte touse = _n == 1
gen byte indiv_touse = 1

mata:
	
	sample = selectindex(st_data(., "touse"))
	indiv_sample = selectindex(st_data(., "indiv_touse"))
	F = indiv_factor("patent_id", "inventor_id", sample, indiv_sample, "mean", 1)

	F.cleanup_before_saving()

	fn = "data0.tmp"
	fh = fopen(fn, "w")
	fputmatrix(fh, F)
	fclose(fh)

	fh = fopen(fn, "r")
	F = fgetmatrix(fh)
	fclose(fh)

	assert(F.aggregate_function == "mean")
	//assert(F.FG.varlist == "patent_id")
	//assert(F.FI.varlist == "inventor_id")

end


clear all

exit
