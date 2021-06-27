* ===========================================================================
* Test that we can save+load our Mata objects without segfaults
* ===========================================================================

	clear all
	cls
	*cap ado uninstall reghdfe

	mata: mata set matastrict on

clear all

*do "reghdfe_debug.mata"
*include "reghdfe_debug.mata"
include "reghdfe.mata", adopath

mata: mata desc

use "../../test/part2/toy-patents-long", clear
bys patent_id: gen byte touse = _n == 1
gen byte indiv_touse = 1

cap erase "data0.tmp"
mata:
unlink_folder("C:\Git\asd\borrar/PARALLEL_1234", 1)

HDFE = FixedEffects()
HDFE.absvars = "inventor_id year"
HDFE.tousevar = "touse"
HDFE.weight_type = ""
HDFE.weight_var = ""
HDFE.technique = "lsmr"
HDFE.transform = "symmetric_kaczmarz"
HDFE.acceleration = "conjugate_gradient"
HDFE.preconditioner = "diagonal"
HDFE.parallel_dir = "C:\Git\asd\borrar/PARALLEL_1234"
HDFE.parallel_opts = `"maxproc(2) id(1234) cores(2) tmp_path("C:\Git\asd\borrar/") "'
HDFE.drop_singletons = 1


HDFE.group_id = `"patent_id"' 
HDFE.individual_id = `"inventor_id"' 
HDFE.indiv_tousevar = `"indiv_touse"' 
HDFE.function_individual = `"mean"' 
//HDFE.parallel_maxproc = 2

HDFE.tolerance = 1.00000000000e-08
HDFE.maxiter = 16000
HDFE.compact = 0
HDFE.poolsize = 10
HDFE.verbose = 0

HDFE.vcetype = "unadjusted"
HDFE.num_clusters = 0
HDFE.clustervars = tokens("")
HDFE.base_clustervars = tokens("")


//estimate_dof(HDFE, tokens("pairwise clusters continuous"), "")

HDFE.init()

HDFE.factors[1].is_individual_fe
HDFE.factors[1].absvar

HDFE.partial_out("citations funding lab_size", 1, 1)

cleanup_for_parallel(HDFE)

HDFE.factors[1].is_individual_fe
HDFE.factors[1].absvar


fn = "data0.tmp"
fh = fopen(fn, "w")
fputmatrix(fh, HDFE)
fclose(fh)
end


mata: mata drop HDFE

mata:
fh = fopen(fn, "r")
HDFE = fgetmatrix(fh)
fclose(fh)

HDFE.factors[1].varlist
HDFE.factors[2].varlist


end


exit


clear all
include "reghdfe.mata", adopath
mata: fn = "C:\Git\asd\borrar/data0.tmp"
mata: fh = fopen(fn, "r")
mata: HDFE = fgetmatrix(fh)
mata: fclose(fh)

exit


mata: mata desc
mata: mata drop HDFE

mata: fn = "C:\Git\asd\borrar/data0.tmp"
mata: fh = fopen(fn, "r")
mata: HDFE = fgetmatrix(fh)
mata: fclose(fh)

