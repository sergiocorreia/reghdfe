clear all
cls
*cap ado uninstall reghdfe


/*include "ftools.mata", adopath

sysuse auto
cap erase "data0.tmp"
mata:
F = Factor(2)
F[1] = factor("turn")
F[2] = factor("trunk")
F[1].extra = .
F[1].vl = .
F[2].extra = .
F[2].vl = .
fn = "data0.tmp"
fh = fopen(fn, "w")
fputmatrix(fh, F)
fclose(fh)

fh = fopen(fn, "r")
F2 = fgetmatrix(fh)
fclose(fh)

F2[1].varlist

end

*exit
*/

clear all

*do "reghdfe_debug.mata"
*include "reghdfe_debug.mata"
include "reghdfe.mata", adopath

mata: mata desc

sysuse auto
gen byte touse = 1

cap erase "data0.tmp"
mata:
HDFE = FixedEffects()
HDFE.absvars = "turn"
HDFE.tousevar = "touse"
HDFE.weight_type = ""
HDFE.weight_var = ""
HDFE.technique = "lsmr"
HDFE.transform = "symmetric_kaczmarz"
HDFE.acceleration = "conjugate_gradient"
HDFE.preconditioner = "block_diagonal"
HDFE.parallel_dir = "C:\Git\asd\borrar/PARALLEL_374598141"
HDFE.parallel_opts = `"maxproc(2) id(374598141) tmp_path("C:\Git\asd\borrar/") "'
HDFE.drop_singletons = 1
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
HDFE.partial_out("price weight length", 1, 1)

cleanup_for_parallel(HDFE)

HDFE.factors[1].varlist

fn = "data0.tmp"
fh = fopen(fn, "w")
fputmatrix(fh, HDFE)
fclose(fh)


fh = fopen(fn, "r")
HDFE = fgetmatrix(fh)
HDFE.factors[1].varlist
fclose(fh)

end


exit


clear all
include "reghdfe.mata", adopath
mata: fn = "C:\Git\asd\borrar/data0.tmp"
mata: fh = fopen(fn, "r")
mata: HDFE = fgetmatrix(fh)
mata: fclose(fh)


exit
