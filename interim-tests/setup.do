OBSOLETE
* ===========================================================================
* Install groupreg.ado
* ===========================================================================

	log close _all
	clear all
	discard
	cls
	pr drop _all

	cap ado uninstall ftools
	net install ftools, from("c:/git/ftools/src/")

	**reghdfe, reload // because we are updating ftools
	**pr drop _all
	**clear all // because reghdfe pollutes mata


	* Uninstall existing versions
	cap ado uninstall groupreg

	* Check that the .mata files don't have compile-time errors
	* EG: inspect by hand for the word "unused"
	set trace on
	set tracedepth 3
	
	cd "../src"
	do groupreg.mata
	cd "../test"

	* Reinstall
	* Note: "net install" requires hardcoded paths so we do a workaround
	mata: st_local("path", pathresolve(pwd(), "../src"))
	mata: assert(direxists("`path'"))
	net install groupreg , from("`path'")


	* Quick test that it works
	set trace off
	sysuse auto
	which groupreg

	groupreg price weight length, a(turn)
	groupreg price weight length, a(turn) noregress
	groupreg price weight length, a(turn) noregress v(1)
	groupreg price weight length, a(turn, save) tech(map)


	do alphas
	do test-error-sample

asd



	cls
	groupreg price weight, a(turn) precond(none) dof(none)
	*groupreg price weight , a(turn##c.gear_ratio trunk##c.(length displacement)) precond(none) pool(1)  tol(1e-12)  v(0) maxiter(1000) keepsing
	groupreg price weight , a(turn##c.gear_ratio trunk##c.(length displacement)) pool(1)  tol(1e-12) v(0) maxiter(1000) keepsing

	* Test that indiv(fe) group(_n) works the same as reghdfe
	*do test_trivial_group
	*stopit


	do "../test/test_simple-absorb.do"
	do "../test/test_aggregation.do"


exit

	do demo_indiv_preconditioner
	stopit

	do demo_groups
	stopit

	groupreg price weight , a(turn trunk)
	groupreg price weight , a(turn trunk) v(2)
	stopit2

	cls
	groupreg price weight , a(turn##c.gear_ratio trunk##c.(length displacement)) precond(none) pool(1)  tol(1e-12)  v(0) maxiter(1000) keepsing
	di e(r2)
	groupreg price weight , a(turn##c.gear_ratio trunk##c.(length displacement)) precond(diag) pool(1)  tol(1e-12)  v(0) maxiter(1000) keepsing
	di e(r2)
	groupreg price weight , a(turn##c.gear_ratio trunk##c.(length displacement)) precond(block) pool(1)  tol(1e-12)  v(0) maxiter(1000) keepsing
	di e(r2)
	qui reg price weight turn##c.gear_ratio trunk##c.(length displacement)
	di e(r2)
	qui reghdfe price weight , a(turn##c.gear_ratio trunk##c.(length displacement)) keepsing
	di e(r2)
	qui reghdfe price weight , a(turn##c.gear_ratio trunk##c.(length displacement)) tol(1e-12) keepsing
	di e(r2)
	stopit

	cls
	set varabbrev off
	cap drop w
	gen double w = runiform()
	//groupreg price weight [aw=w], a(turn) precond(none)
	//reghdfe price weight [aw=w], a(turn)
	//stopit
	//groupreg price weight [aw=w], a(turn trunk##c.(length displacement)) precond(none)
	//reghdfe price weight [aw=w], a(turn trunk##c.(length displacement))
	//stopit


	*keep in 44/60
	*keep in 1/30
	reghdfe price , a(turn##c.length) keepsing
	groupreg price , a(turn##c.length) precond(none) pool(1) keepsing tol(1e-12)
	groupreg price , a(turn##c.length) precond(diag) pool(1) keepsing tol(1e-12)

	//mata: mata matsave "object.tmp" HDFE, replace
	//mata: mata matdesc "object.tmp"
	//mata: mata desc
	//mata: mata matuse "object.tmp", replace

	//groupreg price , a(turn##c.length) precond(no) pool(1) keepsing tol(1e-12)
	//reghdfe  price , a(turn##c.length) keepsing

	* THESE DONT MATCH :( (unless we use 1e-12 tol)
	cls
	groupreg price weight , a(turn trunk##c.(length displacement)) precond(no) pool(1)  tol(1e-12)
	groupreg price weight , a(turn trunk##c.(length displacement)) precond(diag) pool(1)  tol(1e-12)
	groupreg price weight , a(turn trunk##c.(length displacement)) precond(block_diag) pool(1)  tol(1e-12)
	reghdfe price weight , a(turn trunk##c.(length displacement))

	//asd

	groupreg price weight , a(turn trunk##c.(length displacement) foreign#c.mpg) precond(no) pool(1) tol(1e-12)
	reghdfe  price weight , a(turn trunk##c.(length displacement) foreign#c.mpg) tol(1e-12)

	groupreg price weight, a(turn) pool(1) precond(none) keepsing
	reghdfe  price weight, a(turn) keepsing

	groupreg price weight, a(turn)
	groupreg price weight, a(turn trunk)
	groupreg price weight, a(turn)

	gen byte n = 1
	qui groupreg price weight length gear_ratio displacement [fw=n], a(turn)

	groupreg price weight gear_ratio length, a(turn##c.headroom) keepsing preconditioner(none) v(1)

	noi di as text "groupreg has been reinstalled from local source!"

exit
