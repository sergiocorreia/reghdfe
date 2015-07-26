cscript "hdfe" adofile reghdfe

* Setup
	discard
	clear all
	set more off
	* cls

* Convenience: "Trim <size>" will trim e(b) and e(V)
	capture program drop TrimMatrix
	program define TrimMatrix, eclass
	args size q
		assert `size'>0
		matrix trim_b = e(b)
		matrix V = e(V)
		if ("`q'"!="") {
			matrix V = V * `q'
			matrix V2 = V
			ereturn repost V = V2
		}
		matrix trim_V = V
		matrix trim_b = trim_b[1, 1..`size']
		matrix trim_V = trim_V[1..`size',1..`size']

		ereturn matrix trim_b = trim_b
		ereturn matrix trim_V = trim_V
	end
	
* Create fake dataset
	sysuse auto
	bys turn: gen t = _n
	tsset turn t

* [TEST] Cluster
	local lhs price
	local rhs weight length
	local absvars turn
	local clustervar rep78
	fvunab tmp : `rhs'
	local K : list sizeof tmp

	drop if missing(rep)

	* 1. Run -reghdfe- as benchmark
	di as result "reghdfe `lhs' `rhs', absorb(`absvars') vce(cluster `clustervar')"
	reghdfe `lhs' `rhs', absorb(`absvars') vce(cluster `clustervar')
	TrimMatrix `K'
	matrix list e(trim_V)
	storedresults save benchmark e()

	* 2. Run -hdfe-
	hdfe `lhs' `rhs', abs(`absvars') clustervars(`clustervar') clear keepids
	rename __CL1__ `clustervar'
	local df_a = e(df_a)
	return list
	qui regress `lhs' `rhs' , vce(cluster `clustervar') nocons // this will have different DoF so of course different VCE
	local n = e(N)
	local wrong = `n' - e(df_m)
	local correct = `wrong' - `df_a' // - 1
	local q = `wrong' / `correct'
	di as text "wrong=<`wrong'> correct=<`correct'> q=<`q'>"
	TrimMatrix `K' `q'
	`e(cmd)'
	matrix list e(trim_V)

	storedresults compare benchmark e(), tol(1e-12) include( ///
		scalar: N rss df_r ///
		matrix: trim_b trim_V ///
		macros: wexp wtype)

	hdfe `lhs' `rhs', abs(`absvars') clustervars(`clustervar') gen(resid_)
	qui reg resid_* , vce(cluster `clustervar') nocons
		local n = e(N)
		local wrong = `n' - e(df_m)
		local correct = `wrong' - `df_a' // - 1
		local q = `wrong' / `correct'
	TrimMatrix `K' `q'
	`e(cmd)'

	storedresults compare benchmark e(), tol(1e-12) include( ///
		scalar: N rss df_r ///
		matrix: trim_b trim_V ///
		macros: wexp wtype)

	storedresults drop benchmark

* syntax()
	sysuse auto, clear
	hdfe price weight length, a(turn trunk) gen(R_) sample(smpl) clusterv(trunk)
	hdfe price weight length, a(turn trunk) clear keepv(make) keepids clusterv(trunk)
	
	sysuse auto, clear
	hdfe price gear [w=weight], a(turn trunk) clear keepv(make) keepids clusterv(trunk)

cd "D:/Github/reghdfe/test"
exit
