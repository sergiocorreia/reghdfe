cscript "reghdfe with clusters" adofile reghdfe

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
	*gen n = int(uniform()*10+3) // used for weights
	*replace length = 0 if rep==3 // used for DoF adjustment of cont var
	*replace length = 5 if rep==1
	*gen byte one = 1
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
	reghdfe `lhs' `rhs', absorb(`absvars') vce(cluster `clustervar') nocons
	TrimMatrix `K'
	matrix list e(trim_V)
	storedresults save benchmark e()

	* 2. Run -hdfe-
	preserve
	hdfe `lhs' `rhs', abs(`absvars') clustervars(`clustervar') clear
	local df_a = r(df_a)
	return list
	qui regress `lhs' `rhs' , vce(cluster `clustervar') nocons // this will have different DoF so of course different VCE
	local n = e(N)
	local wrong = `n' - e(df_m)
	local correct = `wrong' - `df_a' - 1
	local q = `wrong' / `correct'
	di as text "wrong=<`wrong'> correct=<`correct'> q=<`q'>"
	TrimMatrix `K' `q'
	`e(cmd)'
	matrix list e(trim_V)
	restore

	storedresults compare benchmark e(), tol(1e-12) include( ///
		scalar: N rss df_r ///
		matrix: b V ///
		macros: wexp wtype)

	hdfe `lhs' `rhs', abs(`absvars') clustervars(`clustervar') gen(resid_)
	qui reg resid_* , vce(cluster `clustervar') nocons
		local n = e(N)
		local wrong = `n' - e(df_m)
		local correct = `wrong' - `df_a' - 1
		local q = `wrong' / `correct'
	TrimMatrix `K' `q'
	`e(cmd)'

	storedresults compare benchmark e(), tol(1e-12) include( ///
		scalar: N rss df_r ///
		matrix: b V ///
		macros: wexp wtype)

	storedresults drop benchmark

* partial()
	sysuse auto, clear
	areg price weight length gear disp i.trunk, absorb(turn)
	*reghdfe price weight length gear disp, a(turn trunk)
	TrimMatrix `K'
	storedresults save benchmark e()

	hdfe price weight length, a(turn trunk) partial(gear disp) clear
	*return list
	local df_a = r(df_a) + r(df_partial)

	qui regress price weight length
	local dof = e(df_r) - `df_a'
	regress price weight length, dof(`dof') nocons
	TrimMatrix `K'
	
	storedresults compare benchmark e(), tol(1e-12) include( ///
		scalar: N rss df_r ///
		matrix: trim_b trim_V ///
		macros: wexp wtype)

	storedresults drop benchmark

cd "D:/Github/reghdfe/test"
exit
