noi cscript "reghdfe: allow TS in absorb; i.turn#c.L.gear" adofile reghdfe

* Dataset
	sysuse auto
	bys foreign: gen t = _n
	xtset foreign t

	local included_e ///
		scalar: N rmse tss rss mss r2 r2_a df_r  ll /// ll_0 F df_m
		matrix: trim_b trim_V ///
		macros: wexp wtype

	* Note: cannot test <ll_0  F  df_m> because they are in reghdfe's absorb but as main regressors in areg

* [TEST]

	local lhs price
	local rhs weight length
	local absvars turn foreign#c.(gear L(-1/1).rep)
	fvunab tmp : `rhs'
	local K : list sizeof tmp

	* 1. Run benchmark
	areg `lhs' `rhs' ibn.foreign#c.(gear L(-1/1).rep) , absorb(turn)
	di e(df_a)
	trim_cons `K'
	storedresults save benchmark e()

/*	gen x1 = (foreign==0)*(gear)
	gen x2 = (foreign==1)*(gear)
	gen x3 = (foreign==0)*(L.rep)
	gen x4 = (foreign==1)*(L.rep)
	gen x5 = (foreign==0)*(rep)
	gen x6 = (foreign==1)*(rep)
	gen x7 = (foreign==0)*(F.rep)
	gen x8 = (foreign==1)*(F.rep)
	areg price x* weight length , absorb(turn)*/


	* 2. Run reghdfe and compare
	reghdfe `lhs' `rhs', absorb(`absvars') keepsingletons verbose(-1)
	notrim

	storedresults compare benchmark e(), tol(1e-11) include(`included_e')

storedresults drop benchmark
exit


