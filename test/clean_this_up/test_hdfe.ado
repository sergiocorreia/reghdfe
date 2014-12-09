cap pr drop test_hdfe
pr test_hdfe

* Run and store results
	qui reghdfe `0'
	if !e(alternative_cmd_ok) di as error "Approx command only!"
	***ereturn list
	local rhs = "`e(indepvars)' `e(AvgE_Ws)' `e(endogvars)'"
	
	* Sanity checks
	Assert e(ll)==.
	Assert e(ll_0)==.
	
	* Store scalars
	local scalars N df_m mss rss rank df_r F df_a tss r2 r2_a rmse
	foreach scalar of local scalars {
		local `scalar' = e(`scalar')
	}
	

	
	* Store b & V without constant term (which is not identified wrt absorbed)
	matrix b = e(b)
	matrix V = e(V)
	matrix b = b[1, 1..`df_m']
	matrix V = V[1..`df_m', 1..`df_m']

	
* Run benchmark
	di as text "Benchmark: `e(alternative_cmd)'"
	qui reghdfe, alternative
	matrix bb = e(b)
	matrix VV = e(V)
	matrix bb = bb[1, 1..`df_m']
	matrix VV = VV[1..`df_m', 1..`df_m']
	
* Check same betas
	matrix eps = mreldif(b, bb)
	local eps = eps[1,1]
	Assert `eps'<1e-6 /* epsfloat() */ , msg("b or V not the same (eps=`eps')")
	
	***ereturn list

	* Check same scalars
	Assert e(df_m)==`df_m'+`df_a'

	
	local scalars N mss rss df_r F r2 r2_a rmse
	foreach scalar of local scalars {
		if "`scalar'"=="F" {
			di as text "Test: `rhs'"
			qui testparm `rhs'
			Assert abs(``scalar''/r(F)-1)<1e-4 , rc(88) msg("FTest")
			Assert `df_m'==`r(df)', rc(88) msg("FTest")
			Assert `df_r'==`r(df_r)', rc(88) msg("FTest")
		}
		else if (int(``scalar'')==``scalar'') {
			Assert ``scalar''==`e(`scalar')', msg("`scalar': ``scalar''=/=`e(`scalar')'") rc(88)
		}
		else if abs(``scalar'')<1.0 {
			Assert abs(``scalar''-`e(`scalar')')<1e-4, msg("`scalar': abs(``scalar''-`e(`scalar')')=`=abs(``scalar''-`e(`scalar')')' !<1e-4") rc(88)
		}
		else {
			Assert abs(``scalar''/`e(`scalar')'-1)<1e-4, msg("`scalar': abs(``scalar''/`e(`scalar')'-1)=`=abs(``scalar''/`e(`scalar')'-1)' !<1e-4") rc(88)
		}
	}

	Assert e(rank)==`rank'+`df_a'
	Assert e(df_m)==`df_m'+`df_a'


	noi di in ye "... test passed"
	matrix drop b V
end

exit


* Compare
	ereturn list
	
	exit
	
	
	ereturn list
	
	

	ereturn list

cap pr drop test_one
program test_one
syntax [, CUT]
	matrix bb = e(b)
	matrix VV = e(V)
	
	* Constant is meaningless in areg .. i.i , a(j)
	*di "<`cut'>"
	if ("`cut'"!="") {
		matrix bb = bb[1, 1..2]
		matrix VV = VV[1..2, 1..2]
	}
	*matrix list b
	*matrix list bb
	
	matrix eps = mreldif(b, bb)
	local eps = eps[1,1]
	assert `eps'<epsfloat()
	
	matrix eps = mreldif(V, VV)
	local eps = eps[1,1]
	assert `eps'<epsfloat()
	
	assert N==e(N)
	assert abs(F-e(F))<1e-4
	
	di as text "[HDFE1] Test passed"
end
