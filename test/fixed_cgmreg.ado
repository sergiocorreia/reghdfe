*!Version 3.0.0 5Aug09 (By Jonah B. Gelbach)
*!Version 2.3.1 4Aug09 (By Jonah B. Gelbach)
*!Version 2.3.0 5Feb07 (By Jonah B. Gelbach)
*!Version 2.2.0 24Jan07 (By Jonah B. Gelbach)
*!Version 2.1.0 19Sep06 (By Jonah B. Gelbach)
*!Version 2.0.1 	(By Douglas L. Miller)
*!Version 2.0.0 22May06 (By Jonah B. Gelbach)
*!Version 1.0.1 22May06 (By Jonah B. Gelbach)
*!Version 1.0.0 28Mar06 (By Jonah B. Gelbach)



*************
* CHANGELOG *
*************

*
* 3.0.0: 
*	 
*	 ## Fixed minor issue with weights in ereturn post line, added e(wtype) local macro
*	 ## Added eigenvalue fix for non-psd cases
*	 ## Fixed ereturn behavior for scalars & macros
*
* 2.3.1: small edit by JBG to
*
*		## ensure that observations with missing cluster values are dropped
*		   (this matters b/c cgmreg runs Stata's -regress- without clustering, 
*		    so previous behavior was to include obs with missing cluster values
*		    and then treat "missing" as a cluster in its own right)
*
* 2.3.0: medium edit by JBG to 
*
*		## add treatment of if & in conditions
*		## add treatment of weights
*	
*        (required edit of syntax of to sub_robust subroutine, as well as adding some code on main regress line)
*
*
* 2.2.0: medium edit by JBG to make sure that "robust" option doesn't get passed to regress for line where we obtain (X'X)^(-1) using mse1 option.
*	 (comment: this seems like a stata bug to me -- why should stata allow you to use the robust option when the whole point is to get (X'X)^(-1)????
*
* 2.1.0: medium edit by JBG to move from use of -tab- to -unique- (I just dumped in the text of unique.ado to address this locally)
*
* 2.0.1: minor edit by Doug to unabbreviate "pred"
*
* 2.0.0: major addition: command now handles arbitrary number of grouping vars
*	 		 also, we now use cgmreg to calculate manually when only one clustering variable is used.
*			       this feature helps show that the sub_robust routine is correct

* 1.0.1: corrected error in 1.0.0:
*	I forgot to subtract out the estimate with cluster(bothvars) when `numcvars'==2


*********************
* SCHEMATIC OF CODE *
*********************

/*

	1. Run the regression, with mse1 option (this option sets robust off and variance=1, so that resulting "cov" matrix is (X'X)^-1)

	2. Save some matrices, scalars, and macros from the ereturned regress output

	3. Generate predicted residuals

	4. Set up a giant matrix that has

		* one column for every clustering variable
		* [(2^K) - 1] rows, where K is the number of clustering variables
		* elements equal to either 0 or 1

	   Each row of this matrix corresponds to one subcase, so that it provides a list of clustering vars for which we will calculate the middle matrix.

	   We then add or subtract the middle matrices according to the inclusion/exclusion rule: 

		* when the number of included clustering vars is odd, we add
		* when the number of included clustering vars is even, we subtract

	5. We then iterate over the rows of this matrix, using egen, group() and the list of included clustering variables 
	   to create a grouped var indicating an observation's membership category according to this group

	6. We then use the sub_robust subroutine (which uses _robust) to calculate the appropriate part of the covariance matrix

	7. The resulting part of the covariance matrix is added/subtracted to the matrix `running_mat', as given by the inc/exc rule

	8. The header in the stata output tells us

				Number of obs		[Total number of included observations]
				Num clusvars 		[Number of clustering vars, i.e., dimensions]
				Num combinations	[Total number of possible combinations of the clusvars, i.e., 2^K-1]

	   Followed by a list of the number of distinct categories in each clustering variable.

	9. Then the regression output appears, and we are done.

*/


program define fixed_cgmreg, eclass byable(onecall) sortpreserve

	syntax anything [if] [in] [aweight fweight iweight pweight /], /*
		*/ Cluster(string) [NOEIGenfix *] 

	*NOTE: use "NOEIGenfix" rather than "noEIGenfix" b/c we define a separate macro eigenfix below

	*marksample code added in version 2.3.1, replacing homemade mark that happened after regress
	marksample touse
	markout `touse' `cluster', strok

	local numcvars : word count `cluster'

	di
	while ( regexm("`options'","robust")==1 ) {

		di " -> Removing string 'robust' from your options line: it's unnecessary as an option,"
		di "    but it can cause problems if we leave it in."
		di "    If some variable in your options list contains the string 'robust', you will"
		di "    have to rename it."
		di 
		local options = regexr("`options'", "robust", "")

	} 

	/* deal with weights */
	if "`weight'"~="" {
		local wtype "`weight'"
		local weight "[`weight'=`exp']"
	} 
	else {
		local weight ""
	}

	/* main regression */
	qui regress `anything' if `touse' `weight', `options' mse1
	di "Note: +/- means the corresponding matrix is added/subtracted"
	di

/*commented for version 2.3.1 fix of missing cluster values issue (we make `touse' with marksample, above)
	/* copy some information that regress provides */
	tempvar touse
	qui gen `touse' = e(sample)		/* note that this will take care of zero-weight cases */
end commented for version 2.3.1 fix of missing cluster values issue
*/

	tempname b
	mat `b' = e(b)

	local depname = e(depvar)

	tempname N_clust N df_m F r2 mss rss r2_a ll ll_0 
*	scalar N_clust =	e(N_clust)
	scalar `N' =           	e(N) 
	scalar `df_m' =         e(df_m)
*	scalar `df_r' =         e(df_r)
*	scalar `F' =            e(F)
	scalar `r2' =           e(r2)
*	scalar `rmse' =         e(rmse)
	scalar `mss' =          e(mss)
	scalar `rss' =          e(rss)
	scalar `r2_a' =         e(r2_a)
*	scalar `ll' =           e(ll)
*	scalar `ll_0' =         e(ll_0)
	
	local title =           e(title)
	local depvar =          `depname'
	local cmd =             "cgmreg"
	local properties =      e(properties)
	local predict =         e(predict)
	local model =           e(model)
	local estat_cmd =       e(estat_cmd)
	local vcetype =         e(vcetype)
	local clustvar =        "`cluster'"
	local clusvar  =        "`cluster'"


	/* generate the residuals */
	tempvar resid
	qui predict double `resid' if `touse'==1, residual
	local n = e(N)
	

	*save (x'x)^-1
	tempname xxinv rows
	mat `xxinv' = e(V)
	mat `rows' = rowsof(e(V))
	local rows = `rows'[1,1]
	local cols = `rows'		/* avoid confusion */
	local k = e(df_m) + 1 // plus constant
	di as error "[cgmreg] rows=`rows' k=`k'"

	/* matrix that holds the running sum of covariance matrices as we go through clustering subsets */
	tempname running_sum
	mat `running_sum' = J(`rows',`cols',0)

	/* we will use a_cluster for matrix naming below as our trick to enumerate all clustering combinations */
	tempname Bigmat
	mat `Bigmat' = J(1,1,1)

	*taking inductive approach
	forvalues a=2/`numcvars' { /* inductive loop for Bigmat */

		mat `Bigmat' = J(1,`a',0) \ ( J(2^(`a'-1)-1,1,1) , `Bigmat' ) \ (J(2^(`a'-1)-1,1,0) , `Bigmat' ) 
		mat `Bigmat'[1,1] = 1

	} /* end inductive loop for Bigmat */

	mat colnames `Bigmat' = `cluster'

	local numsubs = 2^`numcvars' - 1
	local S = `numsubs' 			/* for convenience below */

	forvalues s=1/`S' { /* loop over rows of `Bigmat' */

		{	/* initializing */
			local included=0
			local grouplist
		} /* done initializing */

		foreach clusvar in `cluster' { /* checking whether each `clusvar' is included in row `s' of `Bigmat' */

			tempname element
			mat `element' = `Bigmat'[`s',"`clusvar'"] 
			local element = `element'[1,1]


			if `element' == 1 { /* add `clusvar' to grouplist if it's included in row `s' of `Bigmat' */

				local included= `included' + 1
				local grouplist "`grouplist' `clusvar'"

			} /* end add `clusvar' to grouplist if it's included in row `s' of `Bigmat' */
		} /* checking whether each `clusvar' is included in row `s' of `Bigmat' */


		*now we use egen to create the var that groups observations by the clusvars in `grouplist'
		tempname groupvar
		qui egen `groupvar' = group(`grouplist') if `touse'

		*now we get the robust estimate
		local plusminus "+"
		if mod(`included',2)==0 { /* even number */
			local plusminus "-"
		} /* end even number */

                sub_robust `if' `in' `weight', groupvar(`groupvar') xxinv(`xxinv') plusminus(`plusminus') resid(`resid') running_sum(`running_sum') touse(`touse') k(`k')
	
		di "Calculating cov part for variables: `grouplist' (`plusminus')"

	} /* end loop over rows of `Bigmat' */
	

	*checking/fixing non-psd variance estimate
	tempname eigenvalues eigenvectors
	*use mata to get eigenvalues after ensuring that variance matrix is (numerically) symmetric
	mata { 
		B = st_matrix("`running_sum'") 
		A = makesymmetric(B) 
		symeigensystem(A, C=., lamda=.) 
  		st_matrix("`eigenvalues'", lamda) 
		st_matrix("`eigenvectors'", C)
	}

	local rnames  : rownames `running_sum'
	local numcols = colsof(`running_sum')
	local eigenfix "no"
	forvalues col=1/`numcols' { /* column number loop */
		if (`eigenvalues'[1,`col']<0) {

			if "`noeigenfix'"=="noeigenfix" {
			    	di
			    	di " -> NOTE: Raw estimated variance matrix was non positive semi-definite."
			    	di
				di "          Because you used the -noeigenfix- option, -cgmreg- must end."
				di
		    		di "          See Cameron, Gelbach & Miller, "
		    		di "            'Robust Inference with Multi-Way Clustering'."
		    		di	
				di "Program terminated."			
				di
				exit
			}

		    	mat `eigenvalues'[1,`col']=0
		    	local eigenfix "yes"
		}
	} /* end column number loop */

	*now reconstruct variance matrix using spectral decomposition formula (e.g., Def A.16 in Greene, 6th)
	tempname raw_running_sum
	mat `raw_running_sum' = `running_sum'	/* pre eigen-fix variance matrix */
	mat `running_sum' = `eigenvectors'*diag(`eigenvalues')*`eigenvectors''
	mat rownames `running_sum' = `rnames'
	mat colnames `running_sum' = `rnames'
	/* end checking/fixing non-psd variance estimate */


	/* final cleanup and post */
	di
	di _column(50) "Number of obs     =    `n'"
	di _column(50) "Num clusvars      =    `numcvars'"
	di _column(50) "Num combinations  =    `S'"

	if "`if'"~="" di _column(50) "If condition      =    `if'"
	if "`in'"~="" di _column(50)     "In condition      =    `in'"
	if "`weight'"~="" di _column(50) "Weights are       =    `weight'"
	di
	local c 0
	foreach clusvar in `cluster' { /* getting num clusters by cluster var */

		local c = `c' + 1
		qui unique `clusvar' if `touse'
		di _column(50) "G(`clusvar')" _column(68) "=    " _result(18)
		local Gclusvar`c' = _result(18)
	} /* end getting num obs by cluster var */
	di

	ereturn post `b' `running_sum' , e(`touse') depname(`depname') 
	ereturn display

	*scalars
	ereturn scalar N    =              `N'
	ereturn scalar df_m =              `df_m'
*	ereturn scalar df_r =              `df_r'
	ereturn scalar r2 =                `r2'
*	ereturn scalar rmse =              `rmse'
	ereturn scalar mss =               `mss'
	ereturn scalar rss =               `rss'
	ereturn scalar r2_a =              `r2_a'

	local c 0
	foreach clusvar in `cluster' { /* getting num clusters by cluster var */

		local c = `c' + 1
		ereturn scalar N_clus`c' =  `Gclusvar`c''
		ereturn scalar N_clus_`clusvar' =  `Gclusvar`c''
	}

	
	*locals

	ereturn local eigenfix  = 	  "`eigenfix'"
	ereturn local cmdline   = 	  "cgmreg `anything' `if' `in' `weight', cluster(`cluster')"
	ereturn local title =             "Linear regression with CGM covariance estimation"
	ereturn local depvar =            "`depname'"
	ereturn local cmd =               "cgmreg"
	ereturn local properties =        "`properties'"
	ereturn local predict =           "regres_p"
	ereturn local model =             "ols"
	ereturn local estat_cmd =         "regress_estat"
	ereturn local vcetype =           "cgm_robust"
	ereturn local wtype =             "`wtype'"
	ereturn local clustvar =          "`cluster'"
	ereturn local clusvar  =          "`cluster'"

	*matrices
	ereturn matrix rawcovmat = `raw_running_sum'



	if "`eigenfix'"=="yes" {
		    di
		    di " -> NOTE: Raw estimated variance matrix was non-positive semi-definite."
		    di "          -cgmreg- is replacing any/all negative eigenvalues with 0."
		    di
		    di "          See Cameron, Gelbach & Miller, "
		    di "            'Robust Inference with Multi-Way Clustering'."
		    di
		    di "          Raw, non-psd covariance estimate will be available "
		    di "            in e(rawcovmat)."
		    di
		    di "          (If you don't want this behavior, use the 'noeigenfix' option,"
		    di "            in which case -cgmreg- will throw an error)"
		    di
		    di
	}
end



prog define sub_robust

	syntax [if] [in] [aweight fweight iweight pweight /] , groupvar(string) xxinv(string) plusminus(string) resid(string) running_sum(string) touse(string) k(integer)

/*
	local cvar 		"`1'"	/* cluster var, to be fed to us as argument 1 */
	local xxinv 		"`2'"	/* xxinv estimate, to be fed to us as argument 2 */
	local plusminus 	"`3'"	/* whether to add or subtract to `running_sum', argument 3 */
	local resid 		"`4'"	/* name of tempvar with resids in it, arg 4 */
	local running_sum 	"`5'"	/* running_sum estimate, to be fed to us as argument 5 */
	local touse		"`6'"
*/	

	/* deal with weights */
	if "`weight'"~="" {
		local weight "[`weight'=`exp']"
	} 
	else {
		local weight ""
	}

	tempname rows
	*mat `rows' = rowsof(`xxinv')
	*local rows = `rows'[1,1]
	local rows = `k'

	cap mat drop `m'
	tempname m
	mat `m' = `xxinv'

	if "`if'"=="" local if "if 1"
	else          local if "`if' & `touse'"
	qui _robust `resid' `if' `in' `weight', v(`m') minus(`rows') cluster(`groupvar')
	mat `running_sum' = `running_sum' `plusminus' `m'

*	mat li `running_sum'
end
*end program sub_robust

*! version 1.1  mh 15/4/98  arb 20/8/98
*got this from http://fmwww.bc.edu/repec/bocode/u/unique.ado
program define unique
local options "BY(string) GENerate(string) Detail"
local varlist "req ex min(1)"
local if "opt"
local in "opt"
parse "`*'"
tempvar uniq recnum count touse
local sort : sortedby
mark `touse' `if' `in'
qui gen `recnum' = _n
sort `varlist'
summ `touse', meanonly
local N = _result(18)
sort `varlist' `touse'
qui by `varlist': gen byte `uniq' = (`touse' & _n==_N)
qui summ `uniq'
di in gr "Number of unique values of `varlist' is  " in ye _result(18)
di in gr "Number of records is  "in ye "`N'"
if "`detail'" != "" {
	sort `by' `varlist' `touse'
	qui by `by' `varlist' `touse': gen int `count' = _N if _n == 1
	label var `count' "Records per `varlist'"
	if "`by'" == "" {
		summ `count' if `touse', d
	}
	else {
		by `by': summ `count' if `touse', d
	}
}
if "`by'" !="" {
	if "`generate'"=="" {
		cap drop _Unique
		local generat _Unique
	}
	else {
		confirm new var `generate'
	}

        drop `uniq'
	sort `by' `varlist' `touse'
	qui by `by' `varlist': gen byte `uniq' = (`touse' & _n==_N)
	qui by `by': replace `uniq' = sum(`uniq')
	qui by `by': gen `generate' = `uniq'[_N] if _n==1
	di in blu "variable `generate' contains number of unique values of `varlist' by `by'"
	list `by' `generate' if `generate'!=., noobs nodisplay
}
sort `sort' `recnum'
end
