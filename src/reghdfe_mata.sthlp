{smcl}
{* *! version 4.4.0 11sep2017}{...}
{vieweralsosee "reghdfe" "help reghdfe"}{...}
{vieweralsosee "ftools" "help ftools"}{...}
{viewerjumpto "Syntax" "ftools##syntax"}{...}
{viewerjumpto "Creation" "ftools##creation"}{...}
{viewerjumpto "Properties and methods" "ftools##properties"}{...}
{viewerjumpto "Description" "ftools##description"}{...}
{viewerjumpto "Usage" "ftools##usage"}{...}
{viewerjumpto "Example" "ftools##example"}{...}
{viewerjumpto "Remarks" "ftools##remarks"}{...}
{viewerjumpto "Using functions from collapse" "ftools##collapse"}{...}
{viewerjumpto "Experimental/advanced" "ftools##experimental"}{...}
{viewerjumpto "Source code" "ftools##source"}{...}
{viewerjumpto "Author" "ftools##contact"}{...}

{title:Title}

{p2colset 5 22 22 2}{...}
{p2col :{cmd:FixedEffects} {hline 2}}Mata class behind {cmd:reghdfe}{p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{pstd}
{it}To construct the object:


{p 8 16 2}
{it:class FixedEffects}
{cmd:fixed_effects(}{space 1}{it:absvars} [
{cmd:,}
{it:touse}{cmd:,} 
{it:weighttype}{cmd:,} 
{it:weightvar}{cmd:,} 
{it:drop_singletons}{cmd:,} 
{it:verbose}]{cmd:)}

{marker arguments}{...}
{synoptset 38 tabbed}{...}
{synopthdr}
{synoptline}
{p2coldent:* {it:string} absvars}names of variables that identify each set of fixed effects{p_end}
{synopt:{it:string} touse}name of dummy {help mark:touse} variable{p_end}
{synopt:{it:string} weighttype}type of weight (fweight, pweight, aweight, iweight){p_end}
{synopt:{it:string} weightvar}name of weight variable{p_end}
{synopt:{it:string} drop_singletons}whether to drop singleton groups or not{p_end}
{synopt:{it:string} verbose}how much information to report
(0: report warnings, 1 to 4 reports more details, -1 is silent){p_end}
{p2colreset}{...}


{marker usage}{...}
{title:Standard usage}

{pstd}(optional) First, you can declare the FixedEffects object:

{p 8 8 2}
{cmd:class FixedEffects}{it: HDFE}{break}

{pstd}Then, you create the object from categorical variables, categorical-continuous interactions, etc.:

{p 8 8 2}
{it:HDFE }{cmd:=}{bind: }{cmd:fixed_effects(}{it:varnames}{cmd:)}

{pstd}
Then you can modify the object and add important properties:

{p 8 8 2}{it:HDFE.varlist }{cmd:=}{bind: }{it:varlist} // used to report messages about all demeaned variables{p_end}
{p 8 8 2}{it:HDFE.indepvars }{cmd:=}{bind: }{it:indepvars} // used to report messages about demeaned regressors{p_end}
{p 8 8 2}{it:HDFE.num_clusters }{cmd:=}{bind: }{it:#} // Number of clusters{p_end}

{p 8 8 2}
{it: ... see reghdfe.ado for more options and how to combine them}


{marker properties}{...}
{title:Properties and Methods}

{marker arguments}{...}
{synoptset 38 tabbed}{...}

{synopthdr:properties (factors)}
{synoptline}

{synopt:{it:Integer} {cmd:N}}number of obs{p_end}
{synopt:{it:Integer} {cmd:M}}Sum of all possible FE coefs{p_end}
{synopt:{it:Factors} {cmd:factors}}{p_end}
{synopt:{it:Vector} {cmd:sample}}{p_end}
{synopt:{it:Varlist} {cmd:absvars}}{p_end}
{synopt:{it:Varlist} {cmd:ivars}}{p_end}
{synopt:{it:Varlist} {cmd:cvars}}{p_end}
{synopt:{it:Boolean} {cmd:has_intercept}}{p_end}
{synopt:{it:RowVector} {cmd:intercepts}}{p_end}
{synopt:{it:RowVector} {cmd:num_slopes}}{p_end}
{synopt:{it:Integer} {cmd:num_singletons}}{p_end}
{synopt:{it:Boolean} {cmd:save_any_fe}}{p_end}
{synopt:{it:Boolean} {cmd:save_all_fe}}{p_end}
{synopt:{it:Varlist} {cmd:targets}}{p_end}
{synopt:{it:RowVector} {cmd:save_fe}}{p_end}

{synopthdr:properties (optimization options)}
{synoptline}

{synopt:{it:Real} {cmd:tolerance}}{p_end}
{synopt:{it:Integer} {cmd:maxiter}}{p_end}
{synopt:{it:String} {cmd:transform}}Kaczmarz Cimmino Symmetric_kaczmarz (k c s){p_end}
{synopt:{it:String} {cmd:acceleration}}Acceleration method. None/No/Empty is none\{p_end}
{synopt:{it:Integer} {cmd:accel_start}}Iteration where we start to accelerate /set it at 6? 2?3?{p_end}
{synopt:{it:string} {cmd:slope_method}}{p_end}
{synopt:{it:Boolean} {cmd:prune}}Whether to recursively prune degree-1 edges{p_end}
{synopt:{it:Boolean} {cmd:abort}}Raise error if convergence failed?{p_end}
{synopt:{it:Integer} {cmd:accel_freq}}Specific to Aitken's acceleration{p_end}
{synopt:{it:Boolean} {cmd:storing_alphas}}1 if we should compute the alphas/fes{p_end}
{synopt:{it:Real} {cmd:conlim}}specific to LSMR{p_end}
{synopt:{it:Real} {cmd:btol}}specific to LSMR{p_end}

{synopthdr:properties (optimization objects)}
{synoptline}


{synopt:{it:BipartiteGraph} {cmd:bg}}Used when pruning 1-core vertices{p_end}
{synopt:{it:Vector} {cmd:pruned_weight}}temp. weight for the factors that were pruned{p_end}
{synopt:{it:Integer} {cmd:prune_g1}}Factor 1/2 in the bipartite subgraph that gets pruned{p_end}
{synopt:{it:Integer} {cmd:prune_g2}}Factor 2/2 in the bipartite subgraph that gets pruned{p_end}
{synopt:{it:Integer} {cmd:num_pruned}}Number of vertices (levels) that were pruned{p_end}

{synopthdr:properties (misc)}
{synoptline}

{synopt:{it:Integer} {cmd:verbose}}{p_end}
{synopt:{it:Boolean} {cmd:timeit}}{p_end}
{synopt:{it:Boolean} {cmd:store_sample}}{p_end}
{synopt:{it:Real} {cmd:finite_condition}}{p_end}
{synopt:{it:Real} {cmd:compute_rre}}Relative residual error: || e_k - e || / || e ||{p_end}
{synopt:{it:Real} {cmd:rre_depvar_norm}}{p_end}
{synopt:{it:Vector} {cmd:rre_varname}}{p_end}
{synopt:{it:Vector} {cmd:rre_true_residual}}{p_end}

{synopthdr:properties (weight-specific)}
{synoptline}

{synopt:{it:Boolean} {cmd:has_weights}}{p_end}
{synopt:{it:Variable} {cmd:weight}}unsorted weight{p_end}
{synopt:{it:String} {cmd:weight_var}}Weighting variable{p_end}
{synopt:{it:String} {cmd:weight_type}}Weight type (pw, fw, etc){p_end}

{synopthdr:properties (absorbed degrees-of-freedom computations)}
{synoptline}

{synopt:{it:Integer} {cmd:G_extended}}Number of intercepts plus slopes{p_end}
{synopt:{it:Integer} {cmd:df_a_redundant}}e(mobility){p_end}
{synopt:{it:Integer} {cmd:df_a_initial}}{p_end}
{synopt:{it:Integer} {cmd:df_a}}df_a_inital - df_a_redundant{p_end}
{synopt:{it:Vector} {cmd:doflist_M}}{p_end}
{synopt:{it:Vector} {cmd:doflist_K}}{p_end}
{synopt:{it:Vector} {cmd:doflist_M_is_exact}}{p_end}
{synopt:{it:Vector} {cmd:doflist_M_is_nested}}{p_end}
{synopt:{it:Vector} {cmd:is_slope}}{p_end}
{synopt:{it:Integer} {cmd:df_a_nested}}Redundant due to bein nested; used for: r2_a r2_a_within rmse{p_end}

{synopthdr:properties (VCE and cluster variables)}
{synoptline}

{synopt:{it:String} {cmd:vcetype}}{p_end}
{synopt:{it:Integer} {cmd:num_clusters}}{p_end}
{synopt:{it:Varlist} {cmd:clustervars}}{p_end}
{synopt:{it:Varlist} {cmd:base_clustervars}}{p_end}
{synopt:{it:String} {cmd:vceextra}}{p_end}

{synopthdr:properties (regression-specific)}
{synoptline}

{synopt:{it:String} {cmd:varlist}}y x1 x2 x3 x4 z1 z2 z3{p_end}
{synopt:{it:String} {cmd:depvar}}y{p_end}
{synopt:{it:String} {cmd:indepvars}}x1 x2{p_end}
    
{synopt:{it:Boolean} {cmd:drop_singletons}}{p_end}
{synopt:{it:String} {cmd:absorb}}contents of absorb(){p_end}
{synopt:{it:String} {cmd:select_if}}If condition{p_end}
{synopt:{it:String} {cmd:select_in}}In condition{p_end}
{synopt:{it:String} {cmd:model}}ols, iv{p_end}
{synopt:{it:String} {cmd:summarize_stats}}{p_end}
{synopt:{it:Boolean} {cmd:summarize_quietly}}{p_end}
{synopt:{it:StringRowVector} {cmd:dofadjustments}}firstpair pairwise cluster continuous{p_end}
{synopt:{it:Varname} {cmd:groupvar}}{p_end}
{synopt:{it:String} {cmd:residuals}}{p_end}
{synopt:{it:RowVector} {cmd:kept}}1 if the regressors are not deemed as omitted (by partial_out+cholsolve+invsym){p_end}
{synopt:{it:String} {cmd:diopts}}{p_end}

{synopthdr:properties (output)}
{synoptline}

{synopt:{it:String} {cmd:cmdline}}{p_end}
{synopt:{it:String} {cmd:subcmd}}{p_end}
{synopt:{it:String} {cmd:title}}{p_end}
{synopt:{it:Boolean} {cmd:converged}}{p_end}
{synopt:{it:Integer} {cmd:iteration_count}}e(ic){p_end}
{synopt:{it:Varlist} {cmd:extended_absvars}}{p_end}
{synopt:{it:String} {cmd:notes}}{p_end}
{synopt:{it:Integer} {cmd:df_r}}{p_end}
{synopt:{it:Integer} {cmd:df_m}}{p_end}
{synopt:{it:Integer} {cmd:N_clust}}{p_end}
{synopt:{it:Integer} {cmd:N_clust_list}}{p_end}
{synopt:{it:Real} {cmd:rss}}{p_end}
{synopt:{it:Real} {cmd:rmse}}{p_end}
{synopt:{it:Real} {cmd:F}}{p_end}
{synopt:{it:Real} {cmd:tss}}{p_end}
{synopt:{it:Real} {cmd:tss_within}}{p_end}
{synopt:{it:Real} {cmd:sumweights}}{p_end}
{synopt:{it:Real} {cmd:r2}}{p_end}
{synopt:{it:Real} {cmd:r2_within}}{p_end}
{synopt:{it:Real} {cmd:r2_a}}{p_end}
{synopt:{it:Real} {cmd:r2_a_within}}{p_end}
{synopt:{it:Real} {cmd:ll}}{p_end}
{synopt:{it:Real} {cmd:ll_0}}{p_end}

{synopthdr:methods}
{synoptline}

{synopt:{it:Void} {cmd:update_sorted_weights}()}{p_end}
{synopt:{it:Matrix} {cmd:partial_out}()}{p_end}
{synopt:{it:Void} {cmd:_partial_out}()}in-place alternative to {cmd:partial_out()}{p_end}
{synopt:{it:Variables} {cmd:project_one_fe}()}{p_end}
{synopt:{it:Void} {cmd:prune_1core}()}{p_end}
{synopt:{it:Void} {cmd:_expand_1core}()}{p_end}
{synopt:{it:Void} {cmd:estimate_dof}()}{p_end}
{synopt:{it:Void} {cmd:estimate_cond}()}{p_end}
{synopt:{it:Void} {cmd:save_touse}()}{p_end}
{synopt:{it:Void} {cmd:store_alphas}()}{p_end}
{synopt:{it:Void} {cmd:save_variable}()}{p_end}
{synopt:{it:Void} {cmd:post_footnote}()}{p_end}
{synopt:{it:Void} {cmd:post}()}{p_end}
{synopt:{it:Void} {cmd:reload}(copy=0)}{p_end} (run this if e.g. touse changes)

{synopthdr:methods (LSMR-specific)}
{synoptline}

{synopt:{it:Real} {cmd:lsmr_norm}()}{p_end}
{synopt:{it:Vector} {cmd:lsmr_A_mult}()}{p_end}
{synopt:{it:Vector} {cmd:lsmr_At_mult}()}{p_end}


{marker functions}{...}
{title:Additional functions}

{pstd}
Several useful Mata functions are included. For instance,

{p 8 16 2}
{it:void}
{cmd:reghdfe_solve_ols(}{it:HDFE}
{cmd:,}
{it:X}{cmd:,}
{it:...}
{cmd:)}

{pstd}
See {stata "mata: mata desc using lreghdfe"} for full list of functions and classes.


{marker description}{...}
{title:Description}

{pstd}
TBD


{marker example}{...}
{title:Example: OLS regression}

{pstd}
TBD


{inp}
    {hline 60}
	sysuse auto, clear
	local depvar price
	local indepvars weight gear
	mata: HDFE = fixed_effects("turn", "", "fweight", "trunk", 0, 2)
	mata: HDFE.varlist = "`depvar' `indepvars'"
	mata: HDFE.indepvars = "`indepvars'"
	mata: data = HDFE.partial_out("`depvar' `indepvars'")
	mata: reghdfe_solve_ols(HDFE, data, b=., V=., N=., rank=., df_r=., resid=., kept=., "vce_none")
	mata: b
    {hline 60}
{text}


{marker remarks}{...}
{title:Remarks}

{pstd}
TBD

{marker experimental}{...}
{title:Experimental/advanced functions}

{pstd}
TBD (LSMR, Prune, Bipartite?)

{marker source}{...}
{title:Source code}

{pstd}
{view reghdfe.mata, adopath asis:reghdfe.mata};
{view reghdfe_bipartite.mata, adopath asis:reghdfe_bipartite.mata};
{view reghdfe_class.mata, adopath asis:reghdfe_class.mata};
{view reghdfe_constructor.mata, adopath asis:reghdfe_constructor.mata};
{view reghdfe_common.mata, adopath asis:reghdfe_common.mata};
{view reghdfe_projections.mata, adopath asis:reghdfe_projections.mata};
{view reghdfe_transforms.mata, adopath asis:reghdfe_transforms.mata};
{view reghdfe_accelerations.mata, adopath asis:reghdfe_accelerations.mata};
{view reghdfe_lsmr.mata, adopath asis:reghdfe_lsmr.mata}
{p_end}

{pstd}
Also, the latest version is available online: {browse "https://github.com/sergiocorreia/reghdfe/tree/master/src"}


{marker author}{...}
{title:Author}

{pstd}Sergio Correia{break}
{break}
{browse "http://scorreia.com"}{break}
{browse "mailto:sergio.correia@gmail.com":sergio.correia@gmail.com}{break}
{p_end}


{marker project}{...}
{title:More Information}

{pstd}{break}
To report bugs, contribute, ask for help, etc. please see the project URL in Github:{break}
{browse "https://github.com/sergiocorreia/reghdfe"}{break}
{p_end}


{marker acknowledgment}{...}
{title:Acknowledgment}

{pstd}
TBD
