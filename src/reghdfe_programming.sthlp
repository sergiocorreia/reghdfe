{smcl}
{* *! version 6.12.0 26June2021}{...}
{vieweralsosee "reghdfe" "help reghdfe"}{...}
{vieweralsosee "ftools" "help ftools"}{...}
{vieweralsosee "" "--"}{...}
{vieweralsosee "ivreghdfe" "help ivreghdfe"}{...}
{vieweralsosee "ppmlhdfe" "help ppmlhdfe"}{...}
{vieweralsosee "sumhdfe" "help sumhdfe"}{...}
{vieweralsosee "" "--"}{...}
{viewerjumpto "1. Ancillary commands" "reghdfe_internals##description"}{...}
{viewerjumpto "2. Undocumented options" "reghdfe_internals##undocumented"}{...}
{viewerjumpto "3. FixedEffects Mata class" "reghdfe_internals##mata"}{...}
{viewerjumpto "Description" "reghdfe_internals##description"}{...}
{viewerjumpto "Syntax" "ftools##syntax"}{...}
{viewerjumpto "Creation" "ftools##creation"}{...}
{viewerjumpto "Properties and methods" "ftools##properties"}{...}
{viewerjumpto "Usage" "ftools##usage"}{...}
{viewerjumpto "Example" "ftools##example"}{...}
{viewerjumpto "Remarks" "ftools##remarks"}{...}
{viewerjumpto "Using functions from collapse" "ftools##collapse"}{...}
{viewerjumpto "Experimental/advanced" "ftools##experimental"}{...}
{viewerjumpto "Source code" "reghdfe##source"}{...}
{viewerjumpto "Author" "reghdfe##contact"}{...}
{title:Using reghdfe with other commands}

{pstd}
This help file describes how to use reghdfe to within other programs, either in Stata or Mata.
It discusses three types of tools that might be useful for developers:

{phang2}  1. Ancillary commands from {help ftools} that are used by reghdfe, such as {help ms_get_version}.{p_end}
{phang2}  2. Undocumented options of {help reghdfe}.{p_end}
{phang2}  3. The {stata viewsource reghdfe.mata:FixedEffects} Mata class behind reghdfe, which can be used to build efficient Mata estimation programs.{p_end}

{pstd}
These commands are nested in order of integration with reghdfe.
Someone writing a command independent of reghdfe might still benefit from #1.
Someone writing a command that calls reghdfe a few times, such as {help ivreghdfe}, {help sumhdfe}, or {help did_imputation} might benefit from #2.
And someone writing a command that calls reghdfe multiple times, such as {help ppmlhdfe}
might also be interested in #3, due to the increase in efficiency and Mata integration.


{marker ancillary}{...}
{title:1. Ancillary commands}


{pstd}
{bf:{ul:ms_get_version}}

{pstd}
It's possible your command will depend on other user-written commands, in the same way as {cmd:reghdfe} depends on {cmd:ftools}.
To ensure compatibility and reproducibility, you can use the {cmd:ms_get_version} command to ensure that users are not running versions of these programs that are not too old.
For instance, reghdfe version 6.12.3 requires ftools of at least version 2.49.1.

{pstd}
The syntax of {cmd:ms_get_version} is:

{p 8 15 2} {it:ms_get_version}
{it:command}
[{cmd:,} {opt min_version(str)} {opt min_date(str)}]
{p_end}

{pstd}
Note: this command has been superceded by {cmd:require}
{p_end}

{marker options_table}{...}
{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{synopt: {opt min_version(str)}}minimum version.
Only supports {browse "https://semver.org/":semantic versioning} of the form {it:x.y.z}.
Note that versions are defined in the {help which:first line} an ado-file,
and can be verified by typing {it:which <command>}.{p_end}
{synopt : {opt min_date(str)}}(less used) minimum date of the program.
Supports dates of the form 1jan2018 or 01Jan2018.{p_end}

{pstd}
{cmd:min_version} also stores the following local variables:

{synoptset 24 tabbed}{...}
{synopt:`version_number'}version of the requested program{p_end}
{synopt:`version_date'}date of the requested program{p_end}
{synopt:`package_version'}concatenation of `version_number' and `version_string'{p_end}

{pstd}
An sample usage of {it:ms_get_version}, currently used by reghdfe, is:

{phang2}{cmd: ms_get_version ftools, min_version("2.46.0")}{p_end}


{marker undocumented}{...}
{title:2. Undocumented reghdfe options}

{pstd}
Sometimes you might not want to run the entire reghdfe command, but stop at some point and only compute certain objects. There are several objects that allow this.

{pstd}
{bf:A) Compute HDFE Nata object but stop before partialling out variables}

{p 8 15 2} {cmd:reghdfe} {it:...} {cmd:,}
{opt nopartial:out} [{help reghdfe##options_table:options}]{p_end}

{phang2}
This step will parse all inputs and initialize the {it:HDFE} object of the {cmd:FixedEffects} class.
Note that although the regression variables (depvar and indepvar) are not processed, if they have missing values the sample will reflect that.

{phang2}
For instance, the sample ado-file below is enough to create a program that reports the number of
singletons in a regression, without having to actually computed:

{phang2}{hline 16} {cmd: show_singletons.ado} {hline 16}{p_end}
{phang2}{cmd: prog show_singletons}{p_end}
{p 12 12 2}{cmd: qui reghdfe `0' nopartial}{p_end}
{p 12 12 2}{cmd: noi ereturn list}{p_end}
{p 12 12 2}{cmd: mata: st_local("n", strofreal(HDFE.num_singletons))}{p_end}
{p 12 12 2}{cmd: di as text "there are `n' singletons"}{p_end}
{phang2}{cmd: end}{p_end}
{phang2}{cmd: }{p_end}
{phang2}{cmd: qui include "reghdfe.mata", adopath}{p_end}
{phang2}{hline 54}{p_end}


{pstd}
{bf:B) Compute HDFE mata object, partial out the variables, but stop before regressing}

{p 8 15 2} {cmd:reghdfe} {it:...} {cmd:,}
{opt noreg:ress} [{help reghdfe##options_table:options}]{p_end}

{phang2}
This step is as A), but will also partial out the variables wrt. the fixed effects
and save the resulting information in the {it:HDFE.solution} object.
For instance, {it:HDFE.solution.data} will contained the partialled-out data,
and {it:HDFE.solution.depvar} will contain the name of the dependent variable.

{phang2}
This option can be used to (amongst other things) partial out all the variables only once, and then run regressions on the same sample and same regressors but with multiple left-hand-side variables (useful with very large datasets).


{pstd}
{bf:C) Run regression but keep the HDFE Mata object}

{p 8 15 2} {cmd:reghdfe} {it:...} {cmd:,}
{opt keepmata} [{help reghdfe##options_table:options}]{p_end}

{phang2}
By saving the HDFE object, this allows further manipulations of the fixed effects data,
although the data corresponding to the partialled-out variables is not preserved.


{marker mata}{...}
{title:3. FixedEffects Mata class}


{pstd}
In order to use reghdfe's Mata functions without your own ado-file, you need to add the following at the end of your file:

{phang2}{cmd: include "reghdfe.mata", adopath}{p_end}

{pstd}
This dynamically loads all the reghdfe Mata functions and classes, so they are accessible to the ado-file.
This alternative is preferred to sharing precompiled Mata objects, which would require compilation for multiple versions of Stata/Mata (or for the lowest possible version of Stata/Mata).

{pstd}
To construct the object, you can do:

{phang2}{cmd:class FixedEffects HDFE // Optional declaration}{p_end}
{phang2}{cmd:HDFE = FixedEffects() // Note that you can replace "HDFE" with whatever name you choose}{p_end}
{phang2}{cmd:HDFE.absvars = "firm_id year"}{p_end}
{phang2}{cmd:...}{p_end}
{phang2}{cmd:HDFE.init()}{p_end}
{phang2}{cmd:...}{p_end}

{pstd}
For more information, see the code of the {it:Estimate} function of {it:reghdfe.ado}


{marker properties}{...}
{title:Properties and Methods}

{pstd}{it:TODO: update this list}

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


{marker example}{...}
{title:Example: OLS regression}

{pstd}
TODO: Update this example

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
