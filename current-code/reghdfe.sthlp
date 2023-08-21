{smcl}
{* *! version 6.12.3 08aug2023}{...}
{vieweralsosee "[R] areg" "help areg"}{...}
{vieweralsosee "[R] xtreg" "help xtreg"}{...}
{vieweralsosee "" "--"}{...}
{vieweralsosee "reghdfe_programming" "help reghdfe_programming"}{...}
{vieweralsosee "ftools" "help ftools"}{...}
{vieweralsosee "ivreghdfe" "help ivreghdfe"}{...}
{vieweralsosee "ppmlhdfe" "help ppmlhdfe"}{...}
{vieweralsosee "sumhdfe" "help sumhdfe"}{...}
{vieweralsosee "ivreg2" "help ivreg2"}{...}
{vieweralsosee "did_imputation" "help did_imputation"}{...}
{vieweralsosee "parallel" "help parallel"}{...}
{vieweralsosee "" "--"}{...}
{viewerjumpto "Syntax" "reghdfe##syntax"}{...}
{viewerjumpto "Absorb syntax" "reghdfe##absorb"}{...}
{viewerjumpto "Description" "reghdfe##description"}{...}
{viewerjumpto "Group FEs" "reghdfe##group_fes"}{...}
{viewerjumpto "Options" "reghdfe##options"}{...}
{viewerjumpto "Postestimation Syntax" "reghdfe##postestimation"}{...}
{viewerjumpto "Examples" "reghdfe##examples"}{...}
{viewerjumpto "Stored results" "reghdfe##results"}{...}
{viewerjumpto "Authors" "reghdfe##contact"}{...}
{viewerjumpto "Acknowledgements" "reghdfe##acknowledgements"}{...}
{viewerjumpto "References" "reghdfe##references"}{...}
{title:Title}

{p2colset 5 18 20 2}{...}
{p2col :{cmd:reghdfe} {hline 2}}Linear regression with multiple fixed effects. Also supports individual FEs with group-level outcomes{p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{pstd}
{bf:Least-square regressions (no fixed effects):}

{p 8 15 2} {cmd:reghdfe}
{depvar} [{indepvars}] 
{ifin} {it:{weight}}
[{cmd:,} {help reghdfe##options_table:options}]
{p_end}


{pstd}
{bf:Fixed effects regressions:}

{p 8 15 2} {cmd:reghdfe}
{depvar} [{indepvars}] 
{ifin} {it:{weight}} {cmd:,} {opth a:bsorb(reghdfe##absorb:absvars)} [{help reghdfe##options_table:options}]{p_end}


{pstd}
{bf:Fixed effects regressions with group-level outcomes and individual FEs:}

{p 8 15 2} {cmd:reghdfe}
{depvar} [{indepvars}] 
{ifin} {it:{weight}} {cmd:,} {opt a:bsorb(absvars indvar)} {opt g:roup(groupvar)} {opt i:ndividual(indvar)} [{help reghdfe##options_table:options}]{p_end}


{marker options_table}{...}
{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Standard FEs {help reghdfe##opt_absorb:[+]}}
{synopt: {opth a:bsorb(reghdfe##absorb:absvars)}}categorical variables representing the fixed effects to be absorbed{p_end}
{synopt: {cmdab:a:bsorb(}{it:...}{cmd:,} {cmdab:save:fe)}}save all fixed effect estimates with the {it:__hdfe*} prefix{p_end}

{syntab:Group FEs {help reghdfe##opt_group_fes:[+]}}
{synopt : {opth g:roup(reghdfe##opt_groupvar:groupvar)}}categorical variable representing each group (eg: {it:patent_id}){p_end}
{synopt : }{bf:- note:} regression variables (depvar, indepvars) must be constant within each group (eg: {it:patent_citations} must be constant within a {it:patent_id}){p_end}
{synopt : }{bf:- note:} using {cmd:group()} without {cmd:individual()} is equivalent to running the regression on 1 observation per group{p_end}
{synopt : {opth i:ndividual(reghdfe##opt_indvar:indvar)}}categorical variable representing each individual whose fixed effect will be absorbed(eg: {it:inventor_id}){p_end}
{synopt : }{bf:- note:} the {cmd:individual()} option requires the {cmd:group()}{p_end}
{synopt : {opth ag:gregation(help reghdfe##opt_aggregation:str)}}how are the individual FEs aggregated within a group. Valid values are {it:mean} (default) and {it:sum}{p_end}
{synopt : }{bf:- note:} {it:mean} and {it:sum} are equivalent if all groups are of equal size (eg: 11 starting players in a football/soccer team){p_end}

{syntab:Model {help reghdfe##opt_model:[+]}}
{synopt : {opth vce:(reghdfe##model_opts:vcetype)}}{it:vcetype} may be {opt un:adjusted} (default), {opt r:obust} or {opt cl:uster} {help fvvarlist} (allowing two- and multi-way clustering){p_end}
{synopt : {opth res:iduals(newvar)}}save regression residuals{p_end}
{synopt : }{bf:- note:} the postestimation command "{it:predict <varname>, d}" requires this option{p_end}

{syntab:Degrees-of-Freedom Adjustments {help reghdfe##opt_dof:[+]}}
{synopt :{opt dof:adjustments(list)}}allows selecting the desired adjustments for degrees of freedom;
rarely used but changing it can speed-up execution{p_end}
{synopt: {opth groupv:ar(newvar)}}unique identifier for the first mobility group{p_end}

{syntab:Optimization {help reghdfe##opt_optimization:[+]}}
{synopt : {cmdab:tech:nique(map)}} partial out variables using the "method of alternating projections" (MAP) in any of its variants (default){p_end}
{synopt : {cmdab:tech:nique(lsmr)}} Fong and Saunders' LSMR algorithm{p_end}
{synopt : {cmdab:tech:nique(lsqr)}} Page and Saunders' LSQR algorithm{p_end}
{synopt : {cmdab:tech:nique(gt)}} Variation of Spielman et al's graph-theoretical (GT) approach (using spectral sparsification of graphs); currently disabled{p_end}

{synopt :{opt accel:eration(str)}}MAP acceleration method; options are conjugate_gradient ({opt cg}, default), steep_descent ({opt sd}), and {opt a:itken}{p_end}
{synopt :{opt transf:orm(str)}}MAP transform operation; options are {opt kac:zmarz}, {opt cim:mino}, and {opt sym:metric kaczmarz} (default){p_end}
{synopt :{opt prec:onditioner(str)}}LSMR/LSQR preconditioner. options are {opt no:ne}, {opt diag:onal}, and {opt block:_diagonal} (default){p_end}
{synopt :{opt prune}}prune vertices of degree-1; acts as a preconditioner that is useful if the underlying network is very sparse; currently disabled{p_end}

{synopt :{opt tol:erance(#)}}criterion for convergence (default=1e-8, valid values are 1e-1 to 1e-15){p_end}
{synopt :{opt iter:ate(#)}}maximum number of iterations (default=16,000); if set to missing ({cmd:.}) it will run for as long as it takes.{p_end}
{synopt :{opt nosamp:le}}will not create {it:e(sample)}, saving some space and speed{p_end}
{synopt :{opt fastreg:ress}}solve normal equations (X'X b = X'y) instead of the original problem (X=y). Faster but less accurate and less numerically stable. Use carefully{p_end}
{synopt :{opt keepsin:gletons}}do not drop singletons. Use carefully{p_end}

{syntab:Parallel execution {help reghdfe##opt_parallel:[+]}}
{synopt :{opth par:allel(reghdfe##opt_parallel:#)}}partial out variables in {it:#} separate Stata processes, speeding up execution depending on data size and computer characteristics. Requires the {help parallel} package{p_end}
{synopt :{opth par:allel(reghdfe##opt_parallel:#,cores(#2))}}specify that each process will only use #2 cores. More suboptions avalable {help parallel_map:here}{p_end}

{syntab:Memory Usage {help reghdfe##opt_memory:[+]}}
{synopt :{opth pool:size(#)}}apply the within algorithm in groups of {it:#} variables (else, it will run on all variables at the same time). A large pool size is usually faster but uses more memory{p_end}
{synopt :{opt compact}}preserve the dataset and drop variables as much as possible on every step{p_end}

{syntab:Reporting {help reghdfe##opt_reporting:[+]}}
{synopt :{opt l:evel(#)}}set confidence level; default is {cmd:level(95)}{p_end}
{synopt :{it:{help reghdfe##display_options:display_options}}}control columns and column formats, row spacing, line width,
display of omitted variables and base and empty cells, and factor-variable labeling{p_end}
{synopt :}particularly useful are the {opt noomit:ted} and {opt noempty} options to hide regressors omitted due to collinearity{p_end}
{synopt :{opt nohead:er}}suppress output header{p_end}
{synopt :{opt notable}}suppress coefficient table{p_end}
{synopt :{opt nofoot:note}}suppress fixed effects footnote{p_end}
{synopt :{opt nocon:stant}}suppress showing {it:_cons} row{p_end}

{syntab:Diagnostics {help reghdfe##opt_diagnostics:[+]}}
{synopt :{opt v:erbose(#)}}amount of debugging information to show (0=None, 1=Some, 2=More, 3=Parsing/convergence details, 4=Every iteration){p_end}
{synopt :{opt timeit}}show elapsed times by stage of computation{p_end}
{synopt :{opt version(#)}}run previous versions of reghdfe. Valid values are {help reghdfe3:3} (reghdfe 3, circa 2017) and {help reghdfe5:5} (reghdfe 5, circa 2020){p_end}
{synoptline}
{p 4 6 2}{it:depvar} and {it:indepvars} may contain {help tsvarlist:factor variables} and {help tsvarlist:time-series operators}.
{it:depvar} cannot be of the form {it:i.y} though, only {it:#.y} (where # is a number){p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:reghdfe} is a generalization of {help areg} (and {help xtreg:xtreg,fe}, {help xtivreg:xtivreg,fe}) for multiple levels of fixed effects, and multi-way clustering.

{pstd}
For alternative estimators (2sls, gmm2s, liml), as well as additional standard errors (HAC, etc) see {help ivreghdfe}. For nonlinear fixed effects, see {help ppmlhdfe} (Poisson).
For diagnostics on the fixed effects and additional postestimation tables, see {browse "https://github.com/ed-dehaan/sumhdfe":sumhdfe}.

{pstd}Additional features include:{p_end}

{phang2}  1. A novel and robust algorithm to efficiently absorb the fixed effects (extending the work of Guimaraes and Portugal, 2010).{p_end}
{phang2}  2. Can absorb heterogeneous slopes (i.e. regressors with different coefficients for each FE category){p_end}
{phang2}  3. Can absorb {browse "http://scorreia.com/research/individual_fes.pdf":individual fixed effects} where outcomes and regressors are at the group level
(e.g. controlling for inventor fixed effects using patent data where outcomes are at the patent level){p_end}
{phang2}  4. Can save fixed effect point estimates ({it:caveat emptor}: the fixed effects may not be identified, see the {help reghdfe##references:references}).{p_end}
{phang2}  5. Calculates the degrees-of-freedom lost due to the fixed effects
(note: beyond two levels of fixed effects, this is still an open problem, but we provide a conservative approximation).{p_end}
{phang2}  6. Iteratively removes singleton observations, to avoid biasing the standard errors (see ancillary document).{p_end}
{phang2}  7. Coded in Mata, which in most scenarios makes it even faster than {it:areg} and {it:xtreg} for a single fixed effect (see benchmarks on the Github page).{p_end}

{pstd}
For a description of its internal Mata API, as well as options for programmers, see the help file {help reghdfe_programming}.


{marker group_fes}{...}
{title:Description of individual fixed effects in group setting}

{pstd}{cmd:reghdfe} now permits estimations that include individual fixed effects with group-level outcomes.
For instance, a study of innovation might want to estimate patent citations as a function of patent characteristics, standard fixed effects (e.g. year),
and fixed effects for each inventor that worked in a patent.

{pstd}To do so, the data must be stored in a long format (e.g. with each patent spanning as many observations as inventors in the patent.)
Specifically, the individual and group identifiers must uniquely identify the observations
(so for instance the command "isid patent_id inventor_id" will not raise an error).
Note that this allows for groups with a varying number of individuals (e.g. one patent might be solo-authored, another might have 10 authors).

{pstd}Other example cases that highlight the utility of this include:

{pmore}  1. Patents & inventors{p_end}
{pmore}  2. Papers & co-authors{p_end}
{pmore}  3. Time-varying executive boards & board members{p_end}
{pmore}  4. Sports teams & players{p_end}


{pstd}For a more detailed explanation, including examples and technical descriptions, see {browse "http://scorreia.com/research/individual_fes.pdf":Constantine and Correia (2021)}.

{marker links}{...}
{title:Links to online documentation}

{p2col 8 10 10 2: -}{browse "http://scorreia.com/software/reghdfe/":Website}: main reghdfe website (including online help, quickstart, FAQ).{p_end}
{p2col 8 10 10 2: -}{browse "http://scorreia.com/software/reghdfe/":Github page}: code repository, issues/problems/suggestions, and latest news.{p_end}
{p2col 8 10 10 2: -}{browse "http://scorreia.com/research/hdfe.pdf":HDFE paper}: explain the algorithms behind reghdfe.{p_end}
{p2col 8 10 10 2: -}{browse "http://scorreia.com/research/individual_fes.pdf":Individual fixed effects paper}: explain the algorithms behind individual fixed effects in reghdfe.{p_end}
{p2col 8 10 10 2: -}{browse "http://scorreia.com/research/groupfe.pdf":Group FE paper}: illustrate the importance of using individual fixed effects with group-level outcomes.{p_end}


{marker absorb}{...}
{title:Absorb() syntax}

{synoptset 22}{...}
{synopthdr:absvar}
{synoptline}
{synopt:{it:varname}}categorical variable to be absorbed{p_end}
{synopt:{cmd:i.}{it:varname}}categorical variable to be absorbed (same as above; the {cmd:i.} prefix is always implicit){p_end}
{synopt:{cmd:i.}{it:var1}{cmd:#i.}{it:var2}}absorb the interactions of multiple categorical variables{p_end}
{synopt:{cmd:i.}{it:var1}{cmd:#}{cmd:c.}{it:var2}}absorb heterogeneous slopes, where {it:var2} has a different slope estimate depending on {it:var1}. Use carefully (see below!){p_end}
{synopt:{it:var1}{cmd:##}{cmd:c.}{it:var2}}absorb heterogenous intercepts and slopes. Equivalent to "{cmd:i.}{it:var1} {cmd:i.}{it:var1}{cmd:#}{cmd:c.}{it:var2}", but {it:much} faster{p_end}
{synopt:{it:var1}{cmd:##c.(}{it:var2 var3}{cmd:)}}multiple heterogeneous slopes are allowed together. Alternative syntax: {it:var1}{cmd:##(c.}{it:var2} {cmd:c.}{it:var3}{cmd:)}{p_end}
{synopt:{it:v1}{cmd:#}{it:v2}{cmd:#}{it:v3}{cmd:##c.(}{it:v4 v5}{cmd:)}}factor operators can be combined{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}- To save the estimates of specific absvars, write {newvar}{inp:={it:absvar}}.{p_end}
{p 4 6 2}-  However, be aware that estimates for the fixed effects are generally inconsistent and not econometrically identified.{p_end}
{p 4 6 2}- Using categorical interactions (e.g. {it:x}{cmd:#}{it:z}) is easier and faster than running {it:egen group(...)} beforehand.{p_end}
{p 4 6 2}- {browse "http://scorreia.com/research/singletons.pdf":Singleton observations} are dropped iteratively until no more singletons are found (see the linked article for details).{p_end}
{p 4 6 2}- Slope-only absvars ("state#c.time") have poor numerical stability and slow convergence.
If you need those, either i) increase tolerance or
ii) use slope-and-intercept absvars ("state##c.time"), even if the intercept is redundant.
For instance if absvar is "i.zipcode i.state##c.time" then i.state is redundant given i.zipcode, but
convergence will still be {it:much} faster.{p_end}


{marker options}{...}
{title:Options}

{marker opt_absorb}{...}
{dlgtab:Standard FEs}

{phang}
{opth a:bsorb(reghdfe##absorb:absvars)} list of categorical variables (or interactions) representing the fixed effects to be absorbed.
This is equivalent to including an indicator/dummy variable for each category of each {it:absvar}. {cmd:absorb()} is required.

{pmore}
To save a fixed effect, prefix the absvar with "{newvar}{cmd:=}".
For instance, the option {cmd:absorb(firm_id worker_id year_coefs=year_id)} will include firm, worker, and year fixed effects,
but will only save the estimates for the year fixed effects (in the new variable {it:year_coefs}).

{pmore}
If you want to run {help reghdfe##postestimation:predict} afterward but don't particularly care about the names of each fixed effect, use the {cmdab:save:fe} suboption.
This will delete all preexisting variables matching {it:__hdfe*__} and create new ones as required.
Example: {it:reghdfe price weight, absorb(turn trunk, savefe)}.

{marker opt_group_fes}{...}
{dlgtab:Group FEs}

{marker opt_groupvar}{...}
{phang}
{opth g:roup(reghdfe##opt_groupvar:groupvar)} categorical variable representing each group (eg: {it:patent_id}). {cmd:group()} is not required, {it:unless} you specify {cmd:individual()}.

{pmore}If {it:only} {cmd:group()} is specified, the program will run with one observation per group. 

{pmore}Note that group here means whatever aggregation unit at which the outcome is defined.

{marker opt_indvar}{...}

{phang}
{opth i:ndividual(reghdfe##opt_indvar:indvar)} categorical variable representing each individual (eg: {it:inventor_id}).

{pmore}This variable is not automatically added to {cmd:absorb()}, so {it:you must include} it in the absvar list.
This is because the order in which you include it affects the speed of the command, and {cmd:reghdfe} is not smart enough to know the optimal ordering.

{pmore}If {cmd:individual()} is specified you must also call {cmd:group()}.

{marker opt_aggregation}{...}
{phang}
{opt ag:gregation(str)} method of aggregation for the individual components of the group fixed effects.
Valid options are {cmd:mean} (default), and {cmd:sum}.

{pmore}If all groups are of equal size, both options are equivalent and result in identical estimates.

{pmore}Note that both options are econometrically valid, and {cmd:aggregation()} should be determined based on the economics behind each specification.
For instance, adding more authors to a paper or more inventors to an invention might not increase its quality proportionally (i.e. its citations), so using "mean" might be the sensible choice.
In contrast, other production functions might scale linearly in which case "sum" might be the correct choice.

{phang} 
{bf:Combining options:} depending on which of {cmd:absorb()}, {cmd:group()}, and {cmd:individual()} you specify, you will trigger different use cases of reghdfe:

{pmore}  1. If none is specified, reghdfe will run OLS with a constant.{p_end}
{pmore}  2. If only {cmd:absorb()} is present, reghdfe will run a standard fixed-effects regression.{p_end}
{pmore}  3. If {cmd:group()} is specified (but not {cmd:individual()}), this is equivalent to #1 or #2 with only {it:one} observation per group.
That is, running {it:"bysort group: keep if _n == 1"} and then {it:"reghdfe ..."}.{p_end}
{pmore}  3. If all are specified, this is equivalent to a fixed-effects regression at the group level and individual FEs.{p_end}

{marker opt_model}{...}
{dlgtab:Model}

{marker opt_vce}{...}
{phang}
{opth vce:(reghdfe##vcetype:vcetype, subopt)} specifies the type of standard error reported.

{pmore}
{opt un:adjusted}|{opt ols:} estimates conventional standard errors, valid under the assumptions of homoscedasticity and no correlation between observations even in small samples.

{pmore}
{opt r:obust} estimates heteroscedasticity-consistent standard errors (Huber/White/sandwich estimators), which still assume independence between observations.

{pmore}Warning: in a FE panel regression, using {opt r:obust} will
lead to inconsistent standard errors if, for every fixed effect, the {it:other} dimension is fixed.
For instance, in a standard panel with individual and time fixed effects, we require both the number of
individuals and periods to grow asymptotically.
If that is not the case, an alternative may be to use clustered errors,
which as discussed below will still have their own asymptotic requirements.
For a discussion, see
{browse "http://www.princeton.edu/~mwatson/papers/ecta6489.pdf":Stock and Watson, "Heteroskedasticity-robust standard errors for fixed-effects panel-data regression," Econometrica 76 (2008): 155-174}.

{pmore}
{opt cl:uster} {it:clustervars} estimates consistent standard errors even when the observations
are correlated within groups.

{pmore}
Multi-way-clustering is allowed. Thus, you can indicate as many {it:clustervar}s as desired
(e.g. allowing for intragroup correlation across individuals, time, country, etc).
For instance, {it:vce(cluster firm year)} will estimate SEs with firm and year clustering (two-way clustering).

{pmore}
Each {it:clustervar} permits interactions of the type {it:var1{cmd:#}var2}.
This is equivalent to using {it:egen group(var1 var2)} to create a new variable, but more convenient and faster.
For instance, {it:vce(cluster firm#year)} will estimate SEs with {bf:one-way} clustering i.e. where all observations of a given firm and year are clustered together. 

{pmore} {bf:Note:} do not confuse {it:vce(cluster firm#year)} (one-way clustering) with {it:vce(cluster firm year)} (two-way clustering).

{pmore} {bf:Warning:} it is not recommended to run clustered SEs if any of the clustering variables have too few different levels.
A frequent rule of thumb is that each cluster variable must have at least 50 different categories
(the number of categories for each clustervar appears at the top of the regression table).

{pmore} {bf:Note:} More advanced SEs, including autocorrelation-consistent (AC), heteroskedastic and
autocorrelation-consistent (HAC), Driscoll-Kraay, Kiefer, etc. are available in the {help ivreghdfe} package (which uses {help ivreg2} as its back-end).

{phang}
{opth res:iduals(newvar)} saves the regression residuals in a new variable. 

{pmore} {opt res:iduals} (without parenthesis) saves the residuals
in the variable {it:_reghdfe_resid} (overwriting it if it already exists).

{pmore}
This option does not require additional computations and is required for
subsequent calls to {cmd:predict, d}.

{phang}
{opth su:mmarize(tabstat##statname:stats)} this option is now part of {help sumhdfe}

{marker opt_iv}{...}
{dlgtab:IV/2SLS/GMM}

{phang}
The IV functionality of {cmd:reghdfe} has been moved into {help ivreghdfe}.

{marker opt_dof}{...}
{dlgtab:Degrees-of-Freedom Adjustments}

{phang}
{opt dof:adjustments(doflist)} selects how the degrees-of-freedom, as well as e(df_a), are adjusted due to the absorbed fixed effects.

{pmore}
{bf: The problem:} without any adjustment, the degrees-of-freedom (DoF) lost due to the fixed effects is equal to the count of all the fixed effects.
For instance, a regression with {it:absorb(firm_id worker_id)}, and 1000 firms, 1000 workers, would drop 2000 DoF due to the FEs.
This is potentially too aggressive, as many of these fixed effects might be perfectly collinear with each other, and the true number of DoF lost might be lower.
As a consequence, your standard errors might be erroneously too large.

{pmore}
{bf: The solution:} To address this, reghdfe uses several methods to count instances as possible of collinearities of FEs.
In most cases, it will count all instances (e.g. one- and two-way fixed effects), but in others it will only provide a conservative estimate.
Doing this is relatively slow, so reghdfe might be sped up by changing these options.

{pmore}
{opt all} is the default and usually the best alternative. It is equivalent to {opt dof(pairwise clusters continuous)}.
However, an alternative when using many FEs is to run {opt dof(firstpair clusters continuous)}, which is faster and might be almost as good.

{pmore}
{opt none} assumes no collinearity across the fixed effects (i.e. no redundant fixed effects). This is overtly conservative, although it is the faster method by virtue of not doing anything.

{pmore}
{opt first:pair} will exactly identify the number of collinear fixed effects across the first two sets of fixed effects
(i.e. the first absvar and the second absvar).
The algorithm used for this is described in Abowd et al (1999), and relies on results from graph theory
(finding the number of connected sub-graphs in a bipartite graph).
It will not do anything for the third and subsequent sets of fixed effects.

{pmore}
For more than two sets of fixed effects, there are no known results that provide exact degrees-of-freedom as in the case above.
One solution is to ignore subsequent fixed effects (and thus overestimate e(df_a) and underestimate the degrees-of-freedom).
Another solution, described below, applies the algorithm between pairs of fixed effects to obtain a better (but not exact) estimate:

{pmore}
{opt pair:wise} applies the aforementioned connected-subgraphs algorithm between pairs of fixed effects.
For instance, if there are four sets of FEs, the first dimension will usually have no redundant coefficients (i.e. e(M1)==1), since we are running the model without a constant.
For the second FE, the number of connected subgraphs with respect to the first FE will provide an exact estimate of the degrees-of-freedom lost, e(M2).

{pmore}
For the third FE, we do not know exactly.
However, we can compute the number of connected subgraphs between the first and third {it:G(1,3)},
and second and third {it:G(2,3)} fixed effects,
and choose the higher of those as the closest estimate for e(M3).
For the fourth FE, we compute {it:G(1,4)}, {it:G(2,4)}, and {it:G(3,4)} and again choose the highest for e(M4).

{pmore}
Finally, we compute e(df_a) = e(K1) - e(M1) + e(K2) - e(M2) + e(K3) - e(M3) + e(K4) - e(M4);
where e(K#) is the number of levels or dimensions for the #-th fixed effect (e.g. number of individuals or years).
Note that e(M3) and e(M4) are only conservative estimates and thus we will usually be overestimating the standard errors.
However, given the sizes of the datasets typically used with reghdfe, the difference should be small. 

{pmore}
Since the gain from {opt pair:wise} is usually {it:minuscule} for large datasets, and the computation is expensive, it may be a good practice to exclude this option for speedups.

{pmore}
{opt cont:inuous}
Fixed effects with continuous interactions (i.e. individual slopes, instead of individual intercepts) are dealt with differently.
In an i.categorical#c.continuous interaction, we will do one check: we count the number of categories where c.continuous is always zero.
In an i.categorical##c.continuous interaction, we count the number of categories where c.continuos is always the same constant. If that is the case, then the slope is collinear with the intercept.

{pmore}
Additional methods, such as {opt bootstrap} are also possible but not yet implemented.
Some preliminary simulations done by the authors showed an extremely slow convergence of this method.

{phang}
{opth groupv:ar(newvar)} name of the new variable that will contain the first mobility group.
Requires {opt pair:wise}, {opt first:pair}, or the default {opt all}.

{marker opt_optimization}{...}
{dlgtab:Optimization}

{phang}
{opt tech:nique(str)} 

{pmore}{cmd:technique(map)} (default)will partial out variables using the "method of alternating projections" (MAP) in any of its variants. MAP currently does not work with individual & group fixed effects. Fast and stable option

{pmore}{cmd:technique(lsmr)} use the Fong and Saunders LSMR algorithm. Recommended (default) technique when working with individual fixed effects.
LSMR is an iterative method for solving sparse least-squares problems; analytically equivalent to the MINRES method on the normal equations.
For more information on the algorithm, please reference the {browse "https://stanford.edu/group/SOL/software/lsmr/LSMR-SISC-2011.pdf":paper}

{pmore}{cmd:technique(lsqr)} use Paige and Saunders LSQR algorithm. Alternative technique when working with individual fixed effects.
LSQR is an iterative method for solving sparse least-squares problems; analytically equivalent to conjugate gradient method on the normal equations.
Fast, but less precise than LSMR at default tolerance (1e-8).
For more information on the algorithm, please reference the {browse "https://web.stanford.edu/group/SOL/software/lsqr/lsqr-toms82a.pdf":paper}

{pmore}{cmd:technique(gt)} variation of Spielman et al's graph-theoretical (GT) approach (using a spectral sparsification of graphs); currently disabled

{phang}
{opt acceleration(str)} Relevant for {cmd:tech(map)}.
Allows for different acceleration techniques, from the simplest case of no acceleration ({opt no:ne}), to steep descent ({opt st:eep_descent} or {opt sd}), Aitken ({opt a:itken}),
and finally Conjugate Gradient ({opt co:njugate_gradient} or {opt cg}).

{pmore}
Note: Each acceleration is just a plug-in Mata function, so a larger number of acceleration techniques are available, albeit undocumented (and slower).

{phang}
{opt transf:orm(str)} allows for different "alternating projection" transforms.
The classical transform is Kaczmarz ({opt kac:zmarz}), and more stable alternatives are Cimmino ({opt cim:mino}) and Symmetric Kaczmarz ({opt sym:metric_kaczmarz})

{pmore}
Note: The default acceleration is Conjugate Gradient and the default transform is Symmetric Kaczmarz.
Be wary that different accelerations often work better with certain transforms.
For instance, {bf:do not} use conjugate gradient with plain Kaczmarz, as it will not converge (this is because CG requires a symmetric operator in order to converge,
and plain Kaczmarz is not symmetric).

{phang}
{opt prec:onditioner(str)} LSMR/LSQR require a good preconditioner in order to converge efficiently and in few iterations.
reghfe currently supports right-preconditioners of the following types: {opt no:ne}, {opt diag:onal}, and {opt block:_diagonal} (default).

{phang}
{opt prune(str)}prune vertices of degree-1; acts as a preconditioner that is useful if the underlying network is very sparse; currently disabled

{phang}
{opth tol:erance(#)} specifies the tolerance criterion for convergence; default is {cmd:tolerance(1e-8)}.
In general, high tolerances (1e-8 to 1e-14) return more accurate results, but more slowly.
Similarly, low tolerances (1e-7, 1e-6, ...) return faster but potentially inaccurate results.

{pmore}
Note that tolerances higher than 1e-14 might be problematic, not just due to speed, but because they approach the limit of the computer precision (1e-16).
Thus, using e.g. {it:tol(1e15)} might not converge, or take an inordinate amount of time to do so.

{pmore}
At the other end, low tolerances (below 1e-6) are not generally recommended, as the iteration might have been stopped too soon, and thus the reported estimates might be incorrect.
However, with very large datasets, it is sometimes useful to use low tolerances when running preliminary estimates.

{pmore}
Note: detecting perfectly collinear regressors is more difficult with iterative methods (i.e. those used by reghdfe) than with direct methods (i.e. those used by regress).
To spot perfectly collinear regressors that were not dropped, look for extremely high standard errors. In this case, consider using higher tolerances.

{pmore}
Warning: when absorbing heterogeneous slopes without the accompanying heterogeneous intercepts,
convergence is quite poor and a higher tolerance is strongly suggested (i.e. higher than the default).
In other words, an absvar of {it:var1##c.var2} converges easily, but an absvar of {it:var1#c.var2} will converge slowly and may require a higher tolerance.

{phang}
{opth it:erations(#)}
specifies the maximum number of iterations; the default is {cmd:iterations(16000)}; set it to missing ({cmd:.}) to run forever until convergence.

{phang}
{opt nosample} will not create {it:e(sample)}, saving some space and speed.

{marker opt_parallel}{...}
{dlgtab:Parallel execution}

{phang}
{opt par:allel(#1, cores(#2) options)} runs the partialling-out step in #1 separate Stata processeses,
each using #2 cores.
This option requires the {help parallel} package (see {browse "https://github.com/gvegayon/parallel":website}). There are several additional suboptions, discussed {help parallel_map:here}.

{pmore}
Note that {cmd:parallel()} will only speed up execution in certain cases.
First, the dataset needs to be large enough, and/or the partialling-out process needs to be slow enough, that the overhead of opening separate Stata instances will be worth it.
Second, if the computer has only one or a few cores, or limited memory, it might not be able to achieve significant speedups.

{marker opt_memory}{...}
{dlgtab:Memory Usage}
{phang}
{opth pool:size(#)}
Number of variables that are {it:pooled together} into a matrix that will then be transformed.
The default is to pool variables in groups of 10. Larger groups are faster with more than one processor, but may cause out-of-memory errors. In that case, set poolsize to 1.

{phang} 
{opt compact} preserve the dataset and drop variables as much as possible on every step

{marker opt_reporting}{...}
{dlgtab:Reporting}

{phang}
{opt l:evel(#)} sets confidence level; default is {cmd:level(95)}; see {helpb estimation options##level():[R] Estimation options}

INCLUDE help displayopts_list

{phang}
{opt nohead:er} suppresses the display of the table of summary
statistics at the top of the output; only the coefficient table is displayed.
This option is often used in programs and ado-files.

{phang}
{opt notable} suppresses display of the coefficient table.

{phang}
{opt nofoot:note} suppresses display of the footnote table that lists the absorbed fixed effects, including the number of categories/levels of each fixed effect,
redundant categories (collinear or otherwise not counted when computing degrees-of-freedom), and the difference between both.

{phang}
{opt nocon:stant} suppresses display of the {it:_cons} row in the main table. No results or computations change, this is merely a cosmetic option

{marker opt_diagnostics}{...}
{dlgtab:Diagnostic}

{phang}
{opt v:erbose(#)} orders the command to print debugging information.

{pmore}
Possible values are 0 (none), 1 (some information), 2 (even more), 3 (adds dots for each iteration, and reports parsing details), 4 (adds details for every iteration step)

{pmore}
For debugging, the most useful value is 3. For simple status reports, set verbose to 1.

{phang}
{opt time:it} shows the elapsed time at different steps of the estimation. Most time is usually spent on three steps: map_precompute(), map_solve() and the regression step.

{phang}
{opt version(#)} reghdfe has had so far two large rewrites, from version 3 to 4, and version 5 to version 6.
Because the rewrites might have removed certain features (e.g. IV/2SLS was available in version 3 but moved to ivreghdfe on version 4),
this option allows you to run the previous versions without having to install them (they are already included in reghdfe installation).

{pmore}
To use them, just add the options {it:version(3)} or {it:version(5)}. You can check their respective help files here: {help reghdfe3}, {help reghdfe5}.

{pmore}
This option is also useful when replicating older papers, or to verify the correctness of estimates under the latest version.

{pmore}
Tip:To avoid the warning text in red, you can add the undocumented {cmd:nowarn} option.


{marker postestimation}{...}
{title:Postestimation Syntax}


{pstd}
Only {cmd:estat summarize}, {cmd:predict}, and {cmd:test} are currently supported and tested.

{pstd}
For additional postestimation tables specifically tailored to fixed effect models, see the {browse "https://github.com/ed-dehaan/sumhdfe":sumhdfe} package.

{pstd}
The syntax of {it: estat summarize} and {it:predict} is:

{p 8 13 2}
{cmd:estat summarize}
{p_end}{col 23}Summarizes {it:depvar} and the variables described in {it:_b} (i.e. not the excluded instruments)

{p 8 16 2}
{cmd:predict} 
{newvar} 
{ifin}
[{cmd:,} {it:statistic}]
{p_end}{col 23}May require you to previously save the fixed effects (except for option {opt xb}).
{col 23}To see how, see the details of the {help reghdfe##absvar:absorb} option
{col 23}Equation: y = xb + d_absorbvars + e

{synoptset 20 tabbed}{...}
{synopthdr:statistic}
{synoptline}
{syntab :Main}
{p2coldent: {opt xb}}xb fitted values; the default{p_end}
{p2coldent: {opt xbd}}xb + d_absorbvars{p_end}
{p2coldent: {opt d}}d_absorbvars{p_end}
{p2coldent: {opt r:esiduals}}residual{p_end}
{p2coldent: {opt sc:ore}}score; equivalent to {opt residuals}{p_end}
{p2coldent: {opt stdp}}standard error of the prediction (of the xb component){p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}although {cmd:predict} {help data_types:type} {help newvar} is allowed,
the resulting variable will always be of type {it:double}.{p_end}


{p 8 16 2}
{cmd:test} Performs significance test on the parameters, see the {help test:stata help}

{p 8 16 2}
{cmd:suest} Do not use {cmd:suest}. It will run, but the results will be incorrect. See workaround below

{pmore}If you want to perform tests that are usually run with {cmd:suest},
such as non-nested models, tests using alternative specifications of the variables,
or tests on different groups, you can replicate it manually, as described 
{browse "http://www.stata.com/statalist/archive/2009-11/msg01485.html":here}.
{p_end}


{title:Missing Features}

{phang}(If you are interested in discussing these or others, feel free to {help reghdfe##contact:contact us})

{p2colset 8 12 12 2}{...}
{p2col: -}Implement a -bootstrap- option{p_end}
{p2col: -}Improve/reincorporate {opt tech(gt)} and {opt prune} options{p_end}
{p2col: -}Improve DoF adjustments for 3+ HDFEs (e.g. as discussed in the {browse "https://ideas.repec.org/c/boc/bocode/s458181.html":group3hdfe} package){p_end}
{p2col: -}More postestimation commands (lincom? margins?){p_end}
{p2colreset}{...}

{marker examples}{...}
{title:Examples}

{hline}
{pstd}Setup{p_end}
{phang2}{cmd:. sysuse auto}{p_end}

{pstd}Simple case - one fixed effect{p_end}
{phang2}{cmd:. reghdfe price weight length, absorb(rep78)}{p_end}
{hline}

{pstd}As above, but also compute clustered standard errors{p_end}
{phang2}{cmd:. reghdfe price weight length, absorb(rep78) vce(cluster rep78)}{p_end}
{hline}

{pstd}Two and three sets of fixed effects{p_end}
{phang2}{cmd:. webuse nlswork}{p_end}
{phang2}{cmd:. reghdfe ln_w grade age ttl_exp tenure not_smsa south , absorb(idcode year)}{p_end}
{phang2}{cmd:. reghdfe ln_w grade age ttl_exp tenure not_smsa south , absorb(idcode year occ)}{p_end}
{hline}

{title:Advanced examples}

{pstd}Save the FEs as variables{p_end}
{phang2}{cmd:. reghdfe ln_w grade age ttl_exp tenure not_smsa south , absorb(FE1=idcode FE2=year)}{p_end}

{pstd}Interactions in the absorbed variables (notice that only the {it:#} symbol is allowed){p_end}
{phang2}{cmd:. reghdfe ln_w grade age ttl_exp tenure not_smsa , absorb(idcode#occ)}{p_end}

{title:Group Examples}

{hline}
{pstd}Setup{p_end}
{phang2}{cmd:. webuse toy-patents-long}{p_end}

{pstd}Individual (inventor) & group (patent) fixed effects{p_end}
{phang2}{cmd:. reghdfe citations funding, a(inventor_id) group(patent_id) individual(inventor_id)}{p_end}
{hline}

{pstd}Individual & group fixed effects, with an additional standard fixed effects variable{p_end}
{phang2}{cmd:. reghdfe citations funding, a(year inventor_id) group(patent_id) individual(inventor_id)}{p_end}
{hline}

{pstd}Individual & group fixed effects, specifying with a different method of aggregation (sum){p_end}
{phang2}{cmd:. reghdfe citations funding, a(inventor_id) group(patent_id) individual(inventor_id) aggreg(sum)}{p_end}
{pstd}If theory suggests that the effect of multiple authors will enter additively, as opposed to the average effect of the group of authors, this would be the appropriate treatment.
Mean is the default method.{p_end}
{hline}

{pstd}Use one observation per group{p_end}
{phang2}{cmd:. reghdfe citations funding, a(year) group(patent_id)}{p_end}
{hline}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:reghdfe} stores the following in {cmd:e()}:

{synoptset 24 tabbed}{...}
{syntab:Scalars}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(num_singletons)}}number of singleton observations{p_end}
{synopt:{cmd:e(N_full)}}number of observations including singletons{p_end}

{synopt:{cmd:e(N_hdfe)}}number of absorbed fixed-effects{p_end}
{synopt:{cmd:e(tss)}}total sum of squares{p_end}
{synopt:{cmd:e(tss)}}total sum of squares after partialling-out{p_end}
{synopt:{cmd:e(rss)}}residual sum of squares{p_end}
{synopt:{cmd:e(rss)}}model sum of squares (tss-rss){p_end}
{synopt:{cmd:e(r2)}}R-squared{p_end}
{synopt:{cmd:e(r2_a)}}adjusted R-squared{p_end}
{synopt:{cmd:e(r2_within)}}Within R-squared{p_end}
{synopt:{cmd:e(r2_a_within)}}Adjusted Within R-squared{p_end}
{synopt:{cmd:e(df_a)}}degrees of freedom lost due to the fixed effects{p_end}
{synopt:{cmd:e(rmse)}}root mean squared error{p_end}
{synopt:{cmd:e(ll)}}log-likelihood{p_end}
{synopt:{cmd:e(ll_0)}}log-likelihood of fixed-effect-only regression{p_end}
{synopt:{cmd:e(F)}}F statistic{p_end}
{synopt:{cmd:e(rank)}}rank of {cmd:e(V)}{p_end}
{synopt:{cmd:e(N_clustervars)}}number of cluster variables{p_end}
        
{synopt:{cmd:e(N_clust}#{cmd:)}}number of clusters for the #th cluster variable{p_end}
{synopt:{cmd:e(N_clust)}}number of clusters; minimum of {it:e(clust#)}{p_end}

{synopt:{cmd:e(df_m)}}model degrees of freedom{p_end}
{synopt:{cmd:e(df_r)}}residual degrees of freedom{p_end}

{synopt:{cmd:e(sumweights)}}sum of weights{p_end}
{synopt:{cmd:e(ic)}}number of iterations{p_end}
{synopt:{cmd:e(converged)}}{cmd:1} if converged, {cmd:0} otherwise{p_end}
{synopt:{cmd:e(drop_singletons)}}{cmd:1} if singletons were dropped, {cmd:0} otherwise{p_end}
{synopt:{cmd:e(df_a_nested)}}Redundant due to being nested within clustervars{p_end}
{synopt:{cmd:e(report_constant)}}whether _cons was included in the regressions (default)
or as part of the fixed effects{p_end}

{synoptset 24 tabbed}{...}
{syntab:Macros}
{synopt:{cmd:e(cmd)}}{cmd:reghdfe}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(dofmethod)}}dofmethod employed in the regression{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(indepvars)}}names of independent variables{p_end}
{synopt:{cmd:e(absvars)}}name of the absorbed variables or interactions{p_end}
{synopt:{cmd:e(extended_absvars)}}name of the extended absorbed variables (counting intercepts and slopes separately){p_end}
{synopt:{cmd:e(clustvar)}}name of cluster variable{p_end}
{synopt:{cmd:e(clustvar}#{cmd:)}}name of the #th cluster variable{p_end}
{synopt:{cmd:e(vce)}}{it:vcetype} specified in {cmd:vce()}{p_end}
{synopt:{cmd:e(vcetype)}}title used to label Std. Err.{p_end}
{synopt:{cmd:e(properties)}}{cmd:b V}{p_end}
{synopt:{cmd:e(predict)}}program used to implement {cmd:predict}{p_end}
{synopt:{cmd:e(estat_cmd)}}program used to implement {cmd:estat}{p_end}
{synopt:{cmd:e(footnote)}}program used to display footnote{p_end}
{synopt:{cmd:e(dofmethod)}}method(s) used to compute degrees-of-freedom lost due the fixed effects{p_end}
{synopt:{cmd:e(marginsnotok)}}predictions not allowed by {cmd:margins}{p_end}
{synopt:{cmd:e(title)}}title in estimation output{p_end}
{synopt:{cmd:e(title2)}}subtitle in estimation output, indicating how many FEs were being absorbed{p_end}

{synoptset 24 tabbed}{...}
{syntab:Matrices}
{synopt:{cmd:e(b)}}coefficient vector{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix of the estimators{p_end}
{synopt:{cmd:e(dof_table)}}degrees-of-freedom table{p_end}
{synopt:{cmd:r(table)}}main results table{p_end}

{synoptset 24 tabbed}{...}
{syntab:Functions}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}
{p2colreset}{...}


{marker contact}{...}
{title:Authors}

{pstd}Sergio Correia{break}
Board of Governors of the Federal Reserve{break}
Email: {browse "mailto:sergio.correia@gmail.com":sergio.correia@gmail.com}
{p_end}

{pstd}Noah Constantine{break}
Board of Governors of the Federal Reserve{break}
Email: {browse "mailto:noahbconstantine@gmail.com":noahbconstantine@gmail.com}
{p_end}


{marker support}{...}
{title:Support and updates}

{pstd}{cmd:reghdfe} requires the {cmd:ftools} package
({browse "https://github.com/sergiocorreia/ftools":Github repo}).{p_end}


{marker acknowledgements}{...}
{title:Acknowledgements}

{pstd}
This package wouldn't have existed without the invaluable feedback and contributions of Paulo Guimarães, Amine Ouazad, Mark E. Schaffer, Kit Baum, Tom Zylkin, and Matthieu Gomez.
Also invaluable are the great bug-spotting abilities of many users.{p_end}

{pstd}In addition, {it:reghdfe} is built upon important contributions from the Stata community:{p_end}

{phang}{browse "https://ideas.repec.org/c/boc/bocode/s457101.html":reg2hdfe}, from Paulo Guimaraes,
and {browse "https://ideas.repec.org/c/boc/bocode/s456942.html":a2reg} from Amine Ouazad,
were the inspiration and building blocks on which reghdfe was built.{p_end}
 
{phang}{browse "http://www.repec.org/bocode/i/ivreg2.html":ivreg2}, by Christopher F Baum, Mark E Schaffer, and Steven Stillman,
is the package used by default for instrumental-variable regression.{p_end}

{phang}{browse "https://github.com/gvegayon/parallel":parallel} by George Vega Yon and Brian Quistorff, is for parallel processing.{p_end}

{phang}{browse "https://ideas.repec.org/c/boc/bocode/s457689.html":avar} by Christopher F Baum and Mark E Schaffer,
is the package used for estimating the HAC-robust standard errors of ols regressions.{p_end}

{phang}{browse "http://econpapers.repec.org/software/bocbocode/s456797.htm":tuples} by Joseph Lunchman and Nicholas Cox,
is used when computing standard errors with multi-way clustering (two or more clustering variables).{p_end}


{marker references}{...}
{title:References}

{p 0 0 2}
The algorithm underlying reghdfe is a generalization of the works by:

{phang}
Paulo Guimaraes and Pedro Portugal. "A Simple Feasible Alternative Procedure to Estimate
Models with High-Dimensional Fixed Effects".
{it:Stata Journal, 10(4), 628-649, 2010.}
{browse "http://www.stata-journal.com/article.html?article=st0212":[link]}
{p_end}

{phang}
Simen Gaure. "OLS with Multiple High Dimensional Category Dummies".
{it:Memorandum 14/2010, Oslo University, Department of Economics, 2010.}
{browse "https://ideas.repec.org/p/hhs/osloec/2010_014.html":[link]}
{p_end}

{p 0 0 2}
It addresses many of the limitations of previous works, such as possible lack of convergence, arbitrary slow convergence times,
and being limited to only two or three sets of fixed effects (for the first paper).
The paper explaining the specifics of the algorithm is a work-in-progress and available upon request.

{p 0 0 0}
If you use this program in your research, please cite either
the {browse "https://ideas.repec.org/c/boc/bocode/s457874.html":REPEC entry}
or the aforementioned papers.{p_end}


{title:Additional References}

{p 0 0 0}
For details on the Aitken acceleration technique employed, please see "method 3" as described by:

{phang}
Macleod, Allan J. "Acceleration of vector sequences by multi-dimensional Delta-2 methods."
{it:Communications in Applied Numerical Methods 2.4 (1986): 385-392.}
{p_end}
