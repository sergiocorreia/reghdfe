{smcl}
{* *! version 3.0.0  09may2015}{...}
{vieweralsosee "[R] areg" "help areg"}{...}
{vieweralsosee "[R] xtreg" "help xtreg"}{...}
{vieweralsosee "[R] ivregress" "help ivregress"}{...}
{vieweralsosee "" "--"}{...}
{vieweralsosee "ivreg2" "help ivreg2"}{...}
{vieweralsosee "ivregress" "help ivregress"}{...}
{vieweralsosee "reg2hdfe" "help reg2hdfe"}{...}
{vieweralsosee "a2reg" "help a2reg"}{...}
{viewerjumpto "Syntax" "reghdfe##syntax"}{...}
{viewerjumpto "Description" "reghdfe##description"}{...}
{viewerjumpto "Options" "reghdfe##options"}{...}
{viewerjumpto "Postestimation Syntax" "reghdfe##postestimation"}{...}
{viewerjumpto "Remarks" "reghdfe##remarks"}{...}
{viewerjumpto "Examples" "reghdfe##examples"}{...}
{viewerjumpto "Stored results" "reghdfe##results"}{...}
{viewerjumpto "Author" "reghdfe##contact"}{...}
{viewerjumpto "Updates" "reghdfe##updates"}{...}
{viewerjumpto "Acknowledgements" "reghdfe##acknowledgements"}{...}
{viewerjumpto "References" "reghdfe##references"}{...}
{title:Title}

{p2colset 5 18 20 2}{...}
{p2col :{cmd:reghdfe} {hline 2}}Linear and instrumental-variable/GMM regression absorbing any number of fixed effects{p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{p 8 15 2} {cmd:reghdfe}
{depvar} [{indepvars}] [{cmd:(}{it:{help varlist:endogvars}} {cmd:=} {it:{help varlist:iv_vars}}{cmd:)}]
{ifin} {it:{weight}} {cmd:,} {opth a:bsorb(reghdfe##absvar:absvars)} [{help reghdfe##options:options}] {p_end}

{marker opt_summary}{...}
{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Model {help reghdfe##opt_model:[+]}}
{p2coldent:* {opt a:bsorb}{cmd:(}{help reghdfe##absvar:absvar} [...]{cmd:)}   }identifiers of the fixed effects that will be absorbed{p_end}
{synopt: {cmdab:a:bsorb(}{it:...}{cmd:,} {cmdab:save:fe)}}save all fixed effect estimates ({it:__hdfe*} prefix); useful for a subsequent {help reghdfe##postestimation:predict}.
However, see also the {it:resid} option.{p_end}
{synopt :{opth su:mmarize(tabstat##statname:stats)}}equivalent to {help reghdfe##postestimation:estat summarize} after the regression,
but more flexible, compatible with the {opt fast:} option, and saves the matrix of results on {it:e(summarize)}{p_end}
{synopt : {opt sub:options(...)}}additional options that will be passed to the regression command (either {help regress}, {help ivreg2}, or {help ivregress}){p_end}
{synopt : {opth res:iduals(newvar)}}save residuals; more direct and much faster than saving the fixed effects and then running predict{p_end}

{syntab:SE/Robust {help reghdfe##opt_vce:[+]}}
{p2coldent:+ {opt vce}{cmd:(}{help reghdfe##vcetype:vcetype}[, {it:opt}]{cmd:)}}{it:vcetype}
is {opt un:adjusted}/{opt ols} (default), {opt r:obust}, or {opt cl:uster} {it:clustervars};
{it:opt} may be {opt bw(#)}, {opt dkraay(#)}, {opt ker:nel(str)}, {opt kiefer}{p_end}

{syntab:IV/2SLS/GMM {help reghdfe##opt_iv:[+]}}
{synopt :{opt est:imator(str)}}estimator used in the instrumental-variable regression.
The default is {opt 2sls}, valid options are {opt gmm:2s} (two-stage GMM estimator),
{opt liml:} and {opt cue:} (which gives approximate results, as discussed in detail below){p_end}
{synopt :{opt stage:s(stage_list)}}runs and saves additional or alternative regression stages. Possible stages are: {it:ols first acid reduced} (and {it:all}).{p_end}
{synopt :{opth iv:suite(subcmd)}}package used in the IV/GMM regressions;
options are {opt ivreg2} (default; needs installing) and {opt ivregress}{p_end}
{synopt :{opt first}}report first stage regression (but sadly not first-stage summary results){p_end}
{synopt :{opt savefirst}}saves the first-stage regressions results; requires {opt first}{p_end}
{synopt :{opt showraw}}show the raw output of ivreg2 (if that's the ivsuite used); useful to see first-stage summary results{p_end}

{syntab:Diagnostic {help reghdfe##opt_diagnostic:[+]}}
{synopt :{opt v:erbose(#)}}amount of debugging information to show (0=None, 1=Some, 2=More, 3=Parsing/convergence details, 4=Every iteration){p_end}
{synopt :{opt time:it}}show elapsed times by stage of computation{p_end}

{syntab:Optimization {help reghdfe##opt_optimization:[+]}}
{p2coldent:+ {opth tol:erance(#)}}criterion for convergence (default=1e-7){p_end}
{synopt :{opth maxit:erations(#)}}maximum number of iterations (default=10000). if set to 0, it will run for as long as it takes.{p_end}
{synopt :{opth group:size(#)}}apply the within algorithm in groups of {it:#} variables (default 10). a large groupsize is usually faster but uses more memory{p_end}
{synopt :{opt accel:eration(str)}}acceleration method; options are conjugate_gradient (cg), steep_descent (sd), aitken (a), and none (no){p_end}
{synopt :{opt transf:orm(str)}}transform operation that defines the type of alternating projection; options are Kaczmarz (kac), Cimmino (cim), Symmetric Kaczmarz (sym){p_end}
{synopt :{opt check}}if convergence was achieved, the fixed effects should have a 1.0 coeficient in each step{p_end}

{syntab:Speedup Tricks {help reghdfe##opt_speedup:[+]}}
{synopt :{opt fast}}will not create {it:e(sample)}; disabled when saving fixed effects or mobility groups{p_end}
{synopt :{opt save:cache}}compute within transformation but do not regress the variables; useful when comparing alternative specifications (combine it with preserve/restore){p_end}
{synopt :{opth keep:vars(varlist)}}additional variables that will be kept in the demeaned dataset{p_end}
{synopt :{opt use:cache}}required with data previously transformed by {opt save:cache} (limitation: if you use clusters, you need to set the same vce(cluster ...) in both savecache and usecache){p_end}
{p2coldent:X {opt by(groupvar)}}similar to {opt save:cache} but will run the transformations independently for each level/category of {it:varname}. Stores the levels in {it:e(levels)}{p_end}
{p2coldent:X {opt l:evel(value)}}equivalent to regressing "{cmd:if} {it:groupvar}{cmd:==}{it:level}" but faster; needs to be run after {opt by(groupvar)}{p_end}
{p2coldent:X {opt nested}}add each {it:absvar} recursively, reporting the R2 and associated F-test
at each stage (only with ols and unadjusted standard errors){p_end}

{syntab:Degrees-of-Freedom Adjustments {help reghdfe##opt_dof:[+]}}
{synopt :{opt dof:adjustments(list)}}allows selecting the desired adjustments for degrees of freedom;
rarely used{p_end}
{synopt: {opth group(newvar)}}unique identifier for the first mobility group{p_end}

{syntab:Reporting {help reghdfe##opt_reporting:[+]}}
{synopt :{opt version:}}reports the version number and date of reghdfe, and saves it in e(version). standalone option{p_end}
{synopt :{opt l:evel(#)}}set confidence level; default is {cmd:level(95)}{p_end}
{synopt :{it:{help reghdfe##display_options:display_options}}}control column formats, row spacing, line width, display of omitted variables and base and empty cells, and factor-variable labeling{p_end}

{syntab:Undocumented}
{synopt :{opt old}}will call the latest 2.x version of reghdfe instead (see {help reghdfe_old}){p_end}
{synopt :{opt keepsin:gletons}}{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}* {opt absorb(absvars)} is required.{p_end}
{p 4 6 2}+ indicates a recommended or important option.{p_end}
{p 4 6 2}X indicates option is still incomplete or not fully tested; proceed with care.{p_end}
{p 4 6 2}{it:indepvars}, {it:endogvars} and {it:iv_vars} may contain factor variables; see {help fvvarlist}.{p_end}
{p 4 6 2}all the regression variables may contain time-series operators; see {help tsvarlist}.{p_end}
{p 4 6 2}{cmd:fweight}s, {cmd:aweight}s and {cmd:pweight}s are allowed; see {help weight}.{p_end}


{marker absvar}{...}
{title:Absvar Syntax}

{synoptset 22}{...}
{synopthdr:absvar}
{synoptline}
{synopt:{cmd:i.}{it:varname}}categorical variable to be absorbed (the {cmd:i.} prefix is tacit){p_end}
{synopt:{cmd:i.}{it:var1}{cmd:#i.}{it:var2}}absorb the interactions of multiple categorical variables{p_end}
{synopt:{cmd:i.}{it:var1}{cmd:#}{cmd:c.}{it:var2}}absorb "fixed slopes", where {it:var2} has a different slope coef. depending on the category of {it:var1}{p_end}
{synopt:{it:var1}{cmd:##}{cmd:c.}{it:var2}}equivalent to "{cmd:i.}{it:var1} {cmd:i.}{it:var1}{cmd:#}{cmd:c.}{it:var2}", but {it:much} faster{p_end}
{synopt:{it:var1}{cmd:##c.(}{it:var2 var3}{cmd:)}}multiple fixed slopes are allowed together. Alternative syntax: {it:var1}{cmd:##(c.}{it:var2} {cmd:c.}{it:var3}{cmd:)}{p_end}
{synopt:{it:v1}{cmd:#}{it:v2}{cmd:#}{it:v3}{cmd:##c.(}{it:v4 v5}{cmd:)}}factor operators can be combined{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}To save the estimates of only some absvars, write {newvar}{inp:={it:absvar}}.{p_end}
{p 4 6 2}Please be aware that in most cases these estimates are neither consistent nor econometrically identified.{p_end}
{p 4 6 2}Using categorical interactions (e.g. {it:x}{cmd:#}{it:z}) is faster than running {it:egen group(...)} beforehand.{p_end}
{p 4 6 2}Singleton obs. are dropped iteratively until no more singletons are found (see ancilliary article for details).{p_end}

{bf:---> NOTE: below this line, the help file has not been updated for version 3 <---}

{marker description}{...}
{title:Description}

{pstd}
{cmd:reghdfe} fits a linear regression of {depvar} on {indepvars} while absorbing an arbitrary number of fixed effects indicated by the categories of {help reghdfe##absvar:absvars}. It also supports regressing
 on {it:endogvars}, in which case it uses {it:iv_vars} (along with {it:indepvars} and the fixed effects) as instruments for {it:endogvars} (either with 2SLS or GMM/LIML).{p_end}
 
{pstd}No constant is reported as that is absorbed by the fixed effects. However, if you wish to recover it,
after the regression run {cmd:predict dummies, d} followed by {cmd:su dummies}, {cmd:di r(mean)}.{p_end}

{pstd}The estimates for the fixed effects (including those with continous interactions) can be saved, although their standard errors are not recovered. When using multiple highly-dimensional fixed-effects,
the user should be aware of the identification requirements regarding the fixed effects. For instance, the fixed effects cannot form disjoint graphs or else identification would only be possible
within each subgraph (see {help reghdfe##references:references}).{p_end}

{pstd}There are several features generalized from either {cmd: areg} or {cmd: xtreg, fe}, such as:{p_end}

{p2col 6 9 9 2: a)}reporting F-tests on the absorbed variables (except with robust and clustered vce){p_end}
{p2col 6 9 9 2: b)}correlations between the FEs and the other regressors(except with robust and clustered vce){p_end}
{p2col 6 9 9 2: c)}degrees-of-freedom adjustements with clustered data when one absorbed category is contained within the clusters (as in {cmd: xtreg, fe robust}){p_end}

{marker options}{...}
{title:Options}

{marker opt_model}{...}
{dlgtab:Model and Miscellanea}

{phang}
{opth a:bsorb(reghdfe##absvar:absvars)} list of categorical variables (or interactions) representing the fixed effects to be absorbed.
this is equivalent to including a dummy variable for each category of each {it:absvar}. {cmd:absorb()} is required.

{pmore}
To save a fixed effect, prefix the absvar with "{newvar}{cmd:=}".
For instance, the option {input:absorb(firm_id worker_id year_coefs=year_id)} will include firm,
worker and year fixed effects, but will only save the estimates for the year fixed effects (in the new variable {it:year_coefs}).

{pmore}
If you want to {help reghdfe##postestimation:predict} afterwards but don't care about setting the names of each fixed effect, use the {cmdab:save:fe} suboption.
This will delete all variables named {it:__hdfe*__} and create new ones as required.
Example: {it:reghdfe price weight, absorb(turn trunk, savefe)}

{phang}
{opt dropsi:ngletons} will drop the observations that contain values of the fixed effects repeated only once.
For instance, with individual fixed effects, it will drop individuals that appear only once in the sample.
{hi: I strongly recommend you to use this option}.

{pmore}
Benefits of this option include a faster estimation (due to less observations and fixed effects) and
less change of problems when estimating the VCE matrix.

{pmore}
Note that after dropping singletons on the first fixed effect, a value of the second fixed effect that
previously appeared in two observations may now appear in only one and thus will be dropped.
Also note that only one check will be done for each fixed effect.

{phang}
{opt nested} Add each {it:absvar} recursively into the model

{pmore} This option will reporting the R2 at each stage, and also compute the Fisher test of significance
for each set of absorbed variables.

{pmore} Only available in OLS with {it:vce(unadjusted)}.

{phang}
{opth su:mmarize(tabstat##statname:stats)} will report and save a table of summary of statistics of the regression
variables (including the instruments, if applicable), using the same sample as the regression.

{pmore} {opt su:mmarize} (without parenthesis) saves the default set of statistics: {it:mean min max}.

{pmore} The complete list of accepted statistics is available in the {help tabstat##statname:tabstat help}. The most useful are {it:count range sd median p##}.

{pmore} The summary table is saved in {it:e(summarize)}

{pmore} To save the summary table silently (without showing it after the regression table), use the {opt qui:etly} suboption. You can use it by itself ({cmd:summarize(,quietly)}) or with custom statistics ({cmd:summarize(mean, quietly)}).

{phang}
{opt sub:options(...)}
options that will be passed directly to the regression command (either {help regress}, {help ivreg2}, or {help ivregress})

{pmore}
Some options are not allowed and will be silently ignored ({it:nosmall}, {it:hascons}, {it:tsscons})

{marker opt_vce}{...}
{dlgtab:SE/Robust}

{phang}
{opth vce:(reghdfe##vcetype:vcetype, subopt)}
specifies the type of standard error reported.
Note that all the advanced estimators rely on asymptotic theory, and will likely have poor performance with small samples
(but again if you are using reghdfe, that is probably not your case)

{pmore}
{opt un:adjusted}/{opt ols:} estimates conventional standard errors, valid even in small samples
under the assumptions of homoscedasticity and no correlation between observations

{pmore}
{opt r:obust} estimates heteroscedasticity-consistent standard errors (Huber/White/sandwich estimators), but still assuming independence between observations

{pmore}Warning: in a FE panel regression, using {opt r:obust} will
lead to inconsistent standard errors if for every fixed effect, the {it:other} dimension is fixed.
For instance, in an standard panel with individual and time fixed effects, we require both the number of
individuals and time periods to grow asymptotically.
If that is not the case, an alternative may be to use clustered errors,
which as discussed below will still have their own asymptotic requirements.
For a discussion, see
{browse "http://www.princeton.edu/~mwatson/papers/ecta6489.pdf":Stock and Watson, "Heteroskedasticity-robust standard errors for fixed-effects panel-data regression," Econometrica 76 (2008): 155-174}

{pmore}
{opt cl:uster} {it:clustervars} estimates consistent standard errors even when the observations
are correlated within groups.

{pmore}
Multi-way-clustering is allowed. Thus, you can indicate as many {it:clustervar}s as desired
(e.g. allowing for intragroup correlation across individuals, time, country, etc.).

{pmore}
Each {it:clustervar} permits interactions of the type {it:var1{cmd:#}var2}
(this is faster than using {cmd:egen group()} for a one-off regression).

{pmore} Warning: The number of clusters, for all of the cluster variables, must go off to infinity.
A frequent rule of thumb is that each cluster variable must have at least 50 different categories
(the number of categories for each clustervar appears on the header of the regression table).

{pstd}
The following suboptions require either the {help ivreg2} or the {help avar} package from SSC.
For a careful explanation, see the {help ivreg2##s_robust:ivreg help file}, from which the comments below borrow.

{pmore}
{opt u:nadjusted}{cmd:, }{opt bw(#)} (or just {cmd:, }{opt bw(#)}) estimates autocorrelation-consistent standard errors (Newey-West).

{pmore}
{opt r:obust}{cmd:, }{opt bw(#)} estimates autocorrelation-and-heteroscedasticity consistent standard errors (HAC).

{pmore}
{opt cl:uster} {it:clustervars}{cmd:, }{opt bw(#)} estimates standard errors consistent to common autocorrelated disturbances (Driscoll-Kraay). At most two cluster variables can be used in this case.

{pmore}
{cmd:, }{opt kiefer} estimates standard errors consistent under arbitrary intra-group autocorrelation (but not heteroskedasticity) (Kiefer).

{pmore}
{opt kernel(str)} is allowed in all the cases that allow {opt bw(#)}
The default kernel is {it:bar} (Bartlett). Valid kernels are Bartlett (bar); Truncated (tru); Parzen (par);
Tukey-Hanning (thann); Tukey-Hamming (thamm); Daniell (dan); Tent (ten); and Quadratic-Spectral (qua or qs). 

{pstd}
Advanced suboptions:

{pmore}
{cmd:, }{opt suite(default|mwc|avar)} overrides the package chosen by reghdfe to estimate the VCE.
{it:default} uses the default Stata computation (allows unadjusted, robust, and at most one cluster variable).
{it:mwc} allows multi-way-clustering (any number of cluster variables), but without the {it:bw} and {it:kernel} suboptions.
{it:avar} uses the avar package from SSC. Is the same package used by ivreg2, and allows the {it:bw}, {it:kernel}, {it:dkraay} and {it:kiefer} suboptions.
This is useful almost exclusively for debugging.

{pmore}
{cmd:, }{opt twice:robust} will compute robust standard errors not only on the first but on the second step of the gmm2s estimation. Requires {opt ivsuite(ivregress)}, but will not give the exact same results as ivregress.

{pmore}{it:Explanation:} When running instrumental-variable regressions with the {cmd:ivregress} package,
robust standard errors, and a gmm2s estimator,
{opt vce(robust)} is translated into {opt wmatrix(robust)} {opt vce(unadjusted)}.
This maintains compatibility with {cmd:ivreg2} and other packages, but may unadvisable as described in {help ivregress} (technical note). Specifying this option will instead use {opt wmatrix(robust)} {opt vce(robust)}.

{pmore}However, computing the second-step vce matrix requires computing updated estimates (including updated fixed effects).
Since reghdfe currently does not allow this, the resulting standard errors
{hi:will not be exactly the same as with ivregress}.
This issue is similar to applying the CUE estimator, described further below.

{pmore}Note: The above comments are also appliable to clustered standard error.

{pstd}
Options {opt boot:strap} and {opt jack:knife} could be implemented, although their execution would be extremely slow.


{marker opt_iv}{...}
{dlgtab:IV/2SLS/GMM}

{phang}
{opt est:imator}{cmd:(}{opt 2sls}|{opt gmm:2s}|{opt liml}|{opt cue}{cmd:)}
estimator used in the instrumental-variable estimation

{pmore}
{opt 2sls} (two-stage least squares, default), {opt gmm:2s} (two-stage efficient GMM), {opt liml} (limited-information maximum likelihood), and
{opt cue} ("continuously-updated" GMM) are allowed.{p_end}

{pmore}
Warning: {opt cue} will not give the same results as ivreg2. See the discussion in
{browse "http://www.stata-journal.com/sjpdf.html?articlenum=st0030_3": Baum, Christopher F., Mark E. Schaffer, and Steven Stillman. "Enhanced routines for instrumental variables/GMM estimation and testing." Stata Journal 7.4 (2007): 465-506}
(page 484).
Note that even if this is not exactly {opt cue}, it may still be a desirable/useful alternative to standard cue, as explained in the article.

{phang}
{opt stage:s(stage_list)}
adds and saves up to four auxiliary regressions useful when running instrumental-variable regressions:

{phang2}{cmd:ols} an equivalent ols regression (as a benchmark/comparison){p_end}
{phang2}{cmd:first} all the first stage regressions{p_end}
{phang2}{cmd:acid} an "acid" regression that includes both the instruments, exogenous regressors (i.e. included instruments), and endogenous regressors, against the dependent variable.
In this setup, the instruments should not be significant.{p_end}
{phang2}{cmd:reduced} the reduced-form regression, between the instruments and exogenous regressors, and the dependent variable (excludes the endogenous regressor).{p_end}

{phang}
{opth iv:suite(subcmd)}
allows the IV/2SLS regression to be run either using {opt ivregress} or {opt ivreg2}.

{pmore} {opt ivreg2} is the default, but needs to be installed for that option to work.

{phang}
{opt first}
will report first stage regression.

{pmore}
Note that first-stage summary results will not be reported due to the way {it:ivreg2} is called; to view them see the option {it:showraw}

{phang}
{opt savefirst}
will save first stage regressions. Requires the option {opt first} to work.

{phang}
{opt showraw}
shows the entire output of ivreg2 (if that's the ivsuite used); this is used to see the first-stage summary results.

{pmore}
The downside is that it will have temporary names in place of whenever factors variables and time-series operators are used.

{marker opt_diagnostic}{...}
{dlgtab:Diagnostic}

{phang}
{opt v:erbose(#)} orders the command to print debugging information.

{pmore}
Possible values are 0 (none), 1 (some information), 2 (even more), 3 (adds dots for each iteration, and reportes parsing details), 4 (adds details for every iteration step)

{pmore}
For debugging, the most useful value is 3. For simple status reports, set verbose to 1.

{phang}
 {opt check} will regress the variable against the calculated fixed effects. If convergence was indeed achieved, the coefficients should be 1.0 in each step (except under perfect collinearity situations.)
 
{marker opt_dof}{...}
{dlgtab:Degrees-of-Freedom Adjustments}

{phang}
{opt dof:adjustments(doflist)} selects how the degrees-of-freedom, as well as e(df_a), are adjusted due to the absorbed fixed effects.

{pmore}
Without any adjustment, we would assume that the degrees-of-freedom used by the fixed effects is equal to the count of all the fixed effects
(e.g. number of individuals + number of years in a typical panel).
However, in complex setups (e.g. fixed effects by individual, firm, job position, and year),
there may be a huge number of fixed effects collinear with each other, so we want to adjust for that.

{pmore}
Note: changing the default option is rarely needed, except in benchmarks, and to obtain a marginal speed-up by excluding the {opt pair:wise} option.

{pmore}
{opt all} is the default and almost always the best alternative. It is equivalent to {opt dof(pairwise clusters continuous)}

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
One solution is to ignore subsequent fixed effects (and thus oversestimate e(df_a) and understimate the degrees-of-freedom).
Another solution, described below, applies the algorithm between pairs of fixed effects to obtain a better (but not exact) estimate:

{pmore}
{opt pair:wise} applies the aforementioned connected-subgraphs algorithm between pairs of fixed effects.
For instance, if there are four sets of FEs, the first dimension will usually have no redundant coefficients (i.e. e(M1)==1), since we are running the model without a constant.
For the second FE, the number of connected subgraphs with respect to the first FE will provide an exact estimate of the degrees-of-freedom lost, e(M2).

{pmore}
For the third FE, we do not know exactly.
However, we can compute the number of connected subgraphs between the first and third {it:G(1,3)},
and second and third {it:G(2,3)} fixed effects, and choose the higher of those as the closest estimate for e(M3).
For the fourth FE, we compute {it:G(1,4)}, {it:G(2,4)} and {it:G(3,4)} and again choose the highest for e(M4).

{pmore}
Finally, we compute e(df_a) = e(K1) - e(M1) + e(K2) - e(M2) + e(K3) - e(M3) + e(K4) - e(M4);
where e(K#) is the number of levels or dimensions for the #-th fixed effect (e.g. number of individuals or years).
Note that e(M3) and e(M4) are only conservative estimates and thus we will usually be overestimating the standard errors. However, given the sizes of the datasets typically used with reghdfe, the difference should be small. 

{pmore}
Since the gain from {opt pair:wise} is usually {it:minuscule} for large datasets, and the computation is expensive, it may be a good practice to exclude this option for speedups.

{pmore}
{opt cl:usters}
will check if a fixed effect is nested within a {it:clustervar}.
In that case, it will set e(K#)==e(M#) and no degrees-of-freedom will be lost due to this fixed effect.
The rationale is that we are already assuming that the number of effective observations is the number of cluster levels.
This is the same adjustment that {cmd:xtreg, fe} does, but {cmd:areg} does not use it.

{pmore}
{opt cont:inuous}
Fixed effects with continuous interactions (i.e. individual slopes, instead of individual intercepts) are dealt with differently.
In an i.categorical#c.continuous interaction, we will do one check: we count the number of categories where c.continuous is always zero.
In an i.categorical##c.continuous interaction, we do the above check but replace zero for any particular constant.
In the case where continuous is constant for a level of categorical, we know it is collinear with the intercept, so we adjust for it.

{pmore}
Additional methods, such as {opt bootstrap} are also possible but not yet implemented.
Some preliminary simulations done by the author showed a very poor convergence of this method.

{phang}
{opth group(newvar)} name of the new variable that will contain the first mobility group.
Requires {opt pair:wise}, {opt first:pair}, or the default {opt all}.

{marker opt_speedup}{...}
{dlgtab:Speeding Up Estimation}

{phang}
{ul:Parallel Computing}

{phang}
{opt fast} avoids one {it:save}, one {it:use}, and one {it:merge} operation.
It usually makes the regression 30% faster, but comes at a cost:
you cannot save the FEs with predict, you cannot predict the residuals, e(sample) is not saved, and the correlation between the FEs and xb are not reported.
Also, savecache, usecache, nested, check, group options are not compatible with fast.

{pmore}
If you wish to use {opt fast} while reporting {cmd:estat summarize}, see the {opt summarize} option.

{phang}
{opth cores(#)} will run the demeaning algorithm in # parallel instances.

{pmore}
currently, cores() only runs on Windows, although the code is already generalized for Linux and OSX. Volunteers for bugfixing or testing are welcome!

{pmore}
Several Stata processes will be created, and the task of demeaning all the required variables will be distributed amongst them.
This option requires the package {help parallel:parallel} by George Vega Yon (run {it:ssc install parallel} to download it)

{pmore}
Disclaimer: there may still be some rough corners (e.g. sometimes not deleting temporary files)

{pmore}Example:{p_end}
{phang2}{cmd:. sysuse auto}{p_end}
{phang2}{cmd:. reghdfe price weight length, a(turn rep) cores(2)}{p_end}

{phang}
{ul:Precomputing transformations} (useful when testing alternative specifications)

{phang}
{cmd:reghdfe} {varlist} {ifin}{cmd:,} {opt a:bsorb(absvars)} {opth save:cache(filename)} [{it:options}]

{pmore}
This will demean {it:varlist} with respect to {it:absvars}, and save the transformed variables in {it:filename}. Note that if any variable has a missing value, the entire row is dropped.

{pmore}
Options allowed are all optimization options except {it:fast}, and all diagnostic options including {it:cores(#)}

{pmore}
This will create a variable __uid__ in the master data.

{phang}
{opth use:cache(filename)} can be added to a normal command and will load the transformed variables from {it:filename} instead of computing them again.

{pmore}
In order for this to work, i) the filename needs to exist and all the variables need to have been precomputed using that filename,
ii) the variable __uid__ must exist, iii) the precomputed transformation must have been made using the same {it:absvars} and the same observations 
(if it was precomputed with more observations, an error will occur; with less observations, that subset of the dataset will be used for the regression).

{pmore}Example:{p_end}
{phang2}{cmd:. sysuse auto}{p_end}
{phang2}{cmd:. tempfile cache}{p_end}
{phang2}{cmd:. reghdfe price weight length, a(turn rep) savecache("`cache'")}{p_end}
{phang2}{cmd:. reghdfe price weight, a(turn rep) usecache("`cache'")}{p_end}
{phang2}{cmd:. reghdfe price length, a(turn rep) usecache("`cache'")}{p_end}

{phang}
{ul:Running the same regression over different categories of a variable} (equivalent to by:)

{phang} {opth over(varname)} allows regressions over the different levels of {it:varname}, with the advantage that the data needs to be demeaned only once.

{pmore}
To use this, first add this option together with {it:savecache}. This will change the {it:absvars} from e.g. "i.var1 i.var2##c.var3" to "i.over i.over#i.var1 i.over#i.var2##c.var3", 
basically adding one fixed effect and adding the {it:over} variable to every interaction. Transforming the entire dataset with this specification is equivalent to transforming it separately by levels of {it:over}.

{pmore}
This call will return {opt e(over_levels)}, after which the regressions can be called as long as both {it:usecache} and {it:over} are specified. 
The user should be very careful not to change the dataset between calls, as the program will detect only some inconsistencies between the {it:savecache} call and the {it:usecache} call.

{pmore}Example:{p_end}
{phang2}{cmd:. sysuse auto}{p_end}
{phang2}{cmd:. tempfile cache}{p_end}
{phang2}{cmd:. reghdfe price weight length, a(turn rep) savecache("`cache'") over(foreign)}{p_end}
{phang2}{cmd:. local levels `e(levels_over)'}{p_end}
{phang2}{cmd:. foreach level of local levels {c -(}}{p_end}
{phang3}{cmd:. reghdfe price weight length if foreign==`level', a(turn rep) usecache("`cache'") over(foreign)}{p_end}
{phang2}{cmd:. {c )-}}{p_end}

{pmore}This is equivalent to running{p_end}
{phang2}{cmd:. reghdfe price weight length  if foreign==0, a(turn rep)}{p_end}
{phang2}{cmd:. reghdfe price weight length  if foreign==1, a(turn rep)}{p_end}

{marker opt_optimization}{...}
{dlgtab:Optimization}

{phang}
{opth tol:erance(#)} specifies the tolerance criterion for convergence; default is {cmd:tolerance(1e-7)}

{pmore}
Note that for tolerances beyond 1e-14, the limits of the {it:double} precision are reached and the results will most likely not converge.

{pmore}
At the other end, is not tight enough, the regression may not identify perfectly collinear regressors. However, those cases can be easily spotted due to their extremely high standard errors.

{phang}
{opth maxit:erations(#)}
specifies the maximum number of iterations; the default is {cmd:maxiterations(10000)}; 0 means run forever until convergence.

{phang}
{it:Advanced options:}

{phang}
{opt fast} avoids one {it:save}, one {it:use}, and one {it:merge} operation.
Useful if the dataset is very large, but will not save {it:e(sample)} or compute correlations between the fixed effects and Xb.

{pmore} Will not work under {cmd:check}, or if any variable (fixed effect, mobility group) needs to be saved.

{phang}
{opt noaccel:erate} apply fixed point iteration without applying Aitken's acceleration.

{pmore}
This will be much slower but may avoid situations where the acceleration gets stuck.

{phang}
{opth accel_start(#)} indicates how many iterations to wait until the Aitken's acceleration starts. Default is 6.

{phang}
{opth accel_freq(#)} indicates how often the acceleration occurs. The default is 3, meaning that after two steps without acceleration, the third one accelerates. Other common value is 6.

{phang}
{opth bad_loop_threshold(#)} If the acceleration seems stuck # times in a row, pause it. Default is 1.

{phang}
{opth stuck_threshold(#)} Defines when is the acceleration stuck (if the relative improvement is less than #). Default is 5e-3

{phang}
{opth pause_length(#)} indicates how many acceleration steps to pause after the iteration got stuck. Default is 20

{pmore}
Note that this is in terms of acceleration steps, not iterations (so if accel_freq=3 and pause_length=20, 60 iterations will pass until the acceleration resumes)

{marker opt_reporting}{...}
{dlgtab:Reporting}

{phang}
{opt l:evel(#)} sets confidence level; default is {cmd:level(95)}

{marker display_options}{...}
{phang}
{it:display_options}:
{opt noomit:ted},
{opt vsquish},
{opt noempty:cells},
{opt base:levels},
{opt allbase:levels},
{opt nofvlabel},
{opt fvwrap(#)},
{opt fvwrapon(style)},
{opth cformat(%fmt)},
{opt pformat(%fmt)},
{opt sformat(%fmt)}, and
{opt nolstretch};
    see {helpb estimation options##display_options:[R] estimation options}.
    {p_end}
    
{marker postestimation}{...}
{title:Postestimation Syntax}

Only {cmd:estat summarize}, {cmd:predict} and {cmd:test} are currently supported and tested.

{p 8 13 2}
{cmd:estat summarize}
{p_end}{col 23}Summarizes {it:depvar} and the variables described in {it:_b} (i.e. not the excluded instruments)

{p 8 16 2}
{cmd:predict} 
{newvar} 
{ifin}
[{cmd:,} {it:statistic}]
{p_end}{col 23}Requires all set of fixed effects to be previously saved by {cmd:reghdfe} (except for option {opt xb})
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
{synoptline}
{p2colreset}{...}

{p 8 13 2}
{cmd:test}
{p_end}{col 23}Performs significance test on the parameters, see the {help test:stata help}

{p 8 13 2}
{cmd:suest}
{p_end}

{pstd}
If you want to perform tests that are usually run with {cmd:suest},
such as non-nested models, tests using alternative specifications of the variables,
or tests on different groups, you can replicate it manually, as described 
{browse "http://www.stata.com/statalist/archive/2009-11/msg01485.html":here}.
{p_end}

{pstd}Note: do not use {cmd:suest}. It will run, but the results will be incorrect.{p_end}

{marker remarks}{...}
{title:Implementation Details}

{p2col 5 7 7 2: -}This program usually runs at least 10 times faster than related programs ({cmd:areg}, {cmd:xtreg, fe}, {cmd:twfe}, {cmd:a2reg}, {cmd:reg2hdfe}, etc.).{p_end}
{p2col 5 7 7 2: -}The relative gain increases with the number of observations
 and the number of absorbed fixed effects, so for a small dataset the gain could even 
  be negative due to the initial setup of the program).{p_end}
  {p2col 5 7 7 2: -}The reason is that it takes the best of these programs 
 (such as the fixed point iteration from reg2hdfe) while improving on their bottlenecks.{p_end}
{p2col 5 7 7 2: -}In particular, the biggest gain comes from computing averages by group
 without sorting the dataset by each group. Since sorting is an o(n log n) operation, 
  and averages by group is o(n), the gains get larger with the dataset size.{p_end}
  {p2col 5 7 7 2: -}It also allows for efficient estimation of interactions between continuous variables and HDFEs{p_end}
 {p2col 5 7 7 2: -}Finally, also allows for running IV/GMM regressions{p_end}
{p2col 5 7 7 2: -}For more details see the PDF notes.{p_end}

{title:Possible Pitfalls and Common Mistakes}

{p 5 8 2}1. (note: as of version 2.1, the constant is no longer reported) Ignore the constant; it doesn't tell you much. If you want to use descriptive stats, that's what the {opt sum:marize()} and {cmd:estat summ} commands are for.
Even better, use {opt noconstant} to drop it (although it's not really dropped as it never existed on the first place!){p_end}
{p 5 8 2}2. Think twice before saving the fixed effects. They are probably inconsistent / not identified and you will likely be using them wrong.{p_end}
{p 5 8 2}3. It's good practice to drop singletons. {opt dropsi:ngleton} is your friend.{p_end}
{p 5 8 2}4. If you use {opt vce(robust)}, be sure that your {it:other} dimension is not "fixed" but grows with N, or your SEs will be wrong.{p_end}
{p 5 8 2}5. If you use {opt vce(cluster ...)}, check that your number of clusters is high enough (50+ is a rule of thumb). If not, you are making the SEs even worse!{p_end}
{p 5 8 2}6. The panel variables (absvars) should probably be nested within the clusters (clustervars) due to the within-panel correlation induced by the FEs.
(this is not the case for *all* the absvars, only those that are treated as growing as N grows){p_end}
{p 5 8 2}7. If you run analytic or probability weights,
you are responsible for ensuring that the weights stay
constant within each unit of a fixed effect (e.g. individual),
or that it is correct to allow varying-weights for that case.
{p_end}
{p 5 8 2}8. Be aware that adding several HDFEs is not a panacea.
The first limitation is that it only uses within variation (more than acceptable if you have a large enough dataset).
The second and subtler limitation occurs if the fixed effects are themselves outcomes of the variable of interest (as crazy as it sounds).
For instance, imagine a regression where we study the effect of past corporate fraud on future firm performance.
We add firm, CEO and time fixed-effects (standard practice). This introduces a serious flaw: whenever a fraud event is discovered,
i) future firm performance will suffer, and ii) a CEO turnover will likely occur.
Moreover, after fraud events, the new CEOs are usually specialized in dealing with the aftershocks of such events
(and are usually accountants or lawyers).
The fixed effects of these CEOs will also tend to be quite low, as they tend to manage firms with very risky outcomes.
Therefore, the regressor (fraud) affects the fixed effect (identity of the incoming CEO).
Adding particularly low CEO fixed effects will then overstate the performance of the firm,
and thus {it:understate} the negative effects of fraud on future firm performance.{p_end}

{title:Missing Features}

{phang}(If you are interested in discussing these or others, feel free to {help reghdfe##contact:contact me})

{p2col 5 7 7 2: -}Add a more thorough discussion on the possible identification issues{p_end}
{p2col 5 7 7 2: -}Find out a way to use reghdfe iteratively with CUE
(right now only OLS/2SLS/2SGMM/LIML give the exact same results){p_end}
{p2col 5 7 7 2: -}Implement a -bootstrap- option in DoF estimation{p_end}
{p2col 5 7 7 2: -}The interaction with cont vars (i.a#c.b) may suffer from numerical accuracy issues, as we are dividing by a sum of squares{p_end}
{p2col 5 7 7 2: -}Calculate exact DoF adjustment for 3+ HDFEs (note: not a problem with cluster VCE when one FE is nested within the cluster){p_end}
{p2col 5 7 7 2: -}More postestimation commands (lincom? margins?){p_end}
{p2col 5 7 7 2: -}Not sure if I should add an F-test for the absvars in the vce(robust) and vce(cluster) cases.
Discussion on e.g. -areg- (methods and formulas) and textbooks suggests not;
on the other hand, there may be alternatives:
{it:{browse "http://www.socialsciences.manchester.ac.uk/disciplines/economics/research/discussionpapers/pdf/EDP-1124.pdf" :A Heteroskedasticity-Robust F-Test Statistic for Individual Effects}}{p_end}

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

{pstd}Report nested F-tests{p_end}
{phang2}{cmd:. reghdfe ln_w grade age ttl_exp tenure not_smsa south , absorb(idcode year) nested}{p_end}

{pstd}Do AvgE instead of absorb() for one FE{p_end}
{phang2}{cmd:. reghdfe ln_w grade age ttl_exp tenure not_smsa south , absorb(idcode year) avge(occ)}{p_end}
{phang2}{cmd:. reghdfe ln_w grade age ttl_exp tenure not_smsa south , absorb(idcode year) avge(AvgByOCC=occ)}{p_end}

{pstd}Check that FE coefs are close to 1.0{p_end}
{phang2}{cmd:. reghdfe ln_w grade age ttl_exp tenure not_smsa , absorb(idcode year) check}{p_end}

{pstd}Save first mobility group{p_end}
{phang2}{cmd:. reghdfe ln_w grade age ttl_exp tenure not_smsa , absorb(idcode occ) group(mobility_occ)}{p_end}

{pstd}Factor interactions in the independent variables{p_end}
{phang2}{cmd:. reghdfe ln_w i.grade#i.age ttl_exp tenure not_smsa , absorb(idcode occ)}{p_end}

{pstd}Interactions in the absorbed variables (notice that only the {it:#} symbol is allowed){p_end}
{phang2}{cmd:. reghdfe ln_w grade age ttl_exp tenure not_smsa , absorb(idcode#occ)}{p_end}

{pstd}Interactions in both the absorbed and AvgE variables (again, only the {it:#} symbol is allowed){p_end}
{phang2}{cmd:. reghdfe ln_w grade age ttl_exp not_smsa , absorb(idcode#occ) avge(tenure#occ)}{p_end}

{pstd}IV regression{p_end}
{phang2}{cmd:. sysuse auto}{p_end}
{phang2}{cmd:. reghdfe price weight (length=head), absorb(rep78)}{p_end}
{phang2}{cmd:. reghdfe price weight (length=head), absorb(rep78) first}{p_end}
{phang2}{cmd:. reghdfe price weight (length=head), absorb(rep78) ivsuite(ivregress)}{p_end}

{pstd}Factorial interactions{p_end}
{phang2}{cmd:. reghdfe price weight (length=head), absorb(rep78)}{p_end}
{phang2}{cmd:. reghdfe price weight length, absorb(rep78 turn##c.price)}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:reghdfe} stores the following in {cmd:e()}:

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(N_hdfe)}}number of absorbed fixed-effects{p_end}
{synopt:{cmd:e(tss)}}total sum of squares{p_end}
{synopt:{cmd:e(df_m)}}model degrees of freedom{p_end}
{synopt:{cmd:e(rss)}}residual sum of squares{p_end}
{synopt:{cmd:e(df_r)}}residual degrees of freedom{p_end}
{synopt:{cmd:e(r2)}}R-squared{p_end}
{synopt:{cmd:e(r2_a)}}adjusted R-squared{p_end}
{synopt:{cmd:e(r2_within)}}Within R-squared{p_end}
{synopt:{cmd:e(r2_a_within)}}Adjusted Within R-squared{p_end}
{synopt:{cmd:e(df_a)}}degrees of freedom lost due to the fixed effects{p_end}
{synopt:{cmd:e(rmse)}}root mean squared error{p_end}
{synopt:{cmd:e(ll)}}log-likelihood{p_end}
{synopt:{cmd:e(ll_0)}}log-likelihood of fixed-effect-only regression{p_end}
{synopt:{cmd:e(N_clust)}}number of clusters{p_end}
{synopt:{cmd:e(F)}}F statistic{p_end}
{synopt:{cmd:e(F_absorb)}}F statistic for absorbed effect (when
        {cmd:vce(robust)} is not specified){p_end}
        {synopt:{cmd:e(rank)}}rank of {cmd:e(V)}{p_end}
        
{synopt:{cmd:e(K#)}}Number of categories of the #th absorbed FE{p_end}
{synopt:{cmd:e(M#)}}Number of redundant categories of the #th absorbed FE{p_end}
{synopt:{cmd:e(mobility)}}Sum of all {cmd:e(M#)}{p_end}
{synopt:{cmd:e(corr#)}}Correlation between #th absorbed FE and xb{p_end}

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:reghdfe}{p_end}
{synopt:{cmd:e(subcmd)}}either {cmd:regress}, {cmd:ivreg2} or {cmd:ivregress}{p_end}
{synopt:{cmd:e(model)}}Either ols or iv{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(dofmethod)}}dofmethod employed in the regression{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(indepvars)}}names of independent variables{p_end}
{synopt:{cmd:e(endogvars)}}names of endogenous right-hand-side variables{p_end}
{synopt:{cmd:e(instruments)}}names of excluded instruments{p_end}
{synopt:{cmd:e(absvars)}}name of the absorbed variables or interactions{p_end}
{synopt:{cmd:e(title)}}title in estimation output{p_end}
{synopt:{cmd:e(clustvar)}}name of cluster variable{p_end}
{synopt:{cmd:e(vce)}}{it:vcetype} specified in {cmd:vce()}{p_end}
{synopt:{cmd:e(vcetype)}}title used to label Std. Err.{p_end}
{synopt:{cmd:e(properties)}}{cmd:b V}{p_end}

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix of the estimators{p_end}

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}
{p2colreset}{...}

{marker contact}{...}
{title:Author}

{phang}
Sergio Correia

{phang}
Fuqua School of Business, Duke University

{phang}
Email: {browse "mailto:sergio.correia@duke.edu":sergio.correia@duke.edu}
{p_end}

{marker updates}{...}
{title:Latest Updates}

{pstd}
{cmd:reghdfe} is updated frequently, and upgrades or minor bug fixes may not be immediately available in SSC.
To check or contribute to the latest version of reghdfe, explore the
{browse "https://github.com/sergiocorreia/reghdfe":Github repository}.{p_end}

{pstd}
To see your current version and installed dependencies, type {cmd:reghdfe, version}
{p_end}

{marker acknowledgements}{...}
{title:Acknowledgements}

{pstd}
This package wouldn't have existed without the invaluable feedback and contributions of Paulo Guimaraes
and Amine Ouazad. I am also indebted to the guidance of Mark Schaffer, Kit Baum, and Nikolas Mittag;
and to the great bug-spotting abilities of many users.{p_end}

{pstd}In addition, {it:reghdfe} is build upon important contributions from the Stata community:{p_end}

{phang}{browse "https://ideas.repec.org/c/boc/bocode/s457101.html":reg2hdfe}, from Paulo Guimaraes,
and {browse "https://ideas.repec.org/c/boc/bocode/s456942.html":a2reg} from Amine Ouazad,
 were the inspiration and building blocks on which reghdfe was built.{p_end}
 
{phang}{browse "http://www.repec.org/bocode/i/ivreg2.html":ivreg2}, by Christopher F Baum, Mark E Schaffer and Steven Stillman, is the package used by default for instrumental-variable regression.{p_end}

{phang}{browse "https://ideas.repec.org/c/boc/bocode/s457689.html":avar} by Christopher F Baum and Mark E Schaffer, is the package used for estimating the HAC-robust standard errors of ols regressions.{p_end}

{phang}{browse "http://econpapers.repec.org/software/bocbocode/s456797.htm)":tuples} by Joseph Lunchman and Nicholas Cox, is used when computing standard errors with multi-way clustering (two or more clustering variables).{p_end}

{phang}{browse "https://ideas.repec.org/c/boc/bocode/s457527.html":parallel} by George Vega Yon is used when running reghdfe with multiple processors.{p_end}

{marker references}{...}
{title:References}

This program implements an extension on the fixed-point iteration proposed by:

{phang}
Paulo Guimaraes and Pedro Portugal. "A Simple Feasible Alternative Procedure to Estimate
Models with High-Dimensional Fixed Effects".
{it:Stata Journal, 10(4), 628-649, 2010.}
{browse "http://www.stata-journal.com/article.html?article=st0212":[link]}
{p_end}

A technical appendix of the algorithm employed is available in the PDF file below:

{phang}Sergio Correia. "Least Squares Iteration with Several High-Dimensional Fixed Effects". {it:Mimeo, 2014}
{browse "https://raw.githubusercontent.com/sergiocorreia/reghdfe/master/docs/explanation.pdf":[link]}
{p_end}

Note that the algorithm is not the same as the one used in:

{phang}
Torres, Sonia & Portugal, Pedro & Addison, John T. & Guimaraes, Paulo, 2013.
"The Sources of Wage Variation: A Three-Way High-Dimensional Fixed Effects Regression Model".
{it: IZA Discussion Papers 7276, Institute for the Study of Labor (IZA).}
{p_end}

{p 0 0 0}
If you use this program in your research, please cite either
the {browse "https://ideas.repec.org/c/boc/bocode/s457874.html":REPEC entry}
or the first article on this list.{p_end}

For details on the acceleration technique employed, please see "method 3" as described by:

{phang}
Macleod, Allan J. "Acceleration of vector sequences by multi-dimensional Delta-2 methods."
{it:Communications in Applied Numerical Methods 2.4 (1986): 385-392.}
{p_end}

For the rationale behind interacting fixed effects with continuous variables, see:

{phang}
Duflo, Esther. "The medium run effects of educational expansion: Evidence from a large school construction program in Indonesia."
{it:Journal of Development Economics 74.1 (2004): 163-197.}{browse "http://www.sciencedirect.com/science/article/pii/S0304387803001846": [link]}
{p_end}

Also see:

{phang}Abowd, J. M., R. H. Creecy, and F. Kramarz 2002.
Computing person and firm effects using linked longitudinal employer-employee data.
{it:Census Bureau Technical Paper TP-2002-06.}
{p_end}

{phang}
Cameron, A. Colin & Gelbach, Jonah B. & Miller, Douglas L., 2011.
"Robust Inference With Multiway Clustering,"
{it:Journal of Business & Economic Statistics, American Statistical Association, vol. 29(2), pages 238-249.}
{p_end}

{phang}
Gormley, T. & Matsa, D. 2014.
"Common errors: How to (and not to) control for unobserved heterogeneity."
{it:The Review of Financial Studies, vol. 27(2), pages 617-661.}
{p_end}

{phang}
Mittag, N. 2012.
"New methods to estimate models with large sets of fixed effects with an application to matched employer-employee data from Germany."
{it:{browse "http://doku.iab.de/fdz/reporte/2012/MR_01-12_EN.pdf":FDZ-Methodenreport 02/2012}.} 
{p_end}
