{smcl}
{* *! version 3.2.9 21feb2016}{...}
{vieweralsosee "[R] areg" "help areg"}{...}
{vieweralsosee "[R] xtreg" "help xtreg"}{...}
{vieweralsosee "[R] ivregress" "help ivregress"}{...}
{vieweralsosee "" "--"}{...}
{vieweralsosee "ivreg2" "help ivreg2"}{...}
{vieweralsosee "ivregress" "help ivregress"}{...}
{vieweralsosee "reg2hdfe" "help reg2hdfe"}{...}
{vieweralsosee "a2reg" "help a2reg"}{...}
{viewerjumpto "Syntax" "reghdfe3##syntax"}{...}
{viewerjumpto "Description" "reghdfe3##description"}{...}
{viewerjumpto "Options" "reghdfe3##options"}{...}
{viewerjumpto "Postestimation Syntax" "reghdfe3##postestimation"}{...}
{viewerjumpto "Remarks" "reghdfe3##remarks"}{...}
{viewerjumpto "Examples" "reghdfe3##examples"}{...}
{viewerjumpto "Stored results" "reghdfe3##results"}{...}
{viewerjumpto "Author" "reghdfe3##contact"}{...}
{viewerjumpto "Updates" "reghdfe3##updates"}{...}
{viewerjumpto "Acknowledgements" "reghdfe3##acknowledgements"}{...}
{viewerjumpto "References" "reghdfe3##references"}{...}
{title:Title}

{p2colset 5 18 20 2}{...}
{p2col :{cmd:reghdfe} {hline 2}}Linear and instrumental-variable/GMM regression absorbing multiple levels of fixed effects{p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{p 8 15 2} {cmd:reghdfe}
{depvar} [{indepvars}] [{cmd:(}{it:{help varlist:endogvars}} {cmd:=} {it:{help varlist:iv_vars}}{cmd:)}]
{ifin} {it:{weight}} {cmd:,} {opth a:bsorb(reghdfe3##absvar:absvars)} [{help reghdfe3##options:options}] {p_end}

{marker opt_summary}{...}
{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Model {help reghdfe3##opt_model:[+]}}
{p2coldent:* {opth a:bsorb(reghdfe3##absvar:absvars)}}identifiers of the absorbed fixed effects; each {help reghdfe3##absvar:absvar} represents one set of fixed effects{p_end}
{synopt: {cmdab:a:bsorb(}{it:...}{cmd:,} {cmdab:save:fe)}}save all fixed effect estimates ({it:__hdfe*} prefix); useful for a subsequent {help reghdfe3##postestimation:predict}.
However, see also the {it:resid} option.{p_end}
{synopt : {opth res:iduals(newvar)}}save residuals; more direct and much faster than saving the fixed effects and then running predict{p_end}
{synopt :{opth su:mmarize(tabstat##statname:stats)}}equivalent to {help reghdfe3##postestimation:estat summarize} after the regression,
but more flexible, compatible with the {opt fast:} option, and saves results on {it:e(summarize)}{p_end}
{synopt : {opt subopt:ions(...)}}additional options that will be passed to the regression command (either {help regress}, {help ivreg2}, or {help ivregress}){p_end}
	
{syntab:SE/Robust {help reghdfe3##opt_vce:[+]}}
{p2coldent:+ {opt vce}{cmd:(}{help reghdfe3##opt_vce:vcetype} [{cmd:,}{it:opt}]{cmd:)}}{it:vcetype}
may be {opt un:adjusted} (default), {opt r:obust} or {opt cl:uster} {help fvvarlist} (allowing two- and multi-way clustering){p_end}
{synopt :}suboptions {opt bw(#)}, {opt ker:nel(str)}, {opt dkraay(#)} and {opt kiefer} allow for AC/HAC estimates; see the {help avar} package{p_end}

{syntab:Instrumental-Variable/2SLS/GMM {help reghdfe3##opt_iv:[+]}}
{synopt :{opt est:imator(str)}}either {opt 2sls} (default), {opt gmm:2s} (two-stage GMM),
{opt liml} (limited-information maximum likelihood) or {opt cue} (which gives approximate results, see discussion below){p_end}
{synopt :{opt stage:s(list)}}estimate additional regressions; choose any of {opt first} {opt ols} {opt reduced} {opt acid} (or {opt all}){p_end}
{synopt :{opt ff:irst}}compute first-stage diagnostic and identification statistics{p_end}
{synopt :{opth iv:suite(subcmd)}}package used in the IV/GMM regressions;
options are {opt ivreg2} (default; needs installing) and {opt ivregress}{p_end}

{syntab:Diagnostic {help reghdfe3##opt_diagnostic:[+]}}
{synopt :{opt v:erbose(#)}}amount of debugging information to show (0=None, 1=Some, 2=More, 3=Parsing/convergence details, 4=Every iteration){p_end}
{synopt :{opt time:it}}show elapsed times by stage of computation{p_end}

{syntab:Optimization {help reghdfe3##opt_optimization:[+]}}
{p2coldent:+ {opth tol:erance(#)}}criterion for convergence (default=1e-8){p_end}
{synopt :{opth maxit:erations(#)}}maximum number of iterations (default=10,000); if set to missing ({cmd:.}) it will run for as long as it takes.{p_end}
{synopt :{opth pool:size(#)}}apply the within algorithm in groups of {it:#} variables (default 10). a large poolsize is usually faster but uses more memory{p_end}
{synopt :{opt accel:eration(str)}}acceleration method; options are conjugate_gradient (cg), steep_descent (sd), aitken (a), and none (no){p_end}
{synopt :{opt transf:orm(str)}}transform operation that defines the type of alternating projection; options are Kaczmarz (kac), Cimmino (cim), Symmetric Kaczmarz (sym){p_end}

{syntab:Speedup Tricks {help reghdfe3##opt_speedup:[+]}}
{synopt :{cmd: cache(save} [,opt]{cmd:)}}absorb all variables without regressing (destructive; combine it with {help preserve:preserve/restore}){p_end}
{synopt :}suboption {opth keep(varlist)} adds additional untransformed variables to the resulting dataset{p_end}
{synopt :{cmd: cache(use)}}run regressions on cached data; {it:vce()} must be the same as with {cmd: cache(save)}.{p_end}
{synopt :{cmd: cache(clear)}}delete Mata objects to clear up memory; no more regressions can be run after this{p_end}
{synopt :{opt fast}}will not create {it:e(sample)}; disabled when saving fixed effects, residuals or mobility groups{p_end}

{syntab:Degrees-of-Freedom Adjustments {help reghdfe3##opt_dof:[+]}}
{synopt :{opt dof:adjustments(list)}}allows selecting the desired adjustments for degrees of freedom;
rarely used{p_end}
{synopt: {opth groupv:ar(newvar)}}unique identifier for the first mobility group{p_end}

{syntab:Reporting {help reghdfe3##opt_reporting:[+]}}
{synopt :{opt version:}}reports the version number and date of reghdfe, and saves it in e(version). standalone option{p_end}
{synopt :{opt l:evel(#)}}set confidence level; default is {cmd:level(95)}{p_end}
{synopt :{it:{help reghdfe3##display_options:display_options}}}control column formats, row spacing, line width, display of omitted variables and base and empty cells, and factor-variable labeling.{p_end}
{synopt :}particularly useful are the {opt noomit:ted} and {opt noempty} options to hide regressors omitted due to collinearity{p_end}

{syntab:Undocumented}
{synopt :{opt keepsin:gletons}}do not drop singleton groups{p_end}
{synopt :{opt old}}will call the latest 2.x version of reghdfe instead (see the {help reghdfe3:old help file}){p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}* {opt absorb(absvars)} is required.{p_end}
{p 4 6 2}+ indicates a recommended or important option.{p_end}
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
{synopt:{cmd:i.}{it:var1}{cmd:#}{cmd:c.}{it:var2}}absorb heterogeneous slopes, where {it:var2} has a different slope coef. depending on the category of {it:var1}{p_end}
{synopt:{it:var1}{cmd:##}{cmd:c.}{it:var2}}equivalent to "{cmd:i.}{it:var1} {cmd:i.}{it:var1}{cmd:#}{cmd:c.}{it:var2}", but {it:much} faster{p_end}
{synopt:{it:var1}{cmd:##c.(}{it:var2 var3}{cmd:)}}multiple heterogeneous slopes are allowed together. Alternative syntax: {it:var1}{cmd:##(c.}{it:var2} {cmd:c.}{it:var3}{cmd:)}{p_end}
{synopt:{it:v1}{cmd:#}{it:v2}{cmd:#}{it:v3}{cmd:##c.(}{it:v4 v5}{cmd:)}}factor operators can be combined{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}To save the estimates specific absvars, write {newvar}{inp:={it:absvar}}.{p_end}
{p 4 6 2}Please be aware that in most cases these estimates are neither consistent nor econometrically identified.{p_end}
{p 4 6 2}Using categorical interactions (e.g. {it:x}{cmd:#}{it:z}) is faster than running {it:egen group(...)} beforehand.{p_end}
{p 4 6 2}Singleton obs. are dropped iteratively until no more singletons are found (see ancilliary article for details).{p_end}
{p 4 6 2}Slope-only absvars ("state#c.time") have poor numerical stability and slow convergence.
If you need those, either i) increase tolerance or
ii) use slope-and-intercept absvars ("state##c.time"), even if the intercept is redundant.
For instance if absvar is "i.zipcode i.state##c.time" then i.state is redundant given i.zipcode, but
convergence will still be {it:much} faster.{p_end}

{marker description}{...}
{title:Description}

{pstd}
{cmd:reghdfe} is a generalization of {help areg} (and {help xtreg:xtreg,fe}, {help xtivreg:xtivreg,fe}) for multiple levels of fixed effects
(including heterogeneous slopes), alternative estimators (2sls, gmm2s, liml), and additional robust standard errors (multi-way clustering, HAC standard errors, etc).{p_end}

{pstd}Additional features include:{p_end}

{p2col 8 12 12 2: a)}A novel and robust algorithm to efficiently absorb the fixed effects (extending the work of Guimaraes and Portugal, 2010).{p_end}
{p2col 8 12 12 2: b)}Coded in Mata, which in most scenarios makes it even faster than {it:areg} and {it:xtreg} for a single fixed effect (see benchmarks on the Github page).{p_end}
{p2col 8 12 12 2: c)}Can save the point estimates of the fixed effects ({it:caveat emptor}: the fixed effects may not be identified, see the {help reghdfe3##references:references}).{p_end}
{p2col 8 12 12 2: d)}Calculates the degrees-of-freedom lost due to the fixed effects
(note: beyond two levels of fixed effects, this is still an open problem, but we provide a conservative approximation).{p_end}
{p2col 8 12 12 2: e)}Iteratively removes singleton groups by default, to avoid biasing the standard errors (see ancillary document).{p_end}

{marker options}{...}
{title:Options}

{marker opt_model}{...}
{dlgtab:Model and Miscellanea}

{phang}
{opth a:bsorb(reghdfe3##absvar:absvars)} list of categorical variables (or interactions) representing the fixed effects to be absorbed.
this is equivalent to including an indicator/dummy variable for each category of each {it:absvar}. {cmd:absorb()} is required.

{pmore}
To save a fixed effect, prefix the absvar with "{newvar}{cmd:=}".
For instance, the option {cmd:absorb(firm_id worker_id year_coefs=year_id)} will include firm,
worker and year fixed effects, but will only save the estimates for the year fixed effects (in the new variable {it:year_coefs}).

{pmore}
If you want to {help reghdfe3##postestimation:predict} afterwards but don't care about setting the names of each fixed effect, use the {cmdab:save:fe} suboption.
This will delete all variables named {it:__hdfe*__} and create new ones as required.
Example: {it:reghdfe price weight, absorb(turn trunk, savefe)}

{phang}
{opth res:iduals(newvar)} will save the regression residuals in a new variable.

{pmore}
This is a superior alternative than running {cmd:predict, resid} afterwards as it's faster and doesn't require saving the fixed effects.

{phang}
{opth su:mmarize(tabstat##statname:stats)} will report and save a table of summary of statistics of the regression
variables (including the instruments, if applicable), using the same sample as the regression.

{pmore} {opt su:mmarize} (without parenthesis) saves the default set of statistics: {it:mean min max}.

{pmore} The complete list of accepted statistics is available in the {help tabstat##statname:tabstat help}. The most useful are {it:count range sd median p##}.

{pmore} The summary table is saved in {it:e(summarize)}

{pmore} To save the summary table silently (without showing it after the regression table), use the {opt qui:etly} suboption. You can use it by itself ({cmd:summarize(,quietly)}) or with custom statistics ({cmd:summarize(mean, quietly)}).

{phang}
{opt subopt:ions(...)}
options that will be passed directly to the regression command (either {help regress}, {help ivreg2}, or {help ivregress})

{marker opt_vce}{...}
{dlgtab:SE/Robust}

{phang}
{opth vce:(reghdfe3##vcetype:vcetype, subopt)}
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
(e.g. allowing for intragroup correlation across individuals, time, country, etc).

{pmore}
Each {it:clustervar} permits interactions of the type {it:var1{cmd:#}var2}
(this is faster than using {cmd:egen group()} for a one-off regression).

{pmore} Warning: The number of clusters, for all of the cluster variables, must go off to infinity.
A frequent rule of thumb is that each cluster variable must have at least 50 different categories
(the number of categories for each clustervar appears on the header of the regression table).

{pstd}
The following suboptions require either the {help ivreg2} or the {help avar} package from SSC.
For a careful explanation, see the {help ivreg2##s_robust:ivreg2 help file}, from which the comments below borrow.

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
robust standard errors, and a gmm2s estimator, reghdfe will translate
{opt vce(robust)} into {opt wmatrix(robust)} {opt vce(unadjusted)}.
This maintains compatibility with {cmd:ivreg2} and other packages, but may unadvisable as described in {help ivregress} (technical note). Specifying this option will instead use {opt wmatrix(robust)} {opt vce(robust)}.

{pmore}However, computing the second-step vce matrix requires computing updated estimates (including updated fixed effects).
Since reghdfe currently does not allow this, the resulting standard errors
{hi:will not be exactly the same as with ivregress}.
This issue is similar to applying the CUE estimator, described further below.

{pmore}Note: The above comments are also appliable to clustered standard error.

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
{opt stage:s(list)}
adds and saves up to four auxiliary regressions useful when running instrumental-variable regressions:

{phang2}{cmd:first} all first-stage regressions{p_end}
{phang2}{cmd:ols} ols regression (between dependent variable and endogenous variables; useful as a benchmark){p_end}
{phang2}{cmd:reduced} reduced-form regression (ols regression with included and excluded instruments as regressors){p_end}
{phang2}{cmd:acid} an "acid" regression that includes both instruments and endogenous variables as regressors; in this setup, excluded instruments should not be significant.{p_end}

{pmore}
You can pass suboptions not just to the iv command but to all stage regressions with a comma after the list of stages. Example:{break}
{cmd:reghdfe price (weight=length), absorb(turn) subopt(nocollin) stages(first, eform(exp(beta)) )}

{pmore}
By default all stages are saved (see {help estimates dir}).
The suboption {cmd:,nosave} will prevent that.
However, future {cmd:replay}s will only replay the iv regression.

{phang}
{opt ffirst}
compute and report first stage statistics ({help ivreg2##s_relevance:details}); requires the ivreg2 package.

{pmore}
These statistics will be saved on the {it:e(first)} matrix. 
If the first-stage estimates are also saved (with the {cmd:stages()} option), the respective statistics will be copied to {cmd:e(first_*)}.

{phang}
{opth iv:suite(subcmd)}
allows the IV/2SLS regression to be run either using {opt ivregress} or {opt ivreg2}.

{pmore} {opt ivreg2} is the default, but needs to be installed for that option to work.

{marker opt_diagnostic}{...}
{dlgtab:Diagnostic}

{phang}
{opt v:erbose(#)} orders the command to print debugging information.

{pmore}
Possible values are 0 (none), 1 (some information), 2 (even more), 3 (adds dots for each iteration, and reportes parsing details), 4 (adds details for every iteration step)

{pmore}
For debugging, the most useful value is 3. For simple status reports, set verbose to 1.

{phang}
{opt time:it} shows the elapsed time at different steps of the estimation. Most time is usually spent on three steps: map_precompute(), map_solve() and the regression step.
 
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
{opth groupv:ar(newvar)} name of the new variable that will contain the first mobility group.
Requires {opt pair:wise}, {opt first:pair}, or the default {opt all}.

{marker opt_speedup}{...}
{dlgtab:Speeding Up Estimation}

{phang}
{cmd:reghdfe} {varlist} {ifin}{cmd:,} {opt a:bsorb(absvars)} {cmd:save(cache)} [{it:options}]

{pmore}
This will transform {it:varlist}, absorbing the fixed effects indicated by {it:absvars}.
It is useful when running a series of alternative specifications with common variables, as the variables will only be transformed once instead of every time a regression is run.

{pmore}
It replaces the current dataset, so it is a good idea to precede it with a {help preserve} command

{pmore}
To keep additional (untransformed) variables in the new dataset, use the {opth keep(varlist)} suboption.

{phang}
{cmd:cache(use)} is used when running reghdfe after a {it:save(cache)} operation. Both the {it:absorb()} and {it:vce()} options must be the same as when the cache was created (the latter because the degrees of freedom were computed at that point).

{phang}
{cmd:cache(clear)} will delete the Mata objects created by {it:reghdfe} and kept in memory after the {it:save(cache)} operation. These objects may consume a lot of memory, so it is a good idea to clean up the cache. Additionally, if you previously specified {it:preserve}, it may be a good time to {it:restore}.

{pmore}Example:{p_end}
{phang2}{cmd:. sysuse auto}{p_end}
{phang2}{cmd:. preserve}{p_end}
{phang2}{cmd:.}{p_end}
{phang2}{cmd:. * Save the cache}{p_end}
{phang2}{cmd:. reghdfe price weight length, a(turn rep) vce(turn) cache(save, keep(foreign))}{p_end}
{phang2}{cmd:.}{p_end}
{phang2}{cmd:. * Run regressions}{p_end}
{phang2}{cmd:. reghdfe price weight, a(turn rep) cache(use)}{p_end}
{phang2}{cmd:. reghdfe price length, a(turn rep) cache(use)}{p_end}
{phang2}{cmd:.}{p_end}
{phang2}{cmd:. * Clean up}{p_end}
{phang2}{cmd:. reghdfe, cache(clear)}{p_end}
{phang2}{cmd:. restore}{p_end}

{phang}
{opt fast} avoids saving {it:e(sample)} into the regression.
Since saving the variable only involves copying a Mata vector, the speedup is currently quite small.
Future versions of reghdfe may change this as features are added.

{pmore}
Note that {opt fast} will be disabled when adding variables to the dataset (i.e. when saving residuals, fixed effects, or mobility groups), and is incompatible with most postestimation commands.

{pmore}
If you wish to use {opt fast} while reporting {cmd:estat summarize}, see the {opt summarize} option.

{marker opt_optimization}{...}
{dlgtab:Optimization}

{phang}
{opth tol:erance(#)} specifies the tolerance criterion for convergence; default is {cmd:tolerance(1e-8)}

{pmore}
Note that for tolerances beyond 1e-14, the limits of the {it:double} precision are reached and the results will most likely not converge.

{pmore}
At the other end, is not tight enough, the regression may not identify perfectly collinear regressors. However, those cases can be easily spotted due to their extremely high standard errors.

{pmore}
Warning: when absorbing heterogeneous slopes without the accompanying heterogeneous intercepts, convergence is quite poor and a tight tolerance is strongly suggested (i.e. higher than the default). In other words, an absvar of {it:var1##c.var2} converges easily, but an absvar of {it:var1#c.var2} will converge slowly and may require a tighter tolerance.

{phang}
{opth maxit:erations(#)}
specifies the maximum number of iterations; the default is {cmd:maxiterations(10000)}; set it to missing ({cmd:.}) to run forever until convergence.

{phang}
{opth pool:size(#)}
Number of variables that are {it:pooled together} into a matrix that will then be transformed.
The default is to pool variables in groups of 5. Larger groups are faster with more than one processor, but may cause out-of-memory errors. In that case, set poolsize to 1.

{phang}
{it:Advanced options:}

{phang}
{opt acceleration(str)} allows for different acceleration techniques, from the simplest case of 
no acceleration ({opt no:ne}), to steep descent ({opt st:eep_descent} or {opt sd}), Aitken ({opt a:itken}),
and finally Conjugate Gradient ({opt co:njugate_gradient} or {opt cg}).

{pmore}
Note: Each acceleration is just a plug-in Mata function, so a larger number of acceleration techniques are available, albeit undocumented (and slower).

{phang}
{opt transf:orm(str)} allows for different "alternating projection" transforms. The classical transform is Kaczmarz ({opt kac:zmarz}), and more stable alternatives are Cimmino ({opt cim:mino}) and Symmetric Kaczmarz ({opt sym:metric_kaczmarz})

{pmore}
Note: Each transform is just a plug-in Mata function, so a larger number of acceleration techniques are available, albeit undocumented (and slower).

{pmore}
Note: The default acceleration is Conjugate Gradient and the default transform is Symmetric Kaczmarz. Be wary that different accelerations often work better with certain transforms. For instance, do not use conjugate gradient with plain Kaczmarz, as it will not converge.

{phang}
{opt precondition} {it:(currently disabled)}

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
{p_end}{col 23}May require you to previously save the fixed effects (except for option {opt xb}).
{col 23}To see how, see the details of the {help reghdfe3##absvar:absorb} option
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


{col 8}{cmd:test}{col 23}Performs significance test on the parameters, see the {help test:stata help}

{col 8}{cmd:suest}{col 23}Do not use {cmd:suest}. It will run, but the results will be incorrect. See workaround below

{pmore}If you want to perform tests that are usually run with {cmd:suest},
such as non-nested models, tests using alternative specifications of the variables,
or tests on different groups, you can replicate it manually, as described 
{browse "http://www.stata.com/statalist/archive/2009-11/msg01485.html":here}.
{p_end}

{marker remarks}{...}

{title:Possible Pitfalls and Common Mistakes}

{p2col 8 12 12 2: 1.}(note: as of version 2.1, the constant is no longer reported) Ignore the constant; it doesn't tell you much. If you want to use descriptive stats, that's what the {opt sum:marize()} and {cmd:estat summ} commands are for.
Even better, use {opt noconstant} to drop it (although it's not really dropped as it never existed on the first place!){p_end}
{p2col 8 12 12 2: 2.}Think twice before saving the fixed effects. They are probably inconsistent / not identified and you will likely be using them wrong.{p_end}
{p2col 8 12 12 2: 3.}(note: as of version 3.0 singletons are dropped by default) It's good practice to drop singletons. {opt dropsi:ngleton} is your friend.{p_end}
{p2col 8 12 12 2: 4.}If you use {opt vce(robust)}, be sure that your {it:other} dimension is not "fixed" but grows with N, or your SEs will be wrong.{p_end}
{p2col 8 12 12 2: 5.}If you use {opt vce(cluster ...)}, check that your number of clusters is high enough (50+ is a rule of thumb). If not, you are making the SEs even worse!{p_end}
{p2col 8 12 12 2: 6.}The panel variables (absvars) should probably be nested within the clusters (clustervars) due to the within-panel correlation induced by the FEs.
(this is not the case for *all* the absvars, only those that are treated as growing as N grows){p_end}
{p2col 8 12 12 2: 7.}If you run analytic or probability weights,
you are responsible for ensuring that the weights stay
constant within each unit of a fixed effect (e.g. individual),
or that it is correct to allow varying-weights for that case.
{p_end}
{p2col 8 12 12 2: 8.}Be aware that adding several HDFEs is not a panacea.
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

{phang}(If you are interested in discussing these or others, feel free to {help reghdfe3##contact:contact me})

{phang}Code, medium term:

{p2col 8 12 12 2: -}Complete GT preconditioning (v4){p_end}
{p2col 8 12 12 2: -}Improve algorithm that recovers the fixed effects (v5){p_end}
{p2col 8 12 12 2: -}Improve statistics and tests related to the fixed effects (v5){p_end}
{p2col 8 12 12 2: -}Implement a -bootstrap- option in DoF estimation (v5){p_end}

{phang}Code, long term:

{p2col 8 12 12 2: -}The interaction with cont vars (i.a#c.b) may suffer from numerical accuracy issues, as we are dividing by a sum of squares{p_end}
{p2col 8 12 12 2: -}Calculate exact DoF adjustment for 3+ HDFEs (note: not a problem with cluster VCE when one FE is nested within the cluster){p_end}
{p2col 8 12 12 2: -}More postestimation commands (lincom? margins?){p_end}

{phang}Theory:

{p2col 8 12 12 2: -}Add a more thorough discussion on the possible identification issues{p_end}
{p2col 8 12 12 2: -}Find out a way to use reghdfe iteratively with CUE
(right now only OLS/2SLS/GMM2S/LIML give the exact same results){p_end}
{p2col 8 12 12 2: -}Not sure if I should add an F-test for the absvars in the vce(robust) and vce(cluster) cases.
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

{pstd}
{it:Note: it also keeps most e() results placed by the regression subcommands (ivreg2, ivregress)}

{synoptset 24 tabbed}{...}
{syntab:Scalars}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(N_hdfe)}}number of absorbed fixed-effects{p_end}
{synopt:{cmd:e(tss)}}total sum of squares{p_end}
{synopt:{cmd:e(rss)}}residual sum of squares{p_end}
{synopt:{cmd:e(r2)}}R-squared{p_end}
{synopt:{cmd:e(r2_a)}}adjusted R-squared{p_end}
{synopt:{cmd:e(r2_within)}}Within R-squared{p_end}
{synopt:{cmd:e(r2_a_within)}}Adjusted Within R-squared{p_end}
{synopt:{cmd:e(df_a)}}degrees of freedom lost due to the fixed effects{p_end}
{synopt:{cmd:e(rmse)}}root mean squared error{p_end}
{synopt:{cmd:e(ll)}}log-likelihood{p_end}
{synopt:{cmd:e(ll_0)}}log-likelihood of fixed-effect-only regression{p_end}
{synopt:{cmd:e(F)}}F statistic{p_end}
{synopt:{cmd:e(F_absorb)}}F statistic for absorbed effect {it:note: currently disabled}{p_end}
{synopt:{cmd:e(rank)}}rank of {cmd:e(V)}{p_end}
{synopt:{cmd:e(N_clustervars)}}number of cluster variables{p_end}
        
{synopt:{cmd:e(clust}#{cmd:)}}number of clusters for the #th cluster variable{p_end}
{synopt:{cmd:e(N_clust)}}number of clusters; minimum of {it:e(clust#)}{p_end}

{synopt:{cmd:e(K}#{cmd:)}}Number of categories of the #th absorbed FE{p_end}
{synopt:{cmd:e(M}#{cmd:)}}Number of redundant categories of the #th absorbed FE{p_end}
{synopt:{cmd:e(mobility)}}Sum of all {cmd:e(M#)}{p_end}
{synopt:{cmd:e(df_m)}}model degrees of freedom{p_end}
{synopt:{cmd:e(df_r)}}residual degrees of freedom{p_end}

{synoptset 24 tabbed}{...}
{syntab:Macros}
{synopt:{cmd:e(cmd)}}{cmd:reghdfe}{p_end}
{synopt:{cmd:e(subcmd)}}either {cmd:regress}, {cmd:ivreg2} or {cmd:ivregress}{p_end}
{synopt:{cmd:e(model)}}{cmd:ols}, {cmd:iv}, {cmd:gmm2s}, {cmd:liml} or {cmd:cue}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(dofmethod)}}dofmethod employed in the regression{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(indepvars)}}names of independent variables{p_end}
{synopt:{cmd:e(endogvars)}}names of endogenous right-hand-side variables{p_end}
{synopt:{cmd:e(instruments)}}names of excluded instruments{p_end}
{synopt:{cmd:e(absvars)}}name of the absorbed variables or interactions{p_end}
{synopt:{cmd:e(title)}}title in estimation output{p_end}
{synopt:{cmd:e(clustvar)}}name of cluster variable{p_end}
{synopt:{cmd:e(clustvar}#{cmd:)}}name of the #th cluster variable{p_end}
{synopt:{cmd:e(vce)}}{it:vcetype} specified in {cmd:vce()}{p_end}
{synopt:{cmd:e(vcetype)}}title used to label Std. Err.{p_end}
{synopt:{cmd:e(stage)}}stage within an IV-regression; only if {it:stages()} was used{p_end}
{synopt:{cmd:e(properties)}}{cmd:b V}{p_end}

{synoptset 24 tabbed}{...}
{syntab:Matrices}
{synopt:{cmd:e(b)}}coefficient vector{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix of the estimators{p_end}

{synoptset 24 tabbed}{...}
{syntab:Functions}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}
{p2colreset}{...}

{marker contact}{...}
{title:Author}

{pstd}Sergio Correia{break}
Fuqua School of Business, Duke University{break}
Email: {browse "mailto:sergio.correia@gmail.com":sergio.correia@gmail.com}
{p_end}

{marker user_guide}{...}
{title:User Guide}

{pstd}
A copy of this help file, as well as a more in-depth user guide is in development and will be available at {browse "http://scorreia.com/reghdfe"}.{p_end}

{marker updates}{...}
{title:Latest Updates}

{pstd}
{cmd:reghdfe} is updated frequently, and upgrades or minor bug fixes may not be immediately available in SSC.
To check or contribute to the latest version of reghdfe, explore the
{browse "https://github.com/sergiocorreia/reghdfe":Github repository}.
Bugs or missing features can be discussed through email or at the {browse "https://github.com/sergiocorreia/reghdfe/issues":Github issue tracker}.{p_end}

{pstd}
To see your current version and installed dependencies, type {cmd:reghdfe, version}
{p_end}

{marker acknowledgements}{...}
{title:Acknowledgements}

{pstd}
This package wouldn't have existed without the invaluable feedback and contributions of Paulo Guimaraes,  Amine Ouazad, Mark Schaffer and Kit Baum. Also invaluable are the great bug-spotting abilities of many users.{p_end}

{pstd}In addition, {it:reghdfe} is build upon important contributions from the Stata community:{p_end}

{phang}{browse "https://ideas.repec.org/c/boc/bocode/s457101.html":reg2hdfe}, from Paulo Guimaraes,
and {browse "https://ideas.repec.org/c/boc/bocode/s456942.html":a2reg} from Amine Ouazad,
 were the inspiration and building blocks on which reghdfe was built.{p_end}
 
{phang}{browse "http://www.repec.org/bocode/i/ivreg2.html":ivreg2}, by Christopher F Baum, Mark E Schaffer and Steven Stillman, is the package used by default for instrumental-variable regression.{p_end}

{phang}{browse "https://ideas.repec.org/c/boc/bocode/s457689.html":avar} by Christopher F Baum and Mark E Schaffer, is the package used for estimating the HAC-robust standard errors of ols regressions.{p_end}

{phang}{browse "http://econpapers.repec.org/software/bocbocode/s456797.htm":tuples} by Joseph Lunchman and Nicholas Cox, is used when computing standard errors with multi-way clustering (two or more clustering variables).{p_end}

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
It addresses many of the limitation of previous works, such as possible lack of convergence, arbitrary slow convergence times,
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

{p 0 0 0}
For the rationale behind interacting fixed effects with continuous variables, see:

{phang}
Duflo, Esther. "The medium run effects of educational expansion: Evidence from a large school construction program in Indonesia."
{it:Journal of Development Economics 74.1 (2004): 163-197.}{browse "http://www.sciencedirect.com/science/article/pii/S0304387803001846": [link]}
{p_end}

{p 0 0 0}
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
