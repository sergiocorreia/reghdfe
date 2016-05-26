{smcl}
{* *! version 3.3.1 26may2016}{...}
{vieweralsosee "[R] areg" "help areg"}{...}
{vieweralsosee "[R] xtreg" "help xtreg"}{...}
{vieweralsosee "[R] ivregress" "help ivregress"}{...}
{vieweralsosee "" "--"}{...}
{vieweralsosee "reghdfe" "help reg2hdfe"}{...}
{vieweralsosee "a2reg" "help a2reg"}{...}
{viewerjumpto "Syntax" "hdfe##syntax"}{...}
{viewerjumpto "Contact" "hdfe##contact"}{...}
{title:Title}

{p2colset 5 18 20 2}{...}
{p2col :{cmd:hdfe} {hline 2}}Partial-out variables with respect to multiple levels of fixed-effects{p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{pstd}Replace current dataset:

{p 8 15 2}
{cmd:hdfe}
{help varlist}
{it:{weight}}
{cmd:,}
{opth a:bsorb(hdfe##absvar:absvars)}
{opt clear}
[{opth keepv:ars(varlist)} {opt keepid:s}]
[{opt clusterv:ars(varlist)} {help hdfe##options:options}]
{p_end}

{pstd}Keep current dataset and add new variables:

{p 8 15 2}
{cmd:hdfe}
{help varlist}
{it:{weight}}
{cmd:,}
{opth a:bsorb(hdfe##absvar:absvars)}
{opt g:enerate(stubname)}
[{opt sample:(newvarname)}]
[{opt clusterv:ars(varlist)} {help hdfe##options:options}]
{p_end}

{marker opt_summary}{...}
{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:HDFE-Specific}
{synopt :{opt clear:}}will overwrite the dataset; leaving the transformed variables, as well as some ancillary ones (such as the fixed effects, weights, cluster variables, etc.).{p_end}
{synopt:}If you use {cmd:hdfe} with factor variables, you may have trouble relating the old names (e.g. i.turn) to new names.{p_end}
{synopt:}The solution lies in this line: {cmd:mata: asarray(varlist_cache, "i.turn")}{p_end}
{synopt :{opth keepv:ars(varlist)}}keep additional variables{p_end}
{synopt :{opt keepid:s}}keep the temporary variables for the fixed effects (useful if you set them up like id#year){p_end}

{synopt :{opt g:enerate(stubname)}}will not overwrite the variables; instead creating new demeaned variables with the {it:stubname} prefix{p_end}
{synopt: {opt sample:(newvarname)}}will save the equivalent of e(sample) in this variable;
useful when dropping singletons. Used with the {opt g:enerate} option.{p_end}

{synopt: {opth clusterv:ars(varlist)}}list of variables containing cluster categories.
This is used to give more accurate number of degrees of freedom lost due to the fixed effects, as reported on r(df_a).{p_end}

{syntab:Diagnostic {help reghdfe##opt_diagnostic:[+]}}
{synopt :{opt v:erbose(#)}}amount of debugging information to show (0=None, 1=Some, 2=More, 3=Parsing/convergence details, 4=Every iteration){p_end}
{synopt :{opt time:it}}show elapsed times by stage of computation{p_end}

{syntab:Optimization {help reghdfe##opt_optimization:[+]}}
{p2coldent:+ {opth tol:erance(#)}}criterion for convergence (default=1e-8){p_end}
{synopt :{opth maxit:erations(#)}}maximum number of iterations (default=10,000); if set to missing ({cmd:.}) it will run for as long as it takes.{p_end}
{synopt :{opth pool:size(#)}}apply the within algorithm in groups of {it:#} variables (default 10). a large poolsize is usually faster but uses more memory{p_end}
{synopt :{opt accel:eration(str)}}acceleration method; options are conjugate_gradient (cg), steep_descent (sd), aitken (a), and none (no){p_end}
{synopt :{opt transf:orm(str)}}transform operation that defines the type of alternating projection; options are Kaczmarz (kac), Cimmino (cim), Symmetric Kaczmarz (sym){p_end}

{syntab:Degrees-of-Freedom Adjustments {help reghdfe##opt_dof:[+]}}
{synopt :{opt dof:adjustments(list)}}allows selecting the desired adjustments for degrees of freedom;
rarely used{p_end}
{synopt: {opth groupv:ar(newvar)}}unique identifier for the first mobility group{p_end}

{syntab:Reporting {help reghdfe##opt_reporting:[+]}}
{synopt :{opt version:}}reports the version number and date of hdfe, and saves it in e(version). standalone option{p_end}

{syntab:Undocumented}
{synopt :{opt keepsin:gletons}}do not drop singleton groups{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}* {opt absorb(absvars)} is required.{p_end}
{p 4 6 2}+ indicates a recommended or important option.{p_end}
{p 4 6 2}all variables may contain time-series operators and factor variables; see {help tsvarlist} and {help fvvarlist}.{p_end}
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
{p 4 6 2}Using categorical interactions (e.g. {it:x}{cmd:#}{it:z}) is faster than running {it:egen group(...)} beforehand.{p_end}
{p 4 6 2}Singleton obs. are dropped iteratively until no more singletons are found (see ancilliary article for details).{p_end}
{p 4 6 2}Slope-only absvars ("state#c.time") have poor numerical stability and slow convergence.
If you need those, either i) increase tolerance or
ii) use slope-and-intercept absvars ("state##c.time"), even if the intercept is redundant.
For instance if absvar is "i.zipcode i.state##c.time" then i.state is redundant given i.zipcode, but
convergence will still be {it:much} faster.{p_end}


{marker description}{...}
{title:Description}

{pstd}{cmd:hdfe} computes the residuals of a set of variables with respect to multiple levels of fixed effects.
It is a generalization of the {it:within} transformation done by {help areg} and {help xtreg: xtreg,fe} for more than one fixed effect, also allowing for multiple heterogeneous intercepts.{p_end}

{pstd}{cmd:hdfe} is a programmers' routine that serves as a building block to other regression packages so they can support multiple fixed effects
(see for instance {search binscatter}, {browse "https://github.com/matthieugomez/stata-regife":regife} and {search poi2hdfe}).
It contains the same code underlying {help reghdfe} and exposes most of its functionality and options.

{pstd}It also computes the degrees-of-freedom absorbed by the fixed effects and stores them in e(df_a).{p_end}

{pstd}It works well with other building-block packages such as {help avar} (from SSC).{p_end}


{marker examples}{...}
{title:Example Usage}

{pstd}Suppose you want to replicate {cmd:reghdfe}. Then, you would do:{p_end}

{phang2}{cmd:. sysuse auto, clear}{p_end}

{phang2}{cmd:. * Benchmark}{p_end}
{phang2}{cmd:. reghdfe price weight length, a(turn trunk)}{p_end}

{phang2}{cmd:. * Demean variables}{p_end}
{phang2}{cmd:. hdfe price weight length, a(turn trunk) gen(RESID_)}{p_end}
{phang2}{cmd:. local df_a = e(df_a)}{p_end}

{phang2}{cmd:. * Run regression}{p_end}
{phang2}{cmd:. quietly regress RESID_*, nocons}{p_end}

{phang2}{cmd:. * Fix degrees-of-freedom}{p_end}
{phang2}{cmd:. local df_r = e(df_r) - `df_a'}{p_end}
{phang2}{cmd:. matrix b = e(b)}{p_end}
{phang2}{cmd:. matrix V = e(V) * e(df_r) / `df_r'}{p_end}
{phang2}{cmd:. ereturn post b V, dep(price) obs(`c(N)') dof(`df_r')}{p_end}
{phang2}{cmd:. ereturn display}{p_end}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:hdfe} stores the following in {cmd:e()}:

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars}{p_end}
{synopt:{cmd:e(df_a)}}degrees of freedom lost due to the fixed effects (taking into account the cluster structure and whether the FEs are nested within the clusters){p_end}
{synopt:{cmd:e(N_hdfe)}}number of sets of fixed effects{p_end}

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Macros}{p_end}
{synopt:{cmd:e(absvars)}}canonical expansion of the fixed effects{p_end}
{synopt:{cmd:e(extended_absvars)}}expansion of the fixed effects separating heterogeneous slopes (e.g. {it:y##c.z} is expanded to {it:y y#c.z}){p_end}

{marker contact}{...}
{title:Author}

{pstd}Sergio Correia{break}
Fuqua School of Business, Duke University{break}
Email: {browse "mailto:sergio.correia@duke.edu":sergio.correia@duke.edu}
{p_end}

{marker user_guide}{...}
{title:User Guide}

{pstd}
A copy of this help file, as well as a more in-depth user guide is in development and will be available at {browse "http://scorreia.com/reghdfe"}.{p_end}

{marker updates}{...}
{title:Latest Updates}

{pstd}
{cmd:hdfe} is updated frequently, and upgrades or minor bug fixes may not be immediately available in SSC.
To check or contribute to the latest version of hdfe, explore the
{browse "https://github.com/sergiocorreia/reghdfe":Github repository}.
Bugs or missing features can be discussed through email or at the {browse "https://github.com/sergiocorreia/reghdfe/issues":Github issue tracker}.{p_end}

{pstd}
To see your current version and installed dependencies, type {cmd:reghdfe, version}
{p_end}

{marker acknowledgements}{...}
{title:Acknowledgements}

{pstd}
This package wouldn't have existed without the invaluable feedback and contributions of Paulo Guimaraes,  Amine Ouazad, Mark Schaffer and Kit Baum. Also invaluable are the great bug-spotting abilities of many users.{p_end}

