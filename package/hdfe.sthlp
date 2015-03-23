{smcl}
{* *! version 1.0.0  12mar2015}{...}
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
{p2col :{cmd:hdfe} {hline 2}}Partial-out variables with respect to a series of fixed-effects{p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{p 8 15 2}
{cmd:hdfe}
{help varlist}
{it:{weight}}
{cmd:,}
{opth a:bsorb(hdfe##absvar:absvars)}
[{opt g:enerate(stubname)} | clear]
[{opth cluster:vars(varlist)}
{opt dropsi:ngletons}
{opt sample:(newvarname)}
{opth cores(#)}
{opt v:erbose(#)}
{opt tol:erance(#)}
{opt maxi:terations(#)}]
{it: maximize_options}]
{p_end}

{pstd}{it: Notes:}{p_end}
{p 5 7 2}
- this is a programmers' command, like -avar-. For a detailed explanation and comments, see the help and website for the reghdfe package.
{p_end}
{p 5 7 2}
- does not accept time series or factor variables
{p_end}
{p 5 7 2}
- varlist and clustervars MUST BE FULLY spelled out (i.e. you need to use unab beforehand!), but that is not needed at all for the absvars.
{p_end}

{marker absvar}{...}
{it:Absvars:}

{synoptset 22}{...}
{synopthdr:absvar}
{synoptline}
{synopt:{cmd:i.}{it:varname}}indicators for each level of {it:varname} (the {cmd:i.} prefix is tacit and can be omitted).{p_end}
{synopt:{it:var1}{cmd:#}{it:var2}}indicators for each combination of levels of {it:var1} and {it:var2} (same as {cmd:i.}{it:var1}{cmd:#i.}{it:var2}).{p_end}
{synopt:{it:var1}{cmd:#}{cmd:c.}{it:var2}}indicators for each level of {it:var1}, multiplied by {it:var2}{p_end}
{synopt:{it:var1}{cmd:##}{cmd:c.}{it:var2}}equivalent to "{cmd:i.}{it:var1} {cmd:i.}{it:var1}{cmd:#}{cmd:c.}{it:var2}", but {it:much} faster (the two sets of fixed effects are absorbed jointly at each iteration){p_end}
{synoptline}
{p2colreset}{...}

{pstd}{it: Notes:}{p_end}
{p 5 7 2}
- Each {it:absvar} in the {it:absvars} list represents a fixed effect that you wish to absorb (like {it:individual}, {it:firm} or {it:time}).
{p_end}
{p 5 7 2}
- It is good practice to put the absvars with more dimensions first.
{p_end}
{p 5 7 2}
- Interactions (e.g. {it:x{cmd:#}z}) are supported. Using categorical interactions is faster than running {it:egen group(...)} beforehand.
{p_end}
{p 5 7 2}
- To partial-out fixed {it:slopes} (and not just fixed intercepts), use continuous interactions (e.g. {it:x{cmd:#}{cmd:c.}z}).
{p_end}
{p 5 7 2}
- Each {it:absvar} can contain any number of categorical interactions 
(e.g. {cmd:i.}{it:var1}{cmd:#i.}{it:var2}{cmd:#i.}{it:var3})
but at most one continuous interaction
(thus, {cmd:i.}{it:var1}{cmd:#c.}{it:var2}{cmd:#c.}{it:var3} is not allowed).
{p_end}
{p 5 7 2}
- The first {it:absvar} cannot contain a continuous variable ({cmd:i.}{it:var1}{cmd:#c.}{it:var2} is not allowed, although {cmd:i.}{it:var1}{cmd:##c.}{it:var2} is ok).{p_end}
{p 5 7 2}
- When saving fixed effects and using {cmd:##} interactions, remember that {newvar}{cmd:=}{it:varname1}{cmd:##c.}{it:varname2} will be expanded to
"{newvar}{cmd:=}{it:varname1} {newvar:_slope}{cmd:=}{it:varname1}{cmd:#c.}{it:varname2}"
{p_end}

{it:Summary of Options:}

{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Model}
{p2coldent:* {opth a:bsorb(hdfe##absvar:absvars)}   }identifiers of the fixed effects that will be absorbed{p_end}
{synopt :{opt g:enerate(stubname)}}will not overwrite the variables; instead creating new demeaned variables with the {it:stubname} prefix{p_end}
{synopt :{opt clear:}}will overwrite the dataset; leaving the transformed variables, as well as some ancillary ones (such as the fixed effects, weights, cluster variables, etc.).
Use {cmd: char list} to see details of those ancillary variables.
{p_end}
{synopt: {opt cluster:vars(varlist)}}list of variables containing cluster categories. This is used to give more accurate number of degrees of freedom lost due to the fixed effects, as reported on r(df_a).{p_end}
{synopt: {opt dropsi:ngletons}}remove singleton groups from the sample; once per {it:absvar}.{p_end}
{synopt: {opt sample:(newvarname)}}will save the equivalent of e(sample) in this variable;
useful when dropping singletons.
Used with the {opt g:enerate} option.{p_end}
{synopt: {opth cores(#)}}will run the demeaning algorithm in # parallel instances.{p_end}
{synopt :{opt v:erbose(#)}}amount of debugging information to show (0=None, 1=Some, 2=More, 3=Parsing/convergence details, 4=Every iteration){p_end}
{synopt :{opth maxit:erations(#)}}specify maximum number of iterations; default is {cmd:maxiterations(1000)}; 0 means run forever until convergence{p_end}
{synopt :{it:maximize_options}}there are several advanced maximization options, useful for tweaking the iteration. See the {help reghdfe##maximize_options:help for reghdfe} for details.{p_end}

{marker recovering}{...}
{title:Recovering Fixed Effects}

{pstd}You can use {cmd:hdfe} again to recover the fixed effects. For instance, in the least-squares case:{p_end}

{phang2}{cmd:. sysuse auto, clear}{p_end}
{phang2}{cmd:. * Demean variables}{p_end}
{phang2}{cmd:. hdfe price weight length, a(turn trunk) gen(RESID_)}{p_end}
{phang2}{cmd:. * Run regression}{p_end}
{phang2}{cmd:. reg RESID_*}{p_end}
{phang2}{cmd:. * Predict using original variables}{p_end}
{phang2}{cmd:. drop RESID_*}{p_end}
{phang2}{cmd:. rename (price weight length) RESID_=}{p_end}
{phang2}{cmd:. predict double resid, resid}{p_end}
{phang2}{cmd:. rename RESID_* *}{p_end}
{phang2}{cmd:. * Obtain fixed effects}{p_end}
{phang2}{cmd:. hdfe resid, a(FE1=turn FE2=trunk) savefe gen(temp_)}{p_end}

{phang2}{cmd:. * Benchmark and verification}{p_end}
{phang2}{cmd:. reghdfe price weight length, a(BENCH1=turn BENCH2=trunk)}{p_end}
{phang2}{cmd:. gen double delta = abs(BENCH1-FE1) + abs(BENCH2-FE2)}{p_end}
{phang2}{cmd:. su delta}{p_end}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:hdfe} stores the following in {cmd:r()}:

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars}{p_end}
{synopt:{cmd:r(df_a)}}degrees of freedom lost due to the fixed effects (taking into account the cluster structure and whether the FEs are nested within the clusters){p_end}
{synopt:{cmd:r(N_hdfe)}}number of sets of fixed effects{p_end}
{synopt:{cmd:r(df_a#)}}degrees of freedom lost due to the #th fixed effect (excluding those collinear with the #th-1 first FEs){p_end}

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Macros}{p_end}
{synopt:{cmd:r(hdfe#)}}canonical expansion of the fixed effects (e.g. for#turn is expanded into i.foreign#i.turn){p_end}

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
{cmd:reghdfe} and {cmd:hdfe} are updated frequently, and upgrades or minor bug fixes may not be immediately available in SSC.
To check or contribute to the latest version of reghdfe, explore the
{browse "https://github.com/sergiocorreia/reghdfe":Github repository}.
{p_end}

