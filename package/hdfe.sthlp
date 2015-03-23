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
[{opth cluster:vars(hdfe##vcetype:vcetype)}
{opt v:erbose(#)}
{opt tol:erance(#)}
{opt maxi:terations(#)}
{opt g:enerate(stubname)}]
{p_end}

{pstd}{it: Notes:}{p_end}
{p 5 7 2}
- does not accept time series or factor variables
{p_end}
{p 5 7 2}
- variables MUST BE FULLY spelled out (i.e. you need to use unab beforehand!)
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
{synopt: {opt cluster:vars(varlist)}}list of variables containing cluster categories{p_end}
{synopt :{opt v:erbose(#)}}amount of debugging information to show (0=None, 1=Some, 2=More, 3=Parsing/convergence details, 4=Every iteration){p_end}
{synopt :{opth maxit:erations(#)}}specify maximum number of iterations; default is {cmd:maxiterations(1000)}; 0 means run forever until convergence{p_end}
{synopt :{opt g:enerate(stubname)}}will not overwrite the variables; instead creating new demeaned variables with the {it:stubname} prefix{p_end}

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

