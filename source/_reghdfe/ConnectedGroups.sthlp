{smcl}
{title:Title}

{p 4 15}ConnectedGroups - Computes the number of connected groups (i.e. redundant fixed effects) in a two-fixed effect model.


{title:Syntax}

{p 8 15}
{cmd:ConnectedGroups} {it:{help varname1}} {it:{help varname2}} [, {cmdab:gen:erate}{cmd:(}{it:newvarname}{cmd:)}} | {cmd:clear}]

{p}

{title:Description}

{p 4 4}This command calculates the number of redundant restrictions that can be substracted from the number of estimated parameters in a two fixed effects model. The result is saved in {it:r(groups)}.

{p 4 4}For performance and simplicity reasons, there are two key limitations / constrainst:

{p 4 4}1) -if- and -in- arguments cannot be specified. Delete unneeded obs. beforehand.

{p 4 4}2) It is preferable (but not required) for the FE variables to be in consecutive order (from 1 until the number of FEs defined by each variable). This can be achieved with the GenerateID command or with -egen group-. This helps in calculating the total number of fixed effects, which will then be substracted from the redundant fixed effects to obtain the number of estimated parameters and the degrees of freedom.

{p 4 4}This program is heavily based on -group3hdfe- by Paulo Guimaraes and -a2reg- by Amine Quazad. It runs in around a third of the time than the previous programs (mainly because it is less general).


{title:Options}

{p 4 4}{cmdab:gen:erate}{cmd:(}{it:newvarname}{cmd:)} Creates a variable with the identifiers for the connected groups. Not compatible with the {cmd:clear} option.

{p 4 4}{cmd:clear} Do not restore the data. Useful to avoid saving it, and if it will be of no further use afterwards. Not compatible with the {cmd:generate} option.

{title:Examples}

Obtain r(groups) but do not destroy the dataset
{p 8 16}{inp:. ConnectedGroups i j }{p_end}

Destroy the dataset, when we only want r(groups)
{p 8 16}{inp:. ConnectedGroups i j, clear}{p_end}

Merge the new group ID back into the dataset
{p 8 16}{inp:. ConnectedGroups i j, gen(group_id)}{p_end}

Keep the messed up dataset with the group ID. Not very useful:
{p 8 16}{inp:. ConnectedGroups i j, gen(group_id) clear}{p_end}

{title:Stored results}

{pstd}
{cmd:ConnectedGroups} stores the following in {cmd:r()}:

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars}{p_end}
{synopt:{cmd:r(groups)}}number of connected groups{p_end}


{title:Author}

{p}
Sergio Correia, Duke University

{p}
Email: {browse "mailto:sergio.correia@duke.edu":sergio.correia@duke.edu}

{title:Reference}

If you use this program in your research please cite Paulo's work:

Paulo Guimaraes and Pedro Portugal. "A Simple Feasible Alternative Procedure to Estimate Models with 
High-Dimensional Fixed Effects", Stata Journal, 10(4), 628-649, 2010.

{title:Also see}

{p 0 21}
{help hdfe} (if installed).
{p_end}
