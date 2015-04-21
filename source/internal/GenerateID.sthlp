{smcl}
{title:Title}

{p 4 15}GenerateID - Create unique identifiers based on one or more variables provided. This is just a faster alternative to {it:egen group}, but less general: i) it's does not accept variables with missing values, ii) it does destructive sorting


{title:Syntax}

{p 8 15}
{cmd:GenerateID} {it:{help varlist}} , [ {cmd:replace} | {cmdab:gen:erate}{cmd:(}{it:newvarname}{cmd:)}}] [{cmd:nocompress}]

{p}

{title:Description}

{p 4 4}This command creates a variable with values from 1 to N where N is the total number of unique values of the variables in {it:varlist}.

{p 4 4}If {it:varlist} only contains one variable, it can be used to make a more compact identifier (e.g. as an -int- or -long- instead of -double-), and the {cmd:replace} option can be used.

{p 4 4}Otherwise, a new variable name must be provided in the {cmd:generate()} option.

{title:Options}

{p 4 4}{cmdab:gen:erate}{cmd:(}{it:newvarname}{cmd:)} creates a variable containing the identifiers spanned by the {it:varlist}. Not compatible with the {cmd:replace} option{p_end}
{p 4 4}{cmd:replace} replace the provided variable, while trying to compress it. Not compatible with the {cmd:generate} option{p_end}
{p 4 4}{cmd:nocompress} prevents trying to compress the involved variables{p_end}

{title:Examples}

{p 8 16}{inp:. GenerateID id1, replace }{p_end}

Destroy the dataset, when we only want r(groups)
{p 8 16}{inp:. GenerateID id1 id2 id3, gen(new_id)}{p_end}
