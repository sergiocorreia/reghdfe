# Reghdfe Guide

`reghdfe` is a Stata package that runs linear and instrumental-variable regressions with many sets of fixed effects.

It solves regressions with any number of indicator variables by extending the fixed point iteration of Guimarães & Portugal (2010).

Within Stata, it can be viewed as a generalization of `areg`/`xtreg`, with several additional features:

* Any number of absorbed variables / fixed effects.
* Support for *fixed slopes* besides *fixed intercepts*.
* Two-stage least squares / instrumental-variable regression thanks to the `ivreg2` and `ivregress` commands.
* Multi-way clustering for 2+ cluster variables.
* Advanced options for computing standard errors, thanks to the `avar` command.
* Careful estimation of degrees of freedom, taking into account nesting of fixed effects within clusters, as well as many possible sources of collinearity within the fixed effects. 

In addition, it is easy to use and supports most Stata conventions:

* Time series and factor variable notation, even within the absorbing variables and cluster variables.
* Multicore support thanks to the -parallel- suite.
* Frequency weights, analytic weights, and probability weights are  allowed.
* It can cache results in order to run many regressions with the same data, as well as run regressions over several categories.
* It supports most post-estimation commands, such as `test`, `estat summarize`, and `predict`.
