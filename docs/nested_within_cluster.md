
# Singletons and Cluster-Robust Standard Errors

## Summary

There is nothing to gain from including singleton groups in a regression. On the other hand, there are three reasons why excluding those observations is useful:

1. When using clustered standard errors and fixed effects, including singletons with commands such as `xtreg,fe` may underestimate the standard errors and an overestimate the P-Values associated with the parameters of interest.
This is problem is particularly problematic when running regressions with high-dimensional fixed effects (HDFE). In those cases, a large number of observations may be perfectly predicted (in-sample) just by the fixed effects. For instance, in matched CEO-Firm regressions many individuals and firms may be short-lived enough in the sample so that singletons abound.
(TODO) If those fixed effects are nested within clusters, the degrees of freedom will not be adjusted ...
2. Another benefit of removing singletons is that it shows the *effective* number of clusters, which helps in knowing if the cluster size is large enough for the asymptotics to kick in.
3. Finally, excluding singletons speeds up the computation of HDFE regressions, as it reduces the number of ancilliary parameters that are estimated.

### Example

As an extreme but illustrative example of the first problem, consider the following regressions using the sample Stata dataset:

```stata
sysuse auto, clear
gen id = _n
replace id = id-1 if _n<8 & mod(id,2)==0
bys id: gen t = _n
xtset id t
bys id: gen is_singleton = (_N==1)
tab is_singleton

xtreg price weight length, fe vce(cluster id)
xtreg price weight length, fe vce(cluster id) dfadj
drop if is_singleton
xtreg price weight length, fe vce(cluster id)
xtreg price weight length, fe vce(cluster id) dfadj
```

The first regression reports a P-Value of 0.007 for the *weight* regressor, while the subsequent regressions (those dropping singletons and/or subtracting the fixed effects from the degrees-of-freedom) report much higher P-Values ranging from 0.212 to 0.796.

## Singletons and Fixed-Effects

Singletons are individuals that appear in only one observation. If we are using individual fixed-effects, then the fixed effect will perfectly explain the value of the dependent variable for the singleton individuals. Therefore, *singletons do not contribute any useful information when estimating the regressors of interest*.

## Computing Clustered Standard Errors

(TODO)
1. Show formula
2. Mention that cluster is the same as robust when clustervar == _n
3. Show the small sample adjustments and how they differ in xtreg and areg

## How Singletons Affect the Bias of Clustered SEs

Singleton individuals have no effect on the estimates of the coefficients, but they *do* affect the estimates of the standard errors (and thus those of the pvalues).

A simulation on an [extreme scenario](https://github.com/sergiocorreia/reghdfe/blob/master/misc/example_nested_bug.do), with 100 observations distributed between i) 90 singletons individuals, and i) 5 singleton individuals with two obs. each, showed that robust standard errors do differ slightly, but clustered standard errors diverge *massively* from their correct values, accepting the null in many cases.

Extreme Scenario (v1)

|   Estimator  |          Sample         |      S.E.     | % with Pvalue < 5%    | % with Pvalue < 10%   |
|:------------:|:-----------------------:|:-------------:|:---------------------:|:---------------------:|
| areg & xtreg | full & w/out singletons |   unadjusted  |                     % |                  4.0% |
|     areg     |           full          |     robust    |                     % |                  0.0% |
|     xtreg    |           full          |     robust    |                     A |                   N/A |
|     areg     |     w/out singletons    |     robust    |                     % |                 15.0% |
|     xtreg    |     w/out singletons    |     robust    |                     A |                   N/A |
|     areg     |           full          |   clustered   |                     % |                    0% |
|     xtreg    |           full          |   clustered   |                     % |                 33.0% |
|     areg     |     w/out singletons    |   clustered   |                     % |                    5% |
|     xtreg    |     w/out singletons    |   clustered   |                     % |                 14.0% |
|     xtreg    |           full          | cl. bootstrap |                     % |                  5.0% |
|     xtreg    |     w/out singletons    | cl. bootstrap |                     % |                 12.0% |
Notes: Number of iterations = 100. 100 observations. 10 non-singleton observations (5 groups of 2 obs. each).


As we can see, with clustered standard errors and including singleton observations, `areg` will over-reject and `xtreg,fe` will under-reject. On the other hand, excluding singleton observations and/or using bootstrapped standard errors will fix/ameliorate the problem.

(TODO) Moderate Scenario (v2)

|   Estimator  |          Sample         |      S.E.     | % with Pvalue < 5%    | % with Pvalue < 10%   |
|:------------:|:-----------------------:|:-------------:|:---------------------:|:---------------------:|
| areg & xtreg | full & w/out singletons |   unadjusted  |                  3.5% |                 11.5% |
|     areg     |           full          |     robust    |                  0.0% |                  2.0% |
|     xtreg    |           full          |     robust    |                   N/A |                   N/A |
|     areg     |     w/out singletons    |     robust    |                  4.5% |                 12.5% |
|     xtreg    |     w/out singletons    |     robust    |                   N/A |                   N/A |
|     areg     |           full          |   clustered   |                  0.0% |                  0.0% |
|     xtreg    |           full          |   clustered   |                  4.5% |                 12.5% |
|     areg     |     w/out singletons    |   clustered   |                  0.0% |                  2.0% |
|     xtreg    |     w/out singletons    |   clustered   |                  4.5% |                 12.5% |
|     xtreg    |           full          | cl. bootstrap |                  5.0% |                 13.0% |
|     xtreg    |     w/out singletons    | cl. bootstrap |                  4.0% |                 13.0% |
Notes: Number of iterations = 500. 200 observations. 100 non-singleton observations (50 groups of 2 obs. each).


## Why Singletons Affect The Bias of Clustered SEs

(TODO)
1. Use formulas for robust standard errors
2. Explain how nested adjustment messes things up
3. Mention the need for asymptotics in the number of clusters.

## References

[Matsa's post about dealing with fixed effects](http://www.kellogg.northwestern.edu/faculty/matsa/htm/fe.htm)
	
	"XTREG’s approach of not adjusting the degrees of freedom is appropriate
	when the fixed effects swept away by the within-group transformation are 
	nested within clusters (meaning all the observations for any given group 
	are in the same cluster), as is commonly the case (e.g., firm fixed 
	effects are nested within firm, industry, or state clusters). 
	See Wooldridge (2010, Chapter 20)."

[A. Colin Cameron and Douglas L. Miller, "A Practitioner's Guide to Cluster-Robust Inference", Journal of Human Resources, forthcoming, Spring 2015](http://cameron.econ.ucdavis.edu/research/Cameron_Miller_Cluster_Robust_October152013.pdf)

[Mark Schaffer's Statalist post](http://www.stata.com/statalist/archive/2006-07/msg00535.html)

	"In this panel data context, a singleton is a group in which there is
	only one observation.  Since singletons have zero within-group
	information, the within (demeaning) transformation will zap them.

	Stata's official commands that do linear fixed effects estimation
	(xtreg, xtivreg, areg) do not adjust the number of observations for the
	singletons.  Explicitly excluding singletons can therefore affect the
	SEs but will leave the coefficients unchanged

	... it is correct to treat singletons as non-observations, no different 
	from observations that are lost because of missing values ..."

