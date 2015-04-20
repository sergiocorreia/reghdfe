
# Singletons and Cluster-Robust Standard Errors

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
Notes: Number of iterations = 100.

As we can see, with clustered standard errors and including singleton observations, `areg` will over-reject and `xtreg,fe` will under-reject. On the other hand, excluding singleton observations and/or using bootstrapped standard errors will fix/ameliorate the problem.

(TODO) Moderate Scenario


## Why Singletons Affect The Bias of Clustered SEs

(TODO)
1. Use formulas for robust standard errors
2. Explain how nested adjustment messes things up

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

