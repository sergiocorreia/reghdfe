

The finite-sample estimates of the asymptotic variance-covariance matrix of estimates $V$ is equal to $V* = q \times V$, where



There is nothing to gain from including singleton groups in a regression. On the other hand, there are three reasons why excluding those observations is useful:

1. When using clustered standard errors and fixed effects, including singletons with commands such as `xtreg,fe` may underestimate the standard errors and an overestimate the P-Values associated with the parameters of interest.
This is problem is particularly problematic when running regressions with high-dimensional fixed effects (HDFE). In those cases, a large number of observations may be perfectly predicted (in-sample) just by the fixed effects. For instance, in matched CEO-Firm regressions many individuals and firms may be short-lived enough in the sample so that singletons abound.
(TODO) If those fixed effects are nested within clusters, the degrees of freedom will not be adjusted ...
2. Another benefit of removing singletons is that it shows the *effective* number of clusters, which helps in knowing if the cluster size is large enough for the asymptotics to kick in.
3. Finally, excluding singletons speeds up the computation of HDFE regressions, as it reduces the number of ancilliary parameters that are estimated.




## Singletons and Fixed-Effects

Singletons are individuals that appear in only one observation. If we are using individual fixed-effects, then the fixed effect will perfectly explain the value of the dependent variable for the singleton individuals. Therefore, *singletons do not contribute any useful information when estimating the regressors of interest*.

## Computing Clustered Standard Errors

In general,
$$
V = q * {X'X}^{-1} \times (S'S) \times {X'X}^{-1}
q = M / (M-1) \times (N-1) / (N-K-L)
$$

Where $q$ is a finite-sample adjustment, $M$ is the number of clusters, $L$ is the number of fixed effects, and $K$ is the number of other regressors excluding the constant (which is grouped with the fixed effects). Note that if the fixed effects are nested within the clusters, $L$ is excluded from $q$.

Also, $X*$ is the matrix of residuals of X over the fixed effects, $S$ is the score matrix, i.e. the sum of $X* . e$ (dot product between regressors and residuals), grouped over the clustered individuals.

Whenever there is a singleton group and the fixed effect is nested within a cluster, the following happens:
1. The corresponding row of $X*$ becomes zero as the residuals of $X$ wrt. the fixed effects are zero. 
2. The corresponding row of S becomes a zero row-vector due to both $X*$ and the residuals being zero for that obs.
3. Since $L$ is excluded from the denominator of $q$ (because it's nested), then that denominator does not change.
4. Thus, the only thing that changes is $M/(M-1)$ and $N$

In summary, results 1 and 2 mean that both the *bread* and the *meat* of the variance remain unchanged when adding singletons. Thus, only the $q$ (the *wooden stick* of the burger), changes.

A natural solution to the formula would be to remove the singleton clusters from $M$. Calling those $M_S$, then:
$$
q* = (M-M_S) / (M-M_S-1) \times (N-M_S-1) / (N-M_S-K-L*)
$$

(Where $L*$ represents the degrees of freedom not absorbed by the clusters, if any).

We also see that S gives us bounds to the size of the bias in the S.E.s. For instance, if (excluding singletons) M=10, N=20, K=1, M_S=10, then
q=1.09, while q*=1.1728, so standard errors should be 3.7% higher than reported.

TODO: Recheck what assumptions are needed for always excluding L from q!

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

Moderate Scenario (v2)

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


## Finite-Sample Adjustments

(This follows McKinnon & White)

Hinkley (1977) proposed HC1: $q=N/(N-1)$
HHD (1975) proposed HC2, that replaces the meat instead: `cross(X,e:^2,X)` becomes `cross(X,e:^2 :/ (1-diag(P),X)` where $P$ is the projection matrix. This is unbiased with homoskedasticity, which is a very nice property. Note that HC1 is usually biased.


## Thought Experiments

1. Expand the dataset, cloning each obs. by e.g. 10. In that case, both meat and bread stay unchanged

```
sysuse auto, clear
bys turn: gen t = _n
gen index = _n
expand 10
bys turn t: gen clone = _n
egen id = group(clone turn)
xtset id t
xtreg price weight length if clone==1, fe vce(cluster turn)
xtreg price weight length, fe vce(cluster turn)
di 2.083286 * sqrt(73/739 / 71 * 737)
```
As we can see, the finite-sample adjustment used is imperfect and could be improved by exploiting the fact that there is no within-cluster variability in the variables. In an extreme case, you could "cheat" in any regression by changing the unit of analysis from e.g. pair of shoes to single shoes, which wouldn't change anything except the finite-sample adjustment.

Can this and the singleton problems be fixed with an unbiased HC# estimator? (See McKinnon's survey)