---
name: Feature request
about: Suggest an idea for reghdfe
title: ''
labels: ''
assignees: ''

---

# Disclaimer

There are many features that I would definitely love to add to reghdfe. For instance:

Additional standard errors:

- [Conley](http://www.trfetzer.com/conley-spatial-hac-errors-with-fixed-effects/) spatial SEs
- SEs on the fixed effects
- Unbiased SEs under heteroskedasticity and high-dimensional fixed effects (Kline et al's leave-one-out estimation, Verdier's sparsely based data method, etc) that can be computed on a finite amount of time on a computer with non-NSA amounts of memory (see [this](https://scholar.google.com/scholar?rlz=1C1GGRV_enUS751US751&sxsrf=ACYBGNQZSBYZzQtI4PveBvnSTlxYGyLbfg:1575088778893&biw=1707&bih=886&um=1&ie=UTF-8&lr&cites=3283837062430382503) rapidly evolving list).

Econometric models:

- More advanced methods: logit fixed effects, GMM fixed effects
- LASSO and other penalized regression methods

Convenience tools:

- `cache()` options
- More examples in the help file

However, due to lack of "copious free time" I will generally will not be able add a new method/tool unless I can use it on my own research. That said, you are more than welcome to submit a pull request and I would be more than happy to help you on the steps to e.g. add different standard errors.

# Feature request

**Is your feature request related to a problem? Please describe.**

A clear and concise description of what the problem is. For instance, "bootstrap takes too long to run with reghdfe"

**Describe the solution you'd like**

*It would be easier if reghdfe saves temporary results and starts from there instead of from scratch for every bootstrap iteration*

**Additional context**

Add any other context or screenshots about the feature request here.
