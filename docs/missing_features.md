# Missing Features

(as of May 27 latest branch)

### Recovering Fixed Effects and additional functions for those

- Latest `reghdfe` has been optimized for recovering betas (e.g. using cg+sym.kac), but not for recovering the alphas i.e. the FEs (where it's using sd+sym which is even slower than Aitken's as used on v2). A Kaczmarz solver like R's `lfe` would significantly speed things up.
- As I'm not doing research related to the alphas, it has very low priority.
- There are several other useful statistics and results that can be reported/obtained about the alphas (e.g. bootstrapping their SEs, etc.) but again it's low priority for me.

### Optimization

- The latest optimization (GT) is not enabled/coded. So far the four implemented theoretical tricks (FWL+MAP+ACC+SYM) are enough for me, but GT can potentially bring linear time even in worst-case scenarios.
- However at this point other slowdowns may occur, so we may end up in a case where sorting the data is the bottleneck, and there is not much to do about it.
- Also, GT may not be possible with `c.` interactions, so the scope may be limited.

### Benchmarking

- Existing tests (e.g. M.Gomez) are "best case scenarios" as the subspaces spanned by each set of FEs are orthogonal by construction. However, most of the practical difficulties lie in solving the model when the FEs have a high angle (this can be thought in many other ways, such as the ratio of the eigenvalues of the laplacian being high, the underlying graph being poorly connected, etc.)
- Thus, writing a good benchmark that is closer to "real life" would be useful for optimizing those scenarios

### HDFE

- Need to update hdfe.ado to use latest improvements (reghdfe.ado is in v3 but hdfe is still in v2)

### Documentation

Need to update/extend the following:

- Help file
- Write tutorial
- Update website with help files and key explanations (why drop singletons, why no constant, what not to do regarding SEs, why alphas are often not identified, etc.)
- Document how the algorithm works, as it's now quite different from `reg2hdfe` and even `lfe`.

### Maybe

- Explain how to produce simple tables from `reghdfe` (e.g. how to annotate the list of included FEs)
- Links with quipu.ado and markdown/pandoc

