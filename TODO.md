# Pending tasks

## High priority:

- Fix open issues
- ivreg2 example
	- basic demo
	- allow passthrough of optimization options through `absorb(, ...)` (just set the absorb parser to return a `r(rest)`)
	- stages
- Optimization
	- enable prune code for G=2 *and* G>2 (carefully)
	- (?) allow weights with cg+sd, in their scalar adjustements (does hestenes has weights now?)
	- Add LSMR as with Mathieu's implem.
	- Improve the case with id#c.var  and id##c.var . Currently it relies on a lot of calls to invsym. This is unnecessary, as we can i) do it once and save it, or ii) use something better than actually inverting the matrix, and saving the steps (see qrsolve or cholsolve)
- Cleanup
	- Remove code vestiges (ivsuite, vcesuite, method, etc.)
- Docs
	- Update help file
	- Fix html help builder
	- update website
- Compute alphas

## Lower priority:

- Memory usage:
	- run `st_data` and `partial_out` in chunks, with the `pool(#)` option
	- `compact` option: preserve, clear, and restore the dataset to save memory (slower but might be necessary for large datasets)
	- allow `nosample` option
- Fix inference
	- wild bootstrap (boottest) helps if firm and CEO numbers are growing at sqrt(N) rate
	- Verdier's approach (or a feasible implementation of Cattaneo et al) help if firm numbers or CEO numbers grow proportional to N
	- The first case can be seen as an expanding interconnected economy
	- The second, as adding more independent or poorly connected countries to a dataset
- hdfe.ado
- Improve degrees-of-ofreedom with group3hdfe (or Mata implem of it)
- demo of how poi2hdfe could work
- `savecache` and `usecache` (or `hold`) options?




# LSMR References

http://web.stanford.edu/group/SOL/software/lsmr/

http://www.mathworks.com/matlabcentral/fileexchange/27183-lsmr--an-iterative-algorithm-for-least-squares-problems

https://github.com/timtylin/lsmr-SLIM
https://github.com/timtylin/lsmr-SLIM/blob/master/lsmr.m
https://github.com/timtylin/lsmr-SLIM/blob/master/lsmrtest.m

https://github.com/PythonOptimizers/pykrylov/tree/master/pykrylov/lls
https://github.com/PythonOptimizers/pykrylov/blob/master/pykrylov/lls/lsmr.py

http://web.stanford.edu/group/SOL/software/craig/ (we could use CRAIG for the alphas?)
