# REGHDFE: Linear and IV/GMM Regressions With Many Fixed Effects

`reghdfe` is a [Stata](http://www.stata.com/) package that estimates linear regressions with multiple levels of fixed effects. It works as a generalization of the built-in `areg`, `xtreg,fe` and `xtivreg,fe` regression commands. It's objectives are similar to the R package [lfe](http://cran.r-project.org/web/packages/lfe/index.html) by Simen Gaure. It's features include:

- A novel and robust algorithm that efficiently absorbs multiple fixed effects. It improves on the work by [Abowd *et al*, 2002](https://ideas.repec.org/p/cen/tpaper/2002-06.html), [Guimaraes and Portugal, 2010](https://ideas.repec.org/a/tsj/stataj/v10y2010i4p628-649.html) and [Simen Gaure, 2013](http://www.sciencedirect.com/science/article/pii/S0167947313001266). This algorithm works particularly well on "hard cases" that converge very slowly (or fail to converge) with the existing algorithms.
- Extremely fast compared to similar Stata programs. 
  - With one fixed effect and clustered-standard errors, it is 3-4 times faster than `areg` and `xtreg,fe` (see [benchmarks](./misc/Benchmarks/areg_xtreg.log.txt)). Note: speed improvements in Stata 14 have reduced this gap. 
  - With multiple fixed effects, it is at least an order of magnitude faster that the alternatives (`reg2hdfe`, `a2reg`, `felsdvreg`, `res2fe`, etc.). Note: a recent paper by [Somaini and Wolak, 2015](http://web.stanford.edu/group/fwolak/cgi-bin/sites/default/files/jem-2014-0008.pdf) reported that `res2fe` was faster than `reghdfe` on some scenarios (namely, with only two fixed effects, where the second fixed effect was low-dimensional). This is no longer correct for the current version of `reghdfe`, which outperforms `res2fe` even on the authors' benchmark (with a low-dimensional second fixed effect; see the [benchmark results](./misc/Benchmarks/res2fe.log.txt) and the Stata [code](./misc/Benchmarks/res2fe.do)).
- Allows two- and multi-way clustering of standard errors, as described in [Cameron *et al* (2011)](http://amstat.tandfonline.com/doi/abs/10.1198/jbes.2010.07136)
- Allows an extensive list of robust variance estimators (thanks to the [avar](https://ideas.repec.org/c/boc/bocode/s457689.html) package by Kit Baum and Mark Schaffer).
- Works with instrumental-variable and GMM estimators (such as two-step-GMM, LIML, etc.) thanks to the [ivreg2](https://ideas.repec.org/c/boc/bocode/s425401.html) routine by Baum, Schaffer and Stillman.
- Allows multiple heterogeneous slopes (e.g. a separate slope coefficients for each individual).
- Supports all standard Stata features:
  - Frequency, probability, and analytic weights.
  - Time-series and factor variables.
  - Fixed effects and cluster variables can be expressed as factor interactions, for both convenience and speed (e.g. directly using `state#year` instead of previously using `egen group` to generate the state-year combination).
  - Postestimation commands such as `predict` and `test`.
- Allows precomputing results with the `cache()` option, so subsequent regressions are faster.
- If requested, saves the point estimates of the fixed effects (*caveat emptor*: these fixed effects may not be consistent nor identifiable; see the Abowd paper for an introduction to the topic).
- Calculates the degrees-of-freedom lost due to the fixed effects (beyond two levels of fixed effects this is still an open problem, but we provide a conservative upper bound).
- Avoids common pitfalls, by excluding singleton groups (see [notes](scorreia.com/reghdfe/nested_within_cluster.pdf)), computing correct within- adjusted-R-squares ([see initial discussion](http://www.statalist.org/forums/forum/general-stata-discussion/general/1290416-anyone-knows-what-is-an-adjusted-within-r2)), etc.

## Author

[Sergio Correia](http://scorreia.com)
<br>Fuqua School of Business, Duke University
<br>Email: sergio.correia@duke.edu

## Acknowledgments

This package wouldn't have existed without the invaluable feedback and contributions of Paulo Guimaraes, Amine Ouazad, Mark Schaffer and Kit Baum. Also invaluable are the great bug-spotting abilities of many users.

## Citation

`reghdfe` is a free contribution to the research community, like a paper. Please cite it as such:

> Sergio Correia, 2015. *reghdfe: Stata module for linear and instrumental-variable/GMM regression absorbing multiple levels of fixed effects.*
>
> https://ideas.repec.org/c/boc/bocode/s457874.html

## Description

This is the *readme* file for developing the reghdfe project, which is comprised of the reghdfe package and the underlying hdfe package. The help files and tutorials [are available here](http://scorreia.com/reghdfe) (work in progress).

Latest version
* Version 3.1.9
* Date: June 15, 2015

## Installing REGHDFE

The latest stable release (2.1.x) can be downloaded from SSC with

```stata
cap ado uninstall reghdfe
ssc install reghdfe
```

The installation of the latest dev. release (3.1.x) depends on the Stata version:

With Stata 13:

```stata
cap ado uninstall reghdfe
net install reghdfe, from(https://raw.githubusercontent.com/sergiocorreia/reghdfe/master/package/)
```
With Stata 12 or older:

1. Download the [zipfile](/misc/reghdfe.zip?raw=true)
2. Extract it into a folder (e.g. C:\SOMEFOLDER)
3. Run: (changing *SOMEFOLDER* with whatever you picked)

```stata
cap ado uninstall reghdfe
net install reghdfe, from("C:\SOMEFOLDER")
```

To find out which version you have installed, type `reghdfe, version`.

## Installing HDFE

`hdfe` is a routine that facilitates absorbing multiple fixed effects in other Stata packages. It is similar to `avar` in that it is a building-block routine that other packages may call (for instance, see [regife](https://github.com/matthieugomez/stata-regife) and [poi2hdfe](https://ideas.repec.org/c/boc/bocode/s457777.html))

The latest stable release (2.1.x) can be downloaded from SSC with

```stata
cap ado uninstall hdfe
ssc install hdfe
```

The installation of the latest developer release (3.1.x) depends on the Stata version:

With Stata 13:

```stata
cap ado uninstall hdfe
net install hdfe, from(https://raw.githubusercontent.com/sergiocorreia/reghdfe/master/package/)
```
With Stata 12 or older:

1. Download the [zipfile](/misc/reghdfe.zip?raw=true)
2. Extract it into a folder (e.g. C:\SOMEFOLDER)
3. Run: (changing *SOMEFOLDER* with whatever you picked)

```stata
cap ado uninstall hdfe
net install hdfe, from("C:\SOMEFOLDER")
```

## Recent Updates

* 3.1 Improved syntax for the `cache()` and `stage()` options
* 3.0 Three key changes: i) faster underlying algorithm (symmetric transforms and cg acceleration perform much better on "hard" cases); ii) slow parts rewritten in mata, iii) simpler syntax
* 2.2 [internal] murphy-topel (unadjusted, robust, cluster), double-or-nothing IV/control function
* 2.1 removed `_cons`. If you really want to see the constant, run *summarize* on the first fixed effect. The last version that supported constants is available with `net from https://raw.githubusercontent.com/sergiocorreia/reghdfe/866f85551b77fe7fda2af0aafccbbf87f8a01987/package/`

## Future/possible updates

* 4.0 Improve underlying algorithm with GT preconditioning
* 5.0 Increase features for recovering the fixed effects. For instance, bootstrapping the standard errors, a better algorithm (Kaczmarz) for recovering the point estimates, and a wider set of statistics for the standard errors. If you currently require any of those, I recomment the [lfe suite](cran.r-project.org/web/packages/lfe/index.html) by Simen Gaure (for the R language).

## Contributing

Contributors and pull requests are more than welcome. There are a number of extension possibilities, such as estimating standard errors for the fixed effects using bootstrapping, exact computation of degrees-of-freedom for more than two HDFEs, and further improvements in the underlying algorithm.

## Building the package

For clarity reasons, the source code is spread through several files and folders (in the [source](./source) folder). To modify and rebuild the package, do the following:

* Download the entire project to your computer (through the "Clone Desktop" or "Download ZIP" buttons on the right).
* Uninstall any existing versions of *reghdfe* (`ado uninstall reghdfe` in Stata).
* Do any changes that you want on the files in that folder. You can run *reghdfe* without problems as long as the working directory is in that folder.
* To build the package, run the *build.py* file (in the *build* folder), using Python 3.x. This python script will carefully combine all the files and update the version/date.
* Install it using `net install reghdfe, from(PATH_OF_THE_PACKAGE_FOLDER)`
* Finally, you can upload it back to github and submit a pull request.
