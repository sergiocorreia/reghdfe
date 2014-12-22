# REGHDFE: Linear and IV Regressions With Many Fixed Effects

This is the *readme* file for developing the reghdfe project. What you probably want is the tutorial/documentation, available [here](TODO) (soon).

* Version 1.2.1
* Date: December 9, 2014

Reghdfe is a [Stata](http://stata.com) package that estimates linear and instrumental variable regressions while controlling for any number of fixed effects. It is most useful when dealing with two or more highly dimensional fixed effects (HDFE).

It extends the algorithm of *reg2hdfe* (by Paulo Guimaraes and Pedro Portugal) in several aspects, **drastically speeding it up**, as well as **generalizing** it to any number of sets of fixed effects, slope interactions, multicore, etc.

It's objectives are similar to the *lfe* package in R by Simen Gaure, and to several precursor Stata packages, such as *a2reg* by Amine Quazad.

## Features

* Can estimate any number of highly-dimensional fixed effects.
* Extended to allow for intercept and/or slope effects by group (i.e. by absorbing variable).
* Allows instrumental variable regressions by using the *ivregress* and *ivreg2* packages.
* Can easily save the fixed effect estimates (whether or not they are identified is another matter..)
* Very fast by default. Controlling for many fixed effects is usually *extremely* slow and memory consuming. Therefore, the code has been carefully profiled to avoid any bottlenecks within the Stata and Mata code.
* Even faster experimental options: Although Mata does not provide custom multicore support, it can be achieved by opening multiple Stata instances and parallelizing by variable. This has been tested only in Windows, but can be extended to OSX/Linux easily. There are also other optimization options that can further speed it up.
* Attention to the detail. For instance, you can just include firm#year instead of writing *egen firm_year = group(firm year)*. This works in the varlist, in the absvars, and even with the vce(cluster ..) option.
* Many other options such as analytic and frequency weights, detailed debugging options, etc.

## Installing the development version (i.e. Github instead of SSC)

```stata
cap ado uninstall reghdfe
net from https://raw.githubusercontent.com/sergiocorreia/reghdfe/master/package/
net install reghdfe
```
**Note: this requires Stata 13, as older versions of stata cannot handle *https*. If you have an older version, see below.**

If you have an older version of Stata, first copy the contents of the [package folder](/package) to a folder on your computer, and then replace the second line above with `net from My/Path/` where My/Path is the path where you saved the files. Note that with this workaround you will not be able to use `ado update` to automatically update from github. If you don't want to download each file by hand, [download the zip](https://github.com/sergiocorreia/reghdfe/archive/master.zip).

## Contributing

Contributors and pull requests are more than welcome. There are a number of extension possibilities, such as adding two-way clustering, multicore support for OSX and Linux, estimating standard errors for the fixed effects using bootstrapping, and exact computation of degrees-of-freedom for more than two HDFEs.

## Building the package

For sanity reasons, the source code is spread through several files and folders (in the *source* folder). To update the package, do the following:

* Download the entire project to your computer (through the "Clone Desktop" or "Download ZIP" buttons on the right).
* Uninstall any existing versions of *reghdfe* (`ado uninstall reghdfe` in Stata).
* In Stata, change the current working directory to the *source* folder. Do any changes that you want on the files in that folder. You can run *reghdfe* without problems as long as the working directory is in that folder.
* To build the package, run the *build.py* file (in the *build* folder), using either Python 2.7 or Python 3.x (3.x is not tested but should work). This python script will carefully combine all the files and update the version/date.
* Finally, you can upload it back to github and submit a pull request.

## Author

[Sergio Correia](sergio.correia@gmail.com), Duke University

## Acknowledgments

This package wouldn't have existed without the invaluable feedback and contributions of Paulo Guimaraes and Amine Quazad, as well as the great bug-spotting abilities of many users.
