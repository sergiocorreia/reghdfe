## Author

Sergio Correia. _Fuqua School of Business, Duke University_

For any comments or suggestions, feel free to contact me at: [sergio.correia@duke.edu](mailto:sergio.correia@duke.edu)

## Acknowledgments

This package wouldn't have existed without the invaluable feedback and contributions of Paulo Guimarães
 and Amine Quazad. I am also indebted to the guidance of Kit Baum and Mark Schaffer, and to the great bug-spotting abilities of many users.

In addition, `reghdfe` uses several important contributions from the Stata community:

- [reg2hdfe](https://ideas.repec.org/c/boc/bocode/s457101.html), from Paulo Guimaraes, and [a2reg](https://ideas.repec.org/c/boc/bocode/s456942.html) from Amine Quazad,
 were the inspiration and building blocks on which reghdfe was built.
- [ivreg2](http://www.repec.org/bocode/i/ivreg2.html), by Christopher F Baum, Mark E Schaffer and Steven Stillman, is the package used by default for instrumental-variable regression.
- [avar](https://ideas.repec.org/c/boc/bocode/s457689.html) by Christopher F Baum and Mark E Schaffer, is the package used for estimating the HAC-robust standard errors of ols regressions.
- [tuples](http://econpapers.repec.org/software/bocbocode/s456797.htm) by Joseph Lunchman and Nicholas Cox, is used when computing standard errors with multi-way clustering (two or more clustering variables).
- [parallel](https://ideas.repec.org/c/boc/bocode/s457527.html) by George Vega Yon is used when running reghdfe with multiple processors.

## License

The MIT License (MIT)

Copyright (c) 2014 Sergio Correia

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
