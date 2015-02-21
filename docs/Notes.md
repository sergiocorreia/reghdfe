

## Options

...

### Do not add back the constant with `nocons`

After partialling-out the fixed effects, we usually add back the constant term. This gives results consistent with [-areg-](http://stackoverflow.com/questions/14179197/how-to-interpret-the-constant-in-an-areg-output) and -xtreg, fe-.

To disable this, use the option <em><ul>nocon</ul>stant</em>. Note that by construction, all -reghdfe- regressions will include a constant through the absorb() terms. Also note that when using -ivreg2- the constant will always be omitted because otherwise -ivreg2- will complain about the "orthogonality conditions S is not of full rank" (when using clustered standard errors).
