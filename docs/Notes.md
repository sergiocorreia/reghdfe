Notar diferencia con xtivreg2:
Radica en el q que multiplica a la VCV
Ellos dicen "no es necesario ajustar por N/N-k-kk" porque N va al infinito wrt k,kk
Que es como decir que tienes suficientes unidades de tiempo

Subtle difference between xtivreg2 and reghdfe:
Both do the N/(N-k-kk) adjustment
(Note than in a normal ols, kk=1 i.e. the constant)
Also both recognize that since clustervar==absvar, kk=0
However, xtivreg2 replaces kk=1 instead of kk=0


When I explain DROPSINGletons , mention that it changes the estimated value of the constant

All in all, constants just seem to be a HUGE PITA
- conflict with over()
- mess up VCV sp. with ivreg2 and avar with clustering
- Also it's interpretation is not particularly useful (it's not the avg. of depvar, just the expected value when all indepvars are zero)
- Constant also changes when dropping singleton groups!




In the helpfile, carefully explain what goes in vce() and absorb()
and *strongly* recommend people to use the interactions # as it's faster in several parts

## Options

...

### Do not add back the constant with `nocons`

After partialling-out the fixed effects, we usually add back the constant term. This gives results consistent with [-areg-](http://stackoverflow.com/questions/14179197/how-to-interpret-the-constant-in-an-areg-output) and -xtreg, fe-.

To disable this, use the option <em><ul>nocon</ul>stant</em>. Note that by construction, all -reghdfe- regressions will include a constant through the absorb() terms. Also note that when using -ivreg2- the constant will always be omitted because otherwise -ivreg2- will complain about the "orthogonality conditions S is not of full rank" (when using clustered standard errors).


## VCE

unadjusted = conventional =? ols
robust
cluster


### Accepted choices within vce()

unadjusted
robust
cluster varnames // -ivreg2- supports up to two vars, -ivregress- 1, and -regress- unlimited (but be sane about it)
[any of the above], bw(#) kernel(str) // This doesn't work with -ivregress-

Note that the last one can replicate dkraay and kiefer:
dkraay(#) = cluster(tsvar), bw(#) kernel(bartlett)
kiefer = kernel(tru) bw(full time length of panel)

Also, ivsuite can be -avar- or -default-
ivreg2 will ignore ivsuite and use it's own (which is based on avar)
ivregress will also do its own thing
regress will do its own thing but..
- if there are 3+ clustervars it will use our custom program to do mwc
- if there are 2 it can choose from -avar- and custom, choosing the later by default
- for the advanced stuff (bw, kernel), use -avar-
- for everything else, do its own thing



We need to preserve the timevar (and thus the panelvar) if we are using bw()