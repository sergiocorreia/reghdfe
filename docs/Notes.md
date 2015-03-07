We need a pitfalls section. For instance:
1. Ignore the constant; it doesn't tell you much. If you want to use descriptive stats, that's what the -summ- commands are for. Even better, use -noconstant- to drop it (although it's not really dropped as it never existed on the first place!)
2. Do not save the fixed effects. They are probably inconsistent / not identified and you will likely use them wrong.
3. It's good practice to drop singletons. Thus, use -dropsingleton-.
4. If you use vce(robust), hope that your *other* dimension is not fixed, or your SEs will be wrong.
5. If you use vce(cluster ..), check that your number of clusters is high enough! If not, you are making the SEs even worse!
6. Probably also the panel variables (absvars) should be nested within the clusters (clustervars) due to the within-panel correlation induced by the FEs
(this is not the case for *all* the absvars, only those that are treated as growing as N grows)


Remember that in a FE panel regression, using a HC (het-consistent aka white/sandwich) estimate for the variance will 
lead to inconsistent results for the VCV  if the time dimension is fixed (!!!)

The soln in that case is to use clustered errors.

See:
http://www.princeton.edu/~mwatson/papers/ecta6489.pdf
Stock and Watson, "Heteroskedasticity-robust standard errors for fixed-effects panel-data regression," Econometrica 76 (2008): 155-174


This translates to the reghdfe case as follows:
If you use vce(robust), then there should be enough observations within each firm/group/individual (or whatever unit the absvars point to).
Again, we need to be able to use asymptotics within *all* of absvars. This is a *huge* assumption,
so most of the time our best bet is to use vce(cluster clustervars); taking care that again there are enough clusters for each clustervar



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