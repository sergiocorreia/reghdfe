# Potential Pitfalls

(*WORK IN PROGRESS*)

You need to be careful when using a high number of fixed effects, as well as advanced standard error options, as it allows you to "shoot yourself in the foot "

Quick list of pitfalls that need to be discussed:

1. Ignore the constant; it doesn't tell you much. If you want to use descriptive stats, that's what the -summ- commands are for. Even better, use -noconstant- to drop it (although it's not really dropped as it never existed on the first place!)
2. You almost surely do not want tosave the fixed effects. They are probably inconsistent / not identified and you will likely use them wrong.
3. It's good practice to drop singletons. Thus, use -dropsingleton-.
4. If you use vce(robust), be sure that your *other* dimension is not "fixed" but grows with N, or your SEs will be wrong.
5. If you use vce(cluster ..), check that your number of clusters is high enough! If not, you are making the SEs even worse!
6. The panel variables (absvars) should probably be nested within the clusters (clustervars) due to the within-panel correlation induced by the FEs
(this is not the case for *all* the absvars, only those that are treated as growing as N grows)
