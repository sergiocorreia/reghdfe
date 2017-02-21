clear

scalar N        = 11
scalar df_r     = 10
scalar df_m     = 1

scalar rmse     = 3.56753034006338
scalar r2       = 0.999365492298663
scalar mss      = 200457.727272727
scalar F        = 15750.2500000000
scalar rss      = 127.272727272727

scalar bx       = 2.07438016528926
scalar sex      = 0.165289256198347E-01

qui input int (y x)
         130    60
         131    61
         132    62
         133    63
         134    64
         135    65
         136    66
         137    67
         138    68
         139    69
         140    70
end

reg y x, nocons
di "R-squared = " %20.15f e(r2)

assert N    == e(N)
assert df_r == e(df_r)
assert df_m == e(df_m)

lrecomp _b[x] bx () _se[x] sex () /*
*/ e(rmse) rmse e(r2) r2 e(mss) mss e(F) F e(rss) rss
