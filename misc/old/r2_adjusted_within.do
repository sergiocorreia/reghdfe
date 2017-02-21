cls
sysuse auto, clear
rename trunk id
*replace id = 1
bys id: gen t = _n
bys id: gen N = _N
tab N, m
drop if N==1
tsset id t


* local rhs // EMPTY (r2w=r2wa=0)
local rhs weight length

areg price `rhs', a(id)
di e(r2) // R2
di e(r2_a) // R2a


cap noi xtreg price, fe
local wss = e(rss)

cap noi xtreg price `rhs', fe
di e(r2) // R2w
di e(r2_a) // ??
di 1 - (e(rss)/(e(N)-e(df_m)-1)) / (`wss'/(e(N)-1))
di e(df_m)
di e(df_a)


cap noi ivreg2 price `rhs', partial(_cons) // Works with only one cat
cap noi xtivreg2 price `rhs', fe
di e(r2) // R2w
di e(r2_a) // ??
di e(N)
di e(rankxx)
di e(sdofminus)
di e(dofminus)
di 1-(1-e(r2))*(e(N))/(e(N)-e(rankxx)-e(dofminus)-e(sdofminus))


reghdfe price `rhs', a(id)
di e(r2) // R2
di e(r2_a) // R2a

di e(r2_within) // R2w
di e(r2_a_within) // R2wa

di e(df_m)
di e(df_a)
di e(df_r)

di 1 - (e(rss)/(e(N)-e(df_m)-e(df_a))) / (e(tss_within)/(e(N)-e(df_a))) // Replicate reghdfe
di 1 - (e(rss)/(e(N)-e(df_m)-e(df_a))) / (e(tss_within)/(e(N)-1)) // Replicate xtreg
di 1 - (e(rss)/(e(N)-e(df_m)-e(df_a))) / (e(tss_within)/(e(N))) // Replicate xtivreg2

exit
