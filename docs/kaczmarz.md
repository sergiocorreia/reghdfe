# Kaczmarz Algo

If we ever decide to implement Kaczmarz for solving the alphas, this is the sketch of the algorithm:

Given Ax=b, where A is NxM, we start with a guess p=0 (Mx1)
Then, we iterate the following until convergence of p:

For every obs. i,

```algorithm
s = b[i]
for g=1:G
	level = F[g].levels[i]
	if has_intercept[g]: s -= F[g].alphas[level, 1]
	// if has_slope: iterate for every slope; don't forget to add back stdev
next
s = s / denom[i]
for g=1:G
	F[g].alphas[level, 1] += s
	// also update slopes
next
```

Beforehand, we created denom s.t. denom[i] = G (if only intercepts)
or = num_intercepts + cvar[i]'cvar[i]
