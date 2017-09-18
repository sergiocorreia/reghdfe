mata:

// --------------------------------------------------------------------------
// Transformations: Compute RESIDUALS, not projections
// --------------------------------------------------------------------------

`Void' function transform_cimmino(`FixedEffects' S, `Variables' y, `Variables' ans,| `Boolean' get_proj) {
	`Integer' 	g
	if (args()<4 | get_proj==.) get_proj = 0
	ans = S.project_one_fe(y, 1)
	for (g=2; g<=S.G; g++) {
		ans = ans + S.project_one_fe(y, g)
	}
	ans = get_proj ? ans / S.G : y - ans / S.G
}

// --------------------------------------------------------------------------

`Void' function transform_kaczmarz(`FixedEffects' S, `Variables' y, `Variables' ans,| `Boolean' get_proj) {
	`Integer' 	g

	if (args()<4 | get_proj==.) get_proj = 0

	ans = y - S.project_one_fe(y, 1)
	for (g=2; g<=S.G; g++) {
		ans = ans - S.project_one_fe(ans, g)
	}
	if (get_proj) ans = y - ans
}

// --------------------------------------------------------------------------
// This seems slower than kaczmarz (sym kaczmarz!); not used currently
`Void' function transform_rand_kaczmarz(`FixedEffects' S, `Variables' y, `Variables' ans,| `Boolean' get_proj) {
	`Integer' 	g
	`Vector' rand
	if (args()<4 | get_proj==.) get_proj = 0
	rand = sort( ( (1::S.G) , uniform(S.G,1) ) , 2 )[.,1]
	ans = y - S.project_one_fe(y, rand[1])
	for (g=2; g<=S.G; g++) {
		ans = ans - S.project_one_fe(ans, rand[g])
	}
	for (g=S.G-1; g>=1; g--) {
		ans = ans - S.project_one_fe(ans, rand[g])
	}
	if (get_proj) ans = y - ans
}

// --------------------------------------------------------------------------

 `Void' function transform_sym_kaczmarz(`FixedEffects' S, `Variables' y, `Variables' ans,| `Boolean' get_proj) {
	`Integer' 	g
	if (args()<4 | get_proj==.) get_proj = 0
	if (S.timeit) timer_on(72)
	ans = y - S.project_one_fe(y, 1)
	if (S.timeit) timer_off(72)
	for (g=2; g<=S.G; g++) {
		if (S.timeit) timer_on(72)
		ans = ans - S.project_one_fe(ans, g)
		if (S.timeit) timer_off(72)
	}
	for (g=S.G-1; g>=1; g--) {
		if (S.timeit) timer_on(72)
		ans = ans - S.project_one_fe(ans, g)
		if (S.timeit) timer_off(72)
	}
	if (get_proj) ans = y - ans
}
end
