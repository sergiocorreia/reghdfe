* Unchecked assumptions:
* 1) groupid should be of type byte/int/long and take values from ZERO (0) to G (number of grups)
*	(This is what egen group generates)
* 2) meanvar should already exist and be of type float/real

version 12 // Maybe lower works

mata:
mata set matastrict on

function fastmean(string scalar varname, string scalar meanvar, string scalar groupid)
{
	real colvector Var, Newvar, ID // Nx1
	real colvector Mean, Count // Gx1
	real scalar G, N, i, k
	
	// st_data and st_view seem equally fast but st_view saves memory
	// Var = st_data(., varname)
	// ID = st_data(., groupid)	
	st_view(Var=., ., varname)
	st_view(ID=., ., groupid)
	
	G = max(ID)
	N = rows(Var)
	Mean = J(G,1, 0)
	Count = J(G,1, 0)
	Newvar = J(N,1, .)
	
	assert(G<N)
	
	for (i=1; i<=N; i++) {
		k = ID[i]
		Count[k] = Count[k] + 1 // Kinda inneficient but else I get error
		
		// http://www.stata-journal.com/sjpdf.html?articlenum=pr0025 - p557
		// Not entirely sure that this improves precision (remove final division if we uncomment this)
		Mean[k] = Mean[k] + (Var[i]) // - Mean[k]) / Count[k]
	}
	Mean = Mean :/ Count
	
	for (i=1; i<=N; i++) {
		k = ID[i]
		Newvar[i] = Mean[k]
	}
	st_store(., meanvar, Newvar)
}

end
