program define chaos_drop_obs

	_on_colon_parse `0'
	loc cmd `"`s(after)'"'
	loc 0 `s(before)'

	syntax, Filename(string) PROB(real) [maxiter(integer 1000) minobs(integer 10) IGNORE(string asis) ALLOW(string asis)]
	assert inrange(`prob', 1e-6, 1-1e-6)

	*di as error `"BEFORE: `0'"'
	*di as error `"AFTER: `cmd'"'
	

	clear
	use "`filename'"
	cap `cmd'
	_assert c(rc), msg(`"Command ran without error: `cmd'"')

	forval i = 1/`maxiter' {
		qui use "`filename'", clear
		loc n = c(N)
		qui drop if runiform() < `prob'
		cap `cmd'
		if (c(rc) & (c(N)<`n')) {

			if ("`ignore'" != "") {
				if inlist(c(rc), `ignore') {
					di as text "x"
					continue
				}
			}

			if ("`allow'" != "") {
				if !inlist(c(rc), `allow') {
					di as text "x"
					continue
				}
			}
			
			di as text "obs=`c(N)'"
			qui save "`filename'", replace
			if (c(N) < `minobs') {
				continue, break				
			}
		}
		else {
			di as text "."
		}
	}

end
