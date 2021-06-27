program timeit
	_on_colon_parse `0'
	loc n  = s(before)
	loc cmd `"`s(after)'"'
	timer on `n'
	di as input `"timer=[`n'] cmd=[`cmd']"'
	`cmd'
	timer off `n'
end
