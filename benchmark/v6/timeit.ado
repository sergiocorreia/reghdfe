capture program drop timeit
program define timeit

	_on_colon_parse `0'
	loc 0 `s(before)'
	loc cmd `"`s(after)'"'
	args timer

	_assert inrange(`timer', 1, 100)

	timer on `timer'
	di as text "`timer' " _c
	qui `cmd'
	timer off `timer'

	qui timer list `timer'
	loc prev_reps = r(nt`timer') - 1
	global delta1_`timer' = (`prev_reps' * ${delta1_`timer'} + abs($b1 - _b[x1])) / (`prev_reps' + 1)
end
