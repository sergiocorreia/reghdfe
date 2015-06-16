capture program drop Toc
program define Toc
syntax, n(integer) msg(string)
	timer off `n'
	qui timer list `n'
	di as text "[timer]{tab}" as result %8.3f `r(t`n')' as text "{col 20}`msg'{col 77}`n'" 
	timer clear `n'
end
