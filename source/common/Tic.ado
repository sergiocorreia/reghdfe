capture program drop Tic
program define Tic
syntax, n(integer)
	timer clear `n'
	timer on `n'
end
