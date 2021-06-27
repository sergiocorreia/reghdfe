cap pr drop tic
pr tic
syntax, [TIMER(integer 10)]
	timer clear `timer'
	timer on `timer'
end

exit
