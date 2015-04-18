capture program drop SortPreserve
program define SortPreserve, sortpreserve
	_on_colon_parse `0'
	`s(after)'
end
