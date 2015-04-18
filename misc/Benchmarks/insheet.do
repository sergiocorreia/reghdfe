* foreach fn in 2e6 1e7 1e8 	{
foreach fn in 1e5 {
	insheet using `fn'.csv, tab clear
	egen g1=group(id1)
	egen g2=group(id2)
	compress
	keep v3 v2 id4 id5 id6 g1 g2
	save `fn', replace
	clear
}
