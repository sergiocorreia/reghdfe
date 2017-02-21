noi cscript "reghdfe: consider fweights when removing singletons" adofile reghdfe

* Prevent bug regression of https://github.com/sergiocorreia/reghdfe/issues/14

* Dataset
	sysuse auto
	gen w = 10

	reghdfe price weight gear length [fw=w], a(turn)
	storedresults save benchmark e()
	
	reghdfe price weight gear length [fw=w], a(turn) keepsingletons v(-1)
	storedresults compare benchmark e(), exclude(macro: cmdline scalar: drop_singletons)
	storedresults drop benchmark

	reghdfe price weight gear length [aw=w], a(turn)
	assert e(N)==70

	reghdfe price weight gear length [aw=w], a(turn) keepsing v(-1)
	assert e(N)==74

exit
