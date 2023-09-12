* ===========================================================================
* Test that all upstream/downstream reghdfe dependencies work
* ===========================================================================
* Note: this changes PLUS so best to restart Stata afterwards

	clear all


// --------------------------------------------------------------------------
// Test SSC
// --------------------------------------------------------------------------

	tempfile adopath
	loc adopath "C:\Git\asd\fakeado"
	di as text `"tempfile=`adopath'"'
	*mkdir "`adopath'", public

	sysdir set PLUS "`adopath'"

	net install require, from("C:\Git\stata-require\src") // SSC version is still 1.1.1
	require moremata				, install
	require require   >= 1.3		, install
	require ftools    >= 2.49.1		, install
	require reghdfe   >= 6.12.3		, install
	require ppmlhdfe  >= 2.3.0		, install
	*require ivreghdfe >= 1.1.3		, install // problem: SSC still has ivreghdfe 1.0.0

	sysuse auto
	fcollapse (sum) price, by(turn)

	sysuse auto, clear
	reghdfe price weight length, a(turn)
	*ivreghdfe price (weight=length), a(turn)

	set trace on
	set tracedepth 1
	ppmlhdfe turn weight length, a(trunk)

exit
