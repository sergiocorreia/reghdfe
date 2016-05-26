noi cscript "reghdfe postestimation: predict" adofile reghdfe

* Setup
	discard
	clear all
	set more off

* Create dataset
	sysuse auto
	egen tt =  group(turn trunk)
	gen c = 1
	gen d = trunk
	replace d = turn if d > 15

* Test #1
	reghdfe price, a(turn trunk tt) keepsing
	assert e(redundant) == 20

	reghdfe price, a(turn trunk tt) keepsing dof(two)
	assert e(redundant) == 3

	group3hdfe trunk turn tt
	assert r(rest) == 36

	reghdfe price, a(turn trunk tt) keepsing dof(three)
	assert e(redundant) == 36

* Test #2
	reghdfe price, a(tt turn trunk) keepsing
	assert e(redundant) == 36

	group3hdfe tt turn trunk
	assert r(rest) == 36

	reghdfe price, a(tt turn trunk) keepsing dof(three)
	assert e(redundant) == 36

* Test #3
	reghdfe price, a(c turn d) keepsing dof(two)
	assert e(redundant) == 2

	group3hdfe c turn d
	assert r(rest) == 7

	reghdfe price, a(c turn d) keepsing dof(three)
	assert e(redundant) == 7

cd "C:/Git/reghdfe/test"
exit
