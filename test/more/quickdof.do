* Setup
	discard
	clear all
	set more off

cap ado uninstall moresyntax
net install moresyntax, from("C:/git/moresyntax/src")

cap ado uninstall ftools
net install ftools, from("C:/git/ftools/src")
ftools, check
ftools, compile

cap ado uninstall fereg
net install fereg , from(C:/git/reghdfe/src)

fereg, check
mata: mata desc using lreghdfe

fereg, compile


* Create dataset
	sysuse auto
	egen tt =  group(turn trunk)
	gen c = 1
	gen d = trunk
	replace d = turn if d > 15

* Test #1
	fereg price, a(turn trunk tt) keepsing
	reghdfe price, a(turn trunk tt) keepsing
	//assert e(redundant) == 20

	fereg price, a(turn trunk tt) keepsing
	reghdfe price, a(turn trunk tt) keepsing // dof(two)
	//assert e(redundant) == 3

	group3hdfe trunk turn tt
	//assert r(rest) == 36

	fereg price, a(turn trunk tt) keepsing // dof(three)
	reghdfe price, a(turn trunk tt) keepsing // dof(three)
	//assert e(redundant) == 36

* Test #2
	fereg price, a(tt turn trunk) keepsing
	reghdfe price, a(tt turn trunk) keepsing
	//assert e(redundant) == 36

	group3hdfe tt turn trunk
	//assert r(rest) == 36

	fereg price, a(tt turn trunk) keepsing // dof(three)
	reghdfe price, a(tt turn trunk) keepsing // dof(three)
	//assert e(redundant) == 36

* Test #3
	fereg price, a(c turn d) keepsing // dof(two)
	reghdfe price, a(c turn d) keepsing // dof(two)
	//assert e(redundant) == 2

	group3hdfe c turn d
	//assert r(rest) == 7

	fereg price, a(c turn d) keepsing // dof(three)
	reghdfe price, a(c turn d) keepsing // dof(three)
	//assert e(redundant) == 7

exit
