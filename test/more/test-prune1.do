
cls


cap ado uninstall ftools
net install ftools, from("C:/git/ftools/src")
ftools, compile

cap ado uninstall fereg
net install fereg, from("C:/git/fereg/src")
fereg, compile

use c:\git\fereg\test\prune, clear

fereg y x, a(id1 id2) v(3)

exit
reghdfe y x, a(id1 id2)
fereg y x, a(id1 id2) noprune
areg y ibn.id2 x, a(id1)
