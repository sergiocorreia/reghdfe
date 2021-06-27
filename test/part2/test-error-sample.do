* ===========================================================================
* In this example having a sample creates problems when dropping singletons
* ===========================================================================

clear
input byte(group_id indiv_id y x)
1 1 -15 5
1 2 -15 5
1 4 -15 5
2 1  -2 1
2 2  -2 1
2 4  -2 1
3 1   8 2
3 2   8 2
3 3   8 2
3 4   8 2
end


cap reghdfe y x if _n > 1, a(indiv_id) gr(group_id) in(indiv_id) v(3)
assert c(rc) == 2001
