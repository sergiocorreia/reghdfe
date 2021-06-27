* ===========================================================================
* In this example there are no observations left after dropping indiv FE singletons
* ===========================================================================

clear
input long group_id int(indiv_id season) byte(y x opp_team_id)
 4789 198 2018  17 1  8
 6025 137 2019  -7 1 16
 7547  61 2019  -4 1 14
 8887  95 2017  -8 . 11
11513 145 2017  -2 1 28
11863  94 2018 -12 .  9
13557 302 2017 -21 2 10
40489  72 2018 -11 1 27
46459 222 2018   8 0 30
end

drop if mi(x)

cap reghdfe y x, a(indiv_id) gr(group_id) in(indiv_id) v(3)
assert c(rc) == 2001
