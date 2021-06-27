* Test error reported by Noah

set trace off
clear all
cls
adopath + "C:\Git\groupreg\test"

use "C:\Users\Sergio\Downloads\nba_mini_bug.dta"

rename pt_diff 			y
rename home				x1
rename b2b				x2
rename win_3gm_pct		x3
rename avg_pt_diff		x4
rename opp_avg_pt_diff	x5

rename player_season	indiv_id
rename team_game_id		group_id

*gen byte c = 1
*groupreg y x1 x2 x3 x4 x5, a(c indiv_id) gr(group_id) in(indiv_id)
*asd

cap noi groupreg y x3, a(indiv_id) gr(group_id) in(indiv_id)
rename x3 x
drop x?
cap noi groupreg y x, a(indiv_id) gr(group_id) in(indiv_id)
save "C:\Git\groupreg\test\data\error-15feb2021", replace
cls
*/

set trace off
*drop if mi(x) // problem is with MVs in x
chaos_drop_obs, f("C:\Git\groupreg\test\data\error-15feb2021") prob(0.1) maxiter(100) minobs(10) ignore(2001): ///
	groupreg y x, a(indiv_id) gr(group_id) in(indiv_id)

replace x = x *3
rename group_id _group_id
rename indiv_id _indiv_id

gegen long group_id = group(_group_id)
gegen long indiv_id = group(_indiv_id)
drop _*
order group_id indiv_id
sort  group_id indiv_id
keep group_id indiv_id y x
recast double x
compress

dataex

