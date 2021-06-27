* ===========================================================================
* Test options specific to reghdfe with individual FEs
* ===========================================================================
	cd C:\Git\groupreg\test\
	*do setup
	cd C:\Git\groupreg\test\part2
	
	clear all
	* do create_toy_team_dataset

	use "./toy-teams-wide", clear
	reghdfe y x1 x2, a(home_id away_id)
	loc b1v1 = _b[x1]
	loc b2v1 = _b[x2]
	
	reghdfe y x1 x2 dummy*, a(home_id away_id)
	loc b1v2 = _b[x1]
	loc b2v2 = _b[x2]


	use "./toy-teams-long", clear
	
	reghdfe y x1 x2, a(home_id away_id) group(game_year_id)
	assert reldif(`b1v1', _b[x1]) < 1e-8
	assert reldif(`b2v1', _b[x2]) < 1e-8

	reghdfe y x1 x2, a(home_id away_id player_id#c.slope) indiv(player_id) group(game_year_id) aggreg(sum) precond(none)
	assert reldif(`b1v2', _b[x1]) < 1e-8
	assert reldif(`b2v2', _b[x2]) < 1e-8

	//reghdfe y x1 x2, a(home_id away_id player_id#c.slope) indiv(player_id) group(game_year_id) aggreg(sum) precond(diag)
	//assert reldif(`b1v2', _b[x1]) < 1e-8
	//assert reldif(`b2v2', _b[x2]) < 1e-8

exit






	use "./toy-teams-long", clear
	//reghdfe y x1 x2, a(team_id#is_home player_id#c.slope) indiv(player_id) group(game_year_id) aggreg(sum) precond(none)
	

exit
