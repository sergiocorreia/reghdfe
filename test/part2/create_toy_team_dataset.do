* ===========================================================================
* Create toy team dataset
* ===========================================================================

	clear all
	cls


// --------------------------------------------------------------------------
// Test dataset
// --------------------------------------------------------------------------
	* Simulate round-robin tournament of # teams over # years,
	* where each team plays # players from a pool of # available players

	set seed 1234
	loc num_teams 4
	loc num_years 5
	loc num_players 2
	loc num_available 5

	loc num_rounds = (`num_teams' - 1) * 2
	loc games_per_round = `num_teams' / 2
	loc games_per_year = `num_rounds' * `games_per_round'
	loc num_games = `num_years' * `games_per_year'

	set obs `num_teams'
	gen byte team_id0 = _n
	tempfile home
	save "`home'"
	rename team_id0 team_id1
	cross using "`home'"
	drop if team_id0 == team_id1
	sort team_id0 team_id1

	gen long game_id = _n
	expand `num_years'
	bys game_id: gen byte year = _n
	assert c(N) == `num_games'
	sort year game_id
	gisid year game_id

	reshape long team_id, i(year game_id) j(is_home)

	gegen long game_year_id = group(year game_id)
	gisid game_year_id is_home
	gisid game_year_id team_id

	expand `num_available'
	bys game_year_id team_id: gen byte _ = _n
	gegen long player_id = group(team_id _)
	drop _

	gen u = runiform()
	bys game_year_id team_id (u): keep if _n <= `num_players'
	drop u

	gen x1 = 0.3 * rnormal()
	gen x2 = 0.5 * runiform()
	gen year_fe = log(1+year) - sin(year)
	gen team_fe = log(1+team_id)
	gen player_fe = sin(player_id)
	bys team_id: replace team_fe = team_fe[1]
	bys game_year_id team_id: replace team_fe = team_fe[1]
	bys game_year_id: replace x1 = x1[1]
	bys game_year_id: replace x2 = x2[1]
	gen _ = cond(is_home, player_fe, -player_fe)
	gegen sum_player_fe = total(_), by(game_year_id)
	drop _
	gen y = 0.1 + 3 * x1 - 1.5 * x2 + 0.1 * team_fe + 0.5 * sum_player_fe + 1 * rnormal()
	drop team_fe player_fe year_fe sum_player_fe
	bys game_year_id: replace y = y[1]

	gegen byte home_id = max( is_home*team_id), by(year game_id)
	gegen byte away_id = max(!is_home*team_id), by(year game_id)

	preserve
		* Wide-data alternative

		tab player_id, gen(dummy)
		foreach var of varlist dummy* {
			replace `var' = -`var' if !is_home
		}
		sort game_year_id, stable
		gcollapse (first) y x1 x2 home_id away_id (sum) dummy*, by(game_year_id game_id year) fast
		areg y x1 x2 dummy*, a(year)
		save "toy-teams-wide", replace
	restore

	gen byte slope = cond(is_home, 1, -1)
	save "toy-teams-long", replace
