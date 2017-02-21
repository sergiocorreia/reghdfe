global initial_obs .
global mufasa 0
global megaseeds

loc try 1
while (1) {
	di "."
	set seed 30462`try'
	loc ++try
	global megaseed = int(runiform()*100000)

	qui do quadros_demo

	di as text "try=`try'"
	di as text "megaseed=$megaseed"
	di as text "stopit=$mufasa"
	di as text "obs=`c(N)'"

	if ($mufasa) {
		global megaseeds "$megaseeds $megaseed"
		global mufasa = 0
		global initial_obs = c(N)
		di as text "MEGASEEDS UPDATED!"
		di "$megaseeds"
	}
	else {
		//cls
	}
}

exit
