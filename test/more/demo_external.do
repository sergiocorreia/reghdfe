discard
pr drop _all
clear all
sysuse auto

gen byte touse = (_n < 70)
gen fw = 1
mata: HDFE = fixed_effects("turn#trunk foreign", "touse", "fweight", "fw", 0, 1)

mata: HDFE.options.num_clusters = 2
mata: HDFE.options.clustervars = tokens("turn foreign")
mata: HDFE.options.base_clustervars = subinstr(HDFE.options.clustervars, "#", " ")

mata: HDFE.estimate_dof()
//mata: HDFE.save_touse()
mata: y = HDFE.partial_out("price weight")
mata: HDFE.output.post_footnote()
reghdfe_footnote


exit


mata:
	HDFE = fixed_effects("turn trunk")
	y = HDFE.partial_out("price weight")
end
