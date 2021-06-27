
clear all

include ../src/groupreg.mata

cls

clear
input int(y x) byte(patent inventor)
1582 1960 1 1
1582 1960 1 2
1582 1960 1 3
 867  961 2 4
5106 6341 3 1
5106 6341 3 2
5106 6341 3 3
7135 5468 4 5
7135 5468 4 1
7135 5468 4 3
7135 5468 4 4
3253  235 5 2
3253  235 5 3
4235 3255 6 6
 323  255 6 1
 323  205 7 7
end
list
sort patent inventor, stable

gen byte indiv_touse = 1
by patent (inventor): gen byte touse = _n == 1

mata:

GF = fe_factor("patent")
IF = fe_factor("inventor")
bg = BipartiteGraph()

bg.init(&GF, &IF, 5)
(void) bg.init_zigzag(1) // 1 => save subgraphs into bg.subgraph_id
bg.compute_cores()

weight = 1
bg.prune_1core(weight)
num_pruned = bg.N_drop

//GF.levels, IF.levels, bg.mask


sample = selectindex(st_data(., "touse"))
indiv_sample = selectindex(st_data(., "indiv_touse"))


sample
indiv_sample

"TESTING INDIV FACTOR"
F = indiv_factor("patent", "inventor", sample, indiv_sample, "sum", 3)

end

exit
