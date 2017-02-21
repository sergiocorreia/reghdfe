cls
set more off
clear all
sysuse auto
drop if rep==.
clonevar fofo= foreign
gen asd = 10*rep+foreign

GenerateID foreign, replace clustervars(turn price)
char list

GenerateID trunk, replace clustervars(turn trunk)
char list


GenerateID fofo, replace clustervars(turn price foreign)
char list

GenerateID asd, replace clustervars(turn price rep)
char list
