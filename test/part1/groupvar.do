clear
input byte(i j) float y
1 1  -1.705456
1 1  1.0428531
1 2  -.1415196
2 2  .15889134
3 2  .33107585
4 3 -.14910758
1 2 -1.3036556
1 1   .3811844
4 3  -.5003815
4 4  -2.692032
5 5   1.670801
5 6  .14196578
end

gen k = 1
reghdfe y, a(i j) groupvar(mob1) keepsing
reghdfe y, a(i#k j) groupvar(mob2) keepsing
reghdfe y, a(i j) groupvar(benchmark) keepsing version(3)
assert mob1 == benchmark
assert mob2 == benchmark
tab1 mob* benchmark, m
drop mob1 mob2 benchmark

reghdfe y, a(i j) groupvar(mob1)
reghdfe y, a(i#k j) groupvar(mob2)
reghdfe y, a(i j) groupvar(benchmark) version(5)
assert mob1 == benchmark
assert mob2 == benchmark
tab1 mob* benchmark, m

exit
