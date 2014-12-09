cap cd tests
clear all
include ../_mata/reghdfe.mata

mata
group = 1\1\2\2\1\1\3
index = order(group, 1)
freqs = 1\1\1\1\1\1\30
sorted_freqs = freqs[index, 1]
ans = count_by_group(group, index, 0, sorted_freqs)
assert(all(ans==(4\2\30)))
end

* ssc install moremata
mata
x = 1.5 \ 2.3 \ 6.5 \ 4.2 \ 1.1 \ 4.6 \ 0.0
group = 1\1\2\2\1\1\3
freqs = 1\3\1\2\1\10\30

index = order(group, 1)
sorted_freqs = freqs[index, 1]

counts_to = count_by_group(group, index)
sum_count = quadrunningsum(counts_to)

counts_to = count_by_group(group, index, 0 , sorted_freqs)

// x , group , freqs
ans = mean_by_group(x, index, sum_count, counts_to, sorted_freqs)
//check = 3.915385 \ 5.35 \ 0
check = 3.7 \ 4.96666 \ 0
ans , check

assert(sum(abs(ans-check))<0.01)


end



exit
DEBUG:

discard
clear all
cls
cd D:\Dropbox\Projects\stata\hdfe\code
sysuse auto
reghdfe price weight, a(turn##c.length turn#c.head)
assert abs(_b[weight]-5.759)<0.1
areg price weight i.turn#c.length   i.turn#c.head, a(turn)




* ESTA MAL EL REGHDFE!









cd D:\Dropbox\Projects\stata\hdfe\code
sysuse auto
reghdfe price weight, a(turn rep)
