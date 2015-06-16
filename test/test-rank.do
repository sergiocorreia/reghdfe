cscript "reghdfe without full rank" adofile reghdfe

* Setup
	discard
	clear all
	set more off
	* cls

* Convenience: "Trim <size>" will trim e(b) and e(V)
	capture program drop TrimMatrix
	program define TrimMatrix, eclass
	args size
		assert `size'>0
		matrix trim_b = e(b)
		matrix trim_V = e(V)
		matrix trim_b = trim_b[1, 1..`size']
		matrix trim_V = trim_V[1..`size',1..`size']
		ereturn matrix trim_b = trim_b
		ereturn matrix trim_V = trim_V
	end


input x y i j
0 1 1 1
1 1 1 2
1 2 2 2
1 3 3 3
1 3 3 2
end

* NOTE: SEE THE REPLY TO THE ISSUE ON GITHUB

areg y x ibn.j , a(i) // Variables are omitted
di e(df_a)
di e(df_r)
reghdfe y x, absorb(i j)  keepsingletons
di e(df_r)
cap noi reghdfe y x, absorb(i j)
assert c(rc)==2001
di e(df_r)
reghdfe y x, absorb(i j)  keepsingletons
di e(df_r)
assert e(df_r)==0

cd "D:/Github/reghdfe/test"
exit
