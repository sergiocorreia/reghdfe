clear all
cls

sysuse auto, clear
local depvar 	"price"
local indepvars "weight gear"
local absvars 	"turn trunk"


reghdfe `depvar' `indepvars', a(`absvars')

qui include "reghdfe.mata", adopath
mata:
    // 0) Optional declaration
    // class FixedEffects scalar HDFE

    // 1) Create the object
    HDFE = FixedEffects() // Note that you can replace "HDFE" with whatever name you choose

    // 2) Set up options as needed
    HDFE.absvars = "`absvars'"

    // 3) Initialize (validate options)
    HDFE.init()

    // 4) Partial out
    HDFE.partial_out("`depvar' `indepvars'")


    // 5) Solve OLS
    data = HDFE.solution.data
    k = cols(data)
    y = data[., 1]
    x = data[., 2::k]
    b = qrsolve(x, y)

    // Note: we standardized variables when partialling out; need to undo this
    HDFE.solution.stdevs
    stdev_y = HDFE.solution.stdevs[1]
    stdev_x = HDFE.solution.stdevs[2..k]
    b :/ stdev_x' * stdev_y
end

exit
