// -------------------------------------------------------------
// Standard version of the HDFE command
// -------------------------------------------------------------
cap pr drop AlternativeCMD
program define AlternativeCMD

    di in ye _n "[Running standard alternative to the HDFE command]"
    if !e(alternative_ok) di as error "Note: the command will be approximate only (e.g. wrong SDs, ommited AvgE)"
    di as text "[CMD] " as result "`e(alternative_cmd)'"
    local rhs `e(indepvars)' `e(endogvars)' `e(AvgE_Ws)'
    `e(alternative_cmd)'
    
    if ("`rhs'"!="") {
        di in ye _n "(Joint test on non-FE regressors)"
        testparm `rhs'
    }

end
