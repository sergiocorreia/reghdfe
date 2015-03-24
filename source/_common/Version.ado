// -------------------------------------------------------------
// Report HDFE/REGHDFE version
// -------------------------------------------------------------
cap pr drop Version
program define Version, eclass
    local version "VERSION_NUMBER"
    ereturn clear
    di as text "`version'"
    ereturn local version "`version'"

end
