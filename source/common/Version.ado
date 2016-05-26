// -------------------------------------------------------------
// Report HDFE/REGHDFE version
// -------------------------------------------------------------
cap pr drop Version
program define Version, eclass
    local version "VERSION_NUMBER"
    ereturn clear
    di as text "`version'"
    ereturn local version "`version'"

    di as text _n "Dependencies installed?"
    local dependencies ivreg2 avar tuples group3hdfe
    foreach dependency of local dependencies {
    	cap findfile `dependency'.ado
    	if (_rc) {
    		di as text "{lalign 20:- `dependency'}" as error "not"
    	}
    	else {
    		di as text "{lalign 20:- `dependency'}" as result "yes"
    	}
    }

end
