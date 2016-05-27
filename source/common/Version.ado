// -------------------------------------------------------------
// Report HDFE/REGHDFE version
// -------------------------------------------------------------
cap pr drop Version
program define Version, eclass
    syntax , [STABLE DEV DEPENDENCIES]
    local all_dependencies ivreg2 avar tuples group3hdfe

    if ("`stable'" != "") {
        ado uninstall reghdfe
        ssc install reghdfe
        pr drop _all
        exit
    }

    if ("`dev'" != "") {
        ado uninstall reghdfe
        net install reghdfe, from("http://scorreia.com/software/reghdfe")
        pr drop _all
        exit
    }


    if ("`dependencies'" != "") {
        foreach dep of local all_dependencies {
            cap which `dep'
            if (_rc) ssc install `dep'
        }
        exit
    }

    local version "VERSION_NUMBER"
    ereturn clear
    di as text "`version'"
    ereturn local version "`version'"

    di as text _n "Dependencies installed?"
    foreach dep of local all_dependencies {
    	cap findfile `dep'.ado
    	if (_rc) {
    		di as text "{lalign 20: - `dep'}" as result " no" ///
             as text " {stata ssc install `dep':(click to install)}"
    	}
    	else {
    		di as text "{lalign 20: - `dep'}" as result "yes"
    	}
    }

    di as text _n "Updates:"
    di as text " - reghdfe:{stata reghdfe, version stable: update to latest stable version (from ssc)}"
    di as text " - reghdfe:{stata reghdfe, version dev: update to latest development version (from github)}"
    di as text " - dependencies:{stata reghdfe, version dependencies: install all}"
    di as text " - dependencies:{stata adoupdate update `all_dependencies', update: update all if installed}"

end
