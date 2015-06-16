cap pr drop toc
pr toc, rclass
syntax, [TIMER(integer 10)] [REPORT]
	timer off `timer'
	qui timer list `timer'
	return local time = r(t`timer')
	if ("`report'"!="") {
		di as text "Done! (`c(current_time)', " string(r(t`timer'),"%20.1f") " seconds elapsed)" _n
	}
	timer clear `timer'
end

exit
