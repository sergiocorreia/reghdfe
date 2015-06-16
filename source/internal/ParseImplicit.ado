capture program drop ParseImplicit
program define ParseImplicit
* Parse options in the form NAME|NAME(arguments)
	* opt()			name of the option (so if opt=spam, we can have spam or spam(...))
	* default()		default value for the implicit form (in case we don't have a parenthesis)
	* syntax()		syntax of the contents of the parenthesis
	* input()		text to parse (usually `options', the result of a previous syntax .. , .. [*] )
	* inject()		what locals to inject on the caller (depend on -syntax)
	* xor			opt is mandatory (one of the two versions must occur)
	syntax, opt(name local) default(string) syntax(string asis) [input(string asis)] inject(namelist local) [xor]

	* First see if the implicit version is possible
	local lower_opt = lower("`opt'")
	local 0 , `input'
	cap syntax, `opt' [*]
	if ("`xor'"=="") local capture capture
	local rc = _rc
	if (`rc') {
		`capture' syntax, `opt'(string asis) [*]
		if ("`capture'"!="" & _rc) exit
	}
	else {
		local `lower_opt' `default'
	}
	local 0 ``lower_opt''
	syntax `syntax'
	foreach loc of local inject {
		c_local `loc' ``loc''
	}
	c_local options `options'
end
