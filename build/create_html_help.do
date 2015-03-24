foreach fn in reghdfe hdfe {
	cd ../source
	copy `fn'.sthlp ../docs/`fn'.smcl, replace
	* local ufn = upper("`fn'")
	cd ../docs
	log2html `fn'.smcl, replace erase ///
		title("help for `fn'.ado") ///
		line(80) ///
		css(help.css)
}

* scheme(yellow) /// black white yellow // Ugly: blue


cd ..
shell mkdocs build --clean
cd build
		
exit
