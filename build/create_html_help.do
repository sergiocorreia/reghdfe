foreach fn in reghdfe hdfe {
	cd ../source
	copy `fn'.sthlp ../docs/`fn'.smcl, replace
	* local ufn = upper("`fn'")
	cd ../docs
	* Note: log2html kinda sucks and the results are a bit ugly
	* EG: {hline} is replaced by a bunch of "-" instead of <hr>
	log2html `fn'.smcl, replace erase ///
		title("help for `fn'.ado") ///
		line(80) ///
		css(help.css)
}

* scheme(yellow) /// black white yellow // Ugly: blue

* This creates the website scorreia.com/reghdfe
* We still need to manually upload the other git repo
cd ..
shell mkdocs build --clean
shell rmdir ..\sergiocorreia.github.io\reghdfe /S /Q
shell move site ..\sergiocorreia.github.io\
shell rename ..\sergiocorreia.github.io\site reghdfe

cd build
exit
