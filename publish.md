Instructions to publish update on website:

```
cd c:\git\parse-smcl
smcl2html.py C:\git\reghdfe\src\reghdfe.sthlp --adopath=C:\Bin\Stata14\ado\base --view
```

Then copy new file and add the following at the top:

```
---
title: Help File for REHGDFE.ADO
layout: page
sub-menu: reghdfe
is-help: true
---
```