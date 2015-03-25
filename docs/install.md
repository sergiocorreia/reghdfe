# Installation

How to install reghdfe.ado and hdfe.ado

---

There are two versions: the stable version, which contains an older-but-proven version, and the development version, which contains the latest updates and improvements.

## Stable Version

It can be downloaded directly from Stata with the `ssc` command:

```stata
ssc install reghdfe
```

## Development Version

Users with Stata 13 or newer can download it with:

```stata
cap ado uninstall reghdfe
net from https://raw.githubusercontent.com/sergiocorreia/reghdfe/master/package/
net install reghdfe
```

Other users need to install it manually:

1. Download the zip file: https://github.com/sergiocorreia/reghdfe/tree/master/misc/reghdfe.zip
2. Extract it to an empty folder such as "C:\Temp\"
3. Within Stata, type the following (replacing "C:\Temp\" with whatever folder you chose)
```stata
cap ado uninstall reghdfe
net from "C:\Temp\"
net install reghdfe
```

reghdfe is now installed, and you can delete the zip file and the temporary folder.

## Installing hdfe.ado

If you want to install the stable version, use:

```stata
ssc install hdfe
```

The developer version is hosted on Github:

```stata
cap ado uninstall hdfe
net from https://raw.githubusercontent.com/sergiocorreia/reghdfe/master/package/
net install hdfe
```
