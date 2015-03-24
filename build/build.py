# -*- coding: utf-8 -*-

"""
build
~~~~~~~~~~~~~~~

Put together all files for reghdfe.ado and place them in the ../package folder

Note: Wrote in Python 2.7 but should work with Python 3.
"""

# -------------------------------------------------------------
# Imports
# -------------------------------------------------------------

from __future__ import print_function
from __future__ import division

import os, time, re, shutil, zipfile

# -------------------------------------------------------------
# Functions
# -------------------------------------------------------------
def zipdir(path, zip):
    for root, dirs, files in os.walk(path):
        for file in files:
            zip.write(os.path.join(root, file))

# -------------------------------------------------------------
# Main
# -------------------------------------------------------------

# Build html help files
#os.system('C:\Bin\Stata13\StataMP-64.exe /e do "create_html_help"')
#THIS IS NOT WORKING CORRECTLY.. MAYBE DO IT BY HAND??

# Misc
os.chdir(os.path.split(__file__)[0])
fn_mata = ur"_mata/reghdfe.mata"
server_path = u"../package"
source_path = u"../source"

# Update version number
with open(os.path.join(source_path, "version.txt"), 'rb') as fh:
    old_version = fh.read()
    regex_string = ur'^(\d+)\.(\d+)\.(\d+) \d+\w+\d+$'
    regex = re.search(regex_string, old_version)
    today = time.strftime("%d%b%Y").lower() # See http://strftime.net/
    new_version = '{}.{}.{} {}'.format(regex.group(1), regex.group(2), int(regex.group(3))+1, today)
    print("Version updated from [{}] to [{}]".format(old_version, new_version))

# Filenames
output_filenames = ["reghdfe.ado", "reghdfe_estat.ado", "reghdfe_p.ado", "reghdfe_footnote.ado", "hdfe.ado"]


for fn in output_filenames:
    print("parsing file <{}>".format(fn))
    full_fn = os.path.join(source_path, fn)
    data = open(full_fn, "rb").read()
    source_data = None

    # Add Mata
    if ("include _mata/reghdfe.mata" in data):
        mata_data = open(os.path.join(source_path, fn_mata), "rb").read()
        data = data.replace(u"\r\nclear mata", mata_data)
        data = data.replace(u"\r\ninclude _mata/reghdfe.mata", mata_data)

    # Add other includes
    includes = re.findall('^\s*include "([^"]+)"', data, re.MULTILINE)
    for include in includes:
        print("    parsing include <{}>".format(include))
        full_include = os.path.join(source_path, include)
        include_data = open(full_include, "rb").read()
        data = data.replace(u'include "{}"'.format(include), '\r\n' + include_data)

    # Remove cap drop
    capdrops = re.findall('\s^\s*cap[a-z]* pr[a-z]* drop [a-zA-Z0-9_]+\s*$', data, re.MULTILINE)
    for capdrop in capdrops:
        data = data.replace(capdrop, "\n")

    # Update version
    data = data.replace("VERSION_NUMBER", new_version)

    # Save
    new_fn = os.path.join(server_path, fn)
    with open(new_fn, 'wb') as new_fh:
        new_fh.write(data)

# Update hdfe/reghdfe.pkg
for pkgname in ["reghdfe.pkg", "hdfe.pkg"]:
    print("updating date in " + pkgname)
    full_pkg = os.path.join(source_path, pkgname)
    pkg = open(full_pkg, "rb").read()
    today = time.strftime("%Y%m%d")
    pkg = re.sub(ur'Distribution-Date: \d+', ur'Distribution-Date: ' + today, pkg)
    open(full_pkg, 'wb').write(pkg)
    shutil.copy(full_pkg, os.path.join(server_path, pkgname))

# Copy
print("Copying misc files...")
shutil.copy(os.path.join(source_path, u"reghdfe.sthlp"), os.path.join(server_path, u"reghdfe.sthlp"))
shutil.copy(os.path.join(source_path, u"hdfe.sthlp"), os.path.join(server_path, u"hdfe.sthlp"))
shutil.copy(os.path.join(source_path, u"stata.toc"), os.path.join(server_path, u"stata.toc"))

print("Building zip file")
zipf = zipfile.ZipFile('../misc/reghdfe.zip', 'w')
zipdir('../package/', zipf)
zipf.close()

# Update version file now that the deed is done
with open(os.path.join(source_path, "version.txt"), 'wb') as fh:
    fh.write(new_version)

print("Done!")
