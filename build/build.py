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

import os, time, re, shutil


# Filenames
output_filenames = [ur"reghdfe.ado", ur"reghdfe_absorb.ado", 
    ur"reghdfe_estat.ado", ur"reghdfe_p.ado", ur"reghdfe_footnote.ado"]

fn_mata = ur"_mata/reghdfe.mata"
server_path = u"../package"
source_path = u"../source"

for fn in output_filenames:
    print("parsing file <{}>".format(fn))
    full_fn = os.path.join(source_path, fn)
    data = open(full_fn, "rb").read()

    # Change header
    if (fn==ur"reghdfe.ado"):
        regex = re.search(ur'^\*! reghdfe (\d+)\.(\d+)\.(\d+) \d+\w+\d+', data)
        version = '{}.{}.{}'.format(regex.group(1), regex.group(2), int(regex.group(3))+1)
        today = time.strftime("%d%b%Y").lower() # See http://strftime.net/
        header = '*! version {} {}'.format(version, today)
        data = data.replace(regex.group(0), header)

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
    capdrops = re.findall('\s^\s*cap pr drop [a-zA-Z0-9_]+\s*$', data, re.MULTILINE)
    for capdrop in capdrops:
        data = data.replace(capdrop, "")        

    # Save
    new_fn = os.path.join(server_path, fn)
    with open(new_fn, 'wb') as new_fh:
        new_fh.write(data)

# Update reghdfe.pkg
print("updating date in reghdfe.pkg")
full_pkg = os.path.join(source_path, u"reghdfe.pkg")
pkg = open(full_pkg, "rb").read()
today = time.strftime("%Y%m%d")
pkg = re.sub(ur'Distribution-Date: \d+', ur'Distribution-Date: ' + today, pkg)
open(full_pkg, 'wb').write(pkg)
shutil.copy(full_pkg, os.path.join(server_path, u"reghdfe.pkg"))

# Copy
print("Copying misc files...")
shutil.copy(os.path.join(source_path, u"reghdfe.sthlp"), os.path.join(server_path, u"reghdfe.sthlp"))
shutil.copy(os.path.join(source_path, u"stata.toc"), os.path.join(server_path, u"stata.toc"))

print("Done!")
