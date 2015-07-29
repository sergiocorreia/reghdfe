# -*- coding: utf-8 -*-

"""
build
~~~~~~~~~~~~~~~

Put together all files for reghdfe.ado and place them in the ../package folder

Requires Python 3.x
"""

# -------------------------------------------------------------
# Imports
# -------------------------------------------------------------
import os, time, re, shutil, zipfile, glob

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
fn_mata = r"map.mata"
server_path = "../package"
source_path = "../source"

header = """*! reghdfe {}
*! Sergio Correia (sergio.correia@duke.edu)

"""

# Clear -package- folder
files = glob.glob(server_path + '/*')
for f in files:
    os.remove(f)


# Update version number
with open(os.path.join(source_path, "version.txt"), 'rb') as fh:
    old_version = fh.read().decode()
    regex_string = r'^(\d+)\.(\d+)\.(\d+) \d+\w+\d+\s*$'
    regex = re.search(regex_string, old_version)
    today = time.strftime("%d%b%Y").lower() # See http://strftime.net/
    new_version = '{}.{}.{} {}'.format(regex.group(1), regex.group(2), int(regex.group(3))+1, today)
    print("Version updated from [{}] to [{}]".format(old_version, new_version))

# Append Mata includes
print("parsing map.mata")
full_fn = os.path.join(source_path, "mata", fn_mata)
mata_data = "\r" + open(full_fn, "rb").read().decode()
includes = re.findall('^\s*include\s+(\S+).mata', mata_data, re.MULTILINE)
for i, include in enumerate(includes,1):
 
    print("    parsing mata include <{}>".format(include), end="")
    full_include = os.path.join(source_path, "mata", include + ".mata")
    include_data = open(full_include, "rb").read().decode()
    print(": {} lines".format(len(include_data.split('\n'))))
    
    mata = re.findall('^\s*mata:\s*$', include_data, re.MULTILINE)
    if (len(mata)>1): print("mata: appears more than once")
    assert len(mata)==1
    if (i>1): include_data = include_data.replace(mata[0],"")

    stricts = re.findall('^\s*mata set matastrict on\s*$', include_data, re.MULTILINE)
    if (len(stricts)>1): print("matastrict appears more than once")
    assert len(stricts)==1
    if (i>1): include_data = include_data.replace(stricts[0],"")

    ends = re.findall('^\s*end\s*$', include_data, re.MULTILINE)
    if (len(ends)>1): print("end appears more than once")
    assert len(ends)==1
    if (i<len(includes)): include_data = include_data.replace(ends[0],"")
    mata_data = mata_data.replace('include {}.mata'.format(include), '\r\n' + include_data.strip())

# Filenames
output_filenames = ["reghdfe.ado", "reghdfe_estat.ado", "reghdfe_p.ado", "reghdfe_footnote.ado", "hdfe.ado"]

for fn in output_filenames:
    print("parsing file <{}>".format(fn))
    full_fn = os.path.join(source_path, fn)
    data = open(full_fn, "rb").read().decode()
    source_data = None

    # Add Mata
    if ('include "mata/map.mata"' in data):
        data = data.replace("\r\nclear mata", "\r\n")
        data = data.replace('\r\ninclude "mata/map.mata"', mata_data)

    # Add other includes
    includes = re.findall('^\s*include "([^"]+)"', data, re.MULTILINE)
    for include in includes:
        print("    parsing include <{}>".format(include), end="")
        full_include = os.path.join(source_path, include)
        include_data = open(full_include, "rb").read().decode()
        data = data.replace('include "{}"'.format(include), '\r\n' + include_data)
        print(": {} lines".format(len(include_data.split('\n'))))

    # Remove cap drop
    capdrops = re.findall('\s^\s*cap[a-z]* pr[a-z]* drop [a-zA-Z0-9_]+\s*$', data, re.MULTILINE)
    for capdrop in capdrops:
        data = data.replace(capdrop, "\n")

    # Update version
    if fn in ("reghdfe.ado", "hdfe.ado"):
        data = header.format(new_version) + data
    data = data.replace("VERSION_NUMBER", new_version)

    # Save
    new_fn = os.path.join(server_path, fn)
    with open(new_fn, 'wb') as new_fh:
        new_fh.write(data.encode())

# Update hdfe/reghdfe.pkg
for pkgname in ["reghdfe.pkg", "hdfe.pkg"]:
    print("updating date in " + pkgname)
    full_pkg = os.path.join(source_path, pkgname)
    pkg = open(full_pkg, "rb").read().decode()
    today = time.strftime("%Y%m%d")
    pkg = re.sub(r'Distribution-Date: \d+', r'Distribution-Date: ' + today, pkg)
    open(full_pkg, 'wb').write(pkg.encode())
    shutil.copy(full_pkg, os.path.join(server_path, pkgname))

# Copy
print("Copying misc files...")
shutil.copy(os.path.join(source_path, "reghdfe.sthlp"), os.path.join(server_path, "reghdfe.sthlp"))
shutil.copy(os.path.join(source_path, "hdfe.sthlp"), os.path.join(server_path, "hdfe.sthlp"))
shutil.copy(os.path.join(source_path, "stata.toc"), os.path.join(server_path, "stata.toc"))
shutil.copy(os.path.join(source_path, "estfe.ado"), os.path.join(server_path, "estfe.ado"))

shutil.copy(os.path.join(source_path, "reghdfe_old.sthlp"), os.path.join(server_path, "reghdfe_old.sthlp"))
shutil.copy(os.path.join(source_path, "reghdfe_old.ado"), os.path.join(server_path, "reghdfe_old.ado"))
shutil.copy(os.path.join(source_path, "reghdfe_old_p.ado"), os.path.join(server_path, "reghdfe_old_p.ado"))
shutil.copy(os.path.join(source_path, "reghdfe_old_estat.ado"), os.path.join(server_path, "reghdfe_old_estat.ado"))
shutil.copy(os.path.join(source_path, "reghdfe_old_footnote.ado"), os.path.join(server_path, "reghdfe_old_footnote.ado"))

# print("Building zip file")
# zipf = zipfile.ZipFile('../misc/reghdfe.zip', 'w', zipfile.ZIP_DEFLATED)
# zipdir('../package/', zipf)
# zipf.close()

# Update version file now that the deed is done
with open(os.path.join(source_path, "version.txt"), 'wb') as fh:
    fh.write(new_version.encode())

print("Done!")
