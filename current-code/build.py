"""
Build stata ado package
TODO: change README.md to "tests - not run" or "tests - failed" until the script gets run
"""


# ---------------------------
# Imports
# ---------------------------

import re
import shutil
from pathlib import Path
import zipfile


# ---------------------------
# Constants
# ---------------------------

regex_include = re.compile(r'^\s*include "(?P<filename>\w+.mata)"\s?,\s?adopath')
regex_pkg_drop = re.compile(r'^d f (?P<filename>\w+\.(ado|sthlp|mata|mlib|do))$')
regex_pkg_drop_mata = re.compile(r'^f (?P<filename>\w+\.mata)$')
regex_pkg_copy = re.compile(r'^f (?P<filename>\w+\.(ado|sthlp|mata|mlib|do))$')
regex_version = re.compile(r'^\*!\s+version\s+(?P<version>[\d.]+) ')

# ---------------------------
# Functions
# ---------------------------

def copy_mata(input_path, input_fn, output_fn):
	'''
	Copy all Mata "include" directives into a single file
	'''

	assert output_fn.parent.stem == 'src'
	lines = input_fn.open(encoding='utf8', mode='r').read().splitlines()
	output = []

	print(f'\nCombining Mata files')
	print(f' - Input:  "{input_fn}"')
	print(f' - Output: "{output_fn}"')
	for line in lines:
		if (match := regex_include.match(line)):
			fn = match.group('filename')
			print(f' - including {fn}')
			full_fn = input_path / fn
			assert full_fn.is_file()
			code = open(full_fn, encoding='utf8', mode='r').read().splitlines()
			output.extend(code)
		else:
			output.append(line)

	output = '\n'.join(output) + '\n'
	with output_fn.open(encoding='utf8', mode='w', newline='\n') as fh:
		fh.write(output)
		print(f' - {len(output):,} lines written\n')


def copy_pkg(project, input_path, output_path, historical_path=None):
	'''
	Copy all files listed in the .pkg file, as well as a modified .pkg file
	(excludes .mata files)
	'''

	historical_files = []
	if historical_path is not None:
		print(f'Adding historical files:')
		for fn in historical_path.glob('*.*'):
			print(f' - {fn}')
			destination_fn = output_path / fn.name
			historical_files.append(fn.name)
			shutil.copyfile(src=fn, dst=destination_fn)


	pkg_fn = input_path / f'{project}.pkg'
	assert pkg_fn.is_file()

	print(f'Processing {pkg_fn}')
	lines = pkg_fn.open(encoding='utf8', mode='r').read().splitlines()
	output = []

	for line in lines:
		
		if (match := regex_pkg_drop.match(line)):
			fn = match.group('filename')
			print(f' - dropping {fn}')
			continue

		if (match := regex_pkg_drop_mata.match(line)):
			fn = match.group('filename')
			if fn == f'{project}.mata':
				output.append(line)
			print(f' - dropping {fn} (already included)')
			continue

		if (match := regex_pkg_copy.match(line)):
			fn = match.group('filename')
			print(f' - copying {fn}')
			full_fn = input_path / fn
			destination_fn = output_path / fn
			shutil.copyfile(src=full_fn, dst=destination_fn)

		output.append(line)

	for fn in historical_files:
		output.append(f'f {fn}')

	output_fn = output_path / f'{project}.pkg'
	output = '\n'.join(output) + '\n'
	with output_fn.open(encoding='utf8', mode='w', newline='\n') as fh:
		fh.write(output)
		print(f' - {len(output):,} lines written into "{output_fn}"\n')


def clean_build_path(path):
	print(f'Removing contents of {path} folder')
	for fn in path.glob('*.*'):
		if fn.suffix in ('.toc', '.pkg', '.ado', '.do', '.mata', '.mlib', '.sthlp'):
			print(f' - {fn}')
			fn.unlink()


def create_zip(path, project, version):
	zip_fn = path / f'{project}-{version}.zip'
	print(f'Creating zipfile "{zip_fn}"')
	with zipfile.ZipFile(zip_fn, 'w', zipfile.ZIP_DEFLATED) as zh:
		for fn in path.glob('*.*'):
			if fn.suffix in ('.toc', '.pkg', '.ado', '.do', '.mata', '.mlib', '.sthlp'):
				arcname = f'{project}/{fn.name}'
				print(arcname)
				zh.write(fn, arcname)
			else:
				print('FILE IGNORED:', fn.name)


def get_version(fn):
	with fn.open(encoding='utf8', mode='r') as fh:
		first_line = fh.readline()

	if (match := regex_version.match(first_line)):
		version = match.group('version')
	else:
		print('Version number not found in line:')
		print(first_line)
		print(f'of file {fn}')
		raise Exception('Version not found')

	print(f'Version: {version}')
	return version


# ---------------------------
# Main
# ---------------------------

def main():
	project = 'reghdfe'

	src_path = Path(".").absolute()
	assert src_path.stem == 'current-code'
	build_path = src_path.parent / 'src'
	historical_path = src_path.parent / 'historical-code' # reghdfe v3 and v5

	clean_build_path(build_path)

	input_fn = src_path / f'{project}.mata'
	output_fn = build_path / f'{project}.mata'
	copy_mata(input_path=src_path, input_fn=input_fn, output_fn=output_fn)

	copy_pkg(project=project, input_path=src_path, output_path=build_path, historical_path=historical_path)

	print('Copying stata.toc')
	src = src_path / 'stata.toc'
	dst = build_path / 'stata.toc'
	shutil.copyfile(src=src, dst=dst)

	version = get_version(fn=src_path/f'{project}.ado')
	create_zip(path=build_path, project=project, version=version)

	print('\nDone!')


if __name__ == '__main__':
	main()
