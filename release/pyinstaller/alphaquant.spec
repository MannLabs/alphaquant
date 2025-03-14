# -*- mode: python ; coding: utf-8 -*-

import os
import sys
from PyInstaller.building.build_main import Analysis, PYZ, EXE, COLLECT, BUNDLE, TOC
import PyInstaller.utils.hooks


##################### User definitions
exe_name = 'alphaquant_gui'
script_name = 'alphaquant_pyinstaller.py'
if sys.platform[:6] == "darwin":
	icon = '../logos/alpha_logo.icns'
else:
	icon = '../logos/alpha_logo.ico'
block_cipher = None
location = os.getcwd()
project = "alphaquant"
bundle_name = "alphaquant"
#####################


datas, binaries, hidden_imports = PyInstaller.utils.hooks.collect_all(
	project,
	include_py_files=True
)

# add extra packages that don't have pyinstaller hooks
# extra_pkgs = ["alphabase", ] # other alphaX packages would be added here
# for pkg in extra_pkgs:
# 	_datas, _binaries, _hidden_imports = PyInstaller.utils.hooks.collect_all(
# 		pkg,
# 		include_py_files=True
# 	)
# 	datas+=_datas
# 	binaries+=_binaries
# 	hidden_imports+=_hidden_imports

# prepare hidden imports and datas
hidden_imports = [h for h in hidden_imports if "__pycache__" not in h]
# hidden_imports = sorted(
# 		[h for h in hidden_imports if "tests" not in h.split(".")]
# 	)
datas = [d for d in datas if ("__pycache__" not in d[0]) and (d[1] not in [".", "Resources", "scripts"])]

# add certifi to datas, otherwise ssh connections fail when they are triggered from the installer, because the certificates are not available
# In the case of the AlphaQuant repo, AlphaMap needs to download data from GitHub and this fails without certifi
datas.extend(PyInstaller.utils.hooks.collect_data_files('certifi'))

# add matplotlib backends to hidden imports
# When using the GUI with windows installer, runs fail because these matplotlib backends are not available. No issues when
# running from command line on windows. And no issues on macOS.
hidden_imports.extend([
	'matplotlib.backends.backend_pdf',
	'matplotlib.backends.backend_agg'
])

a = Analysis(
	[script_name],
	pathex=[location],
	binaries=binaries,
	datas=datas,
	hiddenimports=hidden_imports,
	hookspath=[],
	runtime_hooks=[],
	excludes=[h for h in hidden_imports if "datashader" in h],
	win_no_prefer_redirects=False,
	win_private_assemblies=False,
	cipher=block_cipher,
	noarchive=False
)
pyz = PYZ(
	a.pure,
	a.zipped_data,
	cipher=block_cipher
)

if sys.platform[:5] == "linux":
	exe = EXE(
		pyz,
		a.scripts,
		a.binaries,
		a.zipfiles,
		a.datas,
		name=bundle_name,
		debug=False,
		bootloader_ignore_signals=False,
		strip=False,
		upx=True,
		console=True,
		upx_exclude=[],
		icon=icon
	)
else: # non-linux
	exe = EXE(
		pyz,
		a.scripts,
		# a.binaries,
		a.zipfiles,
		# a.datas,
		exclude_binaries=True,
		name=exe_name,
		debug=False,
		bootloader_ignore_signals=False,
		strip=False,
		upx=True,
		console=True,
		icon=icon
	)
	coll = COLLECT(
		exe,
		a.binaries,
		# a.zipfiles,
		a.datas,
		strip=False,
		upx=True,
		upx_exclude=[],
		name=exe_name
	)
