#!/bin/bash
set -e -u

# Build the installer for Linux.
# This script must be run from the root of the repository.
rm -rf dist_pyinstaller build_pyinstaller build_pyinstaller_data

# Find the wheel file in dist directory
WHL_NAME=$(cd dist && ls ./*.whl && cd ..)
pip install "dist/${WHL_NAME}[stable,gui-stable,dask-stable]"

INSTALLER_DATA_DIR=build_pyinstaller_data release/common/download_installer_data.sh

# Creating the stand-alone pyinstaller folder
pyinstaller release/pyinstaller/alphaquant.spec --distpath dist_pyinstaller --workpath build_pyinstaller -y
