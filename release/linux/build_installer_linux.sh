#!/bin/bash
set -e -u

# Build the installer for Linux.
# This script must be run from the root of the repository.

rm -rf dist build *.egg-info
rm -rf dist_pyinstaller build_pyinstaller

# Creating the wheel
python -m build

# substitute X.Y.Z-devN with X.Y.Z.devN
WHL_NAME=$(echo "alphaquant-0.1.1-py3-none-any.whl" | sed 's/\([0-9][0-9]*\.[0-9][0-9]*\.[0-9][0-9]*\)-dev\([0-9][0-9]*\)/\1.dev\2/g')
pip install "dist/${WHL_NAME}[stable,gui-stable]"

# Creating the stand-alone pyinstaller folder
pip install pyinstaller
pyinstaller release/pyinstaller/alphaquant.spec --distpath dist_pyinstaller --workpath build_pyinstaller -y

