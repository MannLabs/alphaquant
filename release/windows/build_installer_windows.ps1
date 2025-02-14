# Build the installer for Windows.
# This script must be run from the root of the repository.

Remove-Item -Recurse -Force -ErrorAction SilentlyContinue ./build_pyinstaller
Remove-Item -Recurse -Force -ErrorAction SilentlyContinue ./dist_pyinstaller


# substitute X.Y.Z-devN with X.Y.Z.devN
$WHL_NAME = (Get-ChildItem -Path "dist" -Filter "*.whl").Name

## install alphamap
#git clone alphamap@2345rf
#cd alphamap
#pip install .
#remove alphamap from alphaquant requirements.txt
#
## re-build wheel
#rm -rf dist ./*.egg-info
#python -m build


pip install "dist/$WHL_NAME[stable,gui-stable]"

# Creating the stand-alone pyinstaller folder
pip install pyinstaller
pyinstaller release/pyinstaller/alphaquant.spec --distpath dist_pyinstaller --workpath build_pyinstaller -y
