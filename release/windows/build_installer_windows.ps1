# Build the installer for Windows.
# This script must be run from the root of the repository.

Remove-Item -Recurse -Force -ErrorAction SilentlyContinue ./build
Remove-Item -Recurse -Force -ErrorAction SilentlyContinue ./dist
Remove-Item -Recurse -Force -ErrorAction SilentlyContinue ./dist_pyinstaller


# substitute X.Y.Z-devN with X.Y.Z.devN
$WHL_NAME = (Get-ChildItem -Path "dist" -Filter "*.whl").Name
pip install "dist/$WHL_NAME[stable,gui-stable]"

# Creating the stand-alone pyinstaller folder
pip install pyinstaller
pyinstaller release/pyinstaller/alphaquant.spec --distpath dist_pyinstaller --workpath build_pyinstaller -y
