# Build the installer for Windows.
# This script must be run from the root of the repository.

Remove-Item -Recurse -Force -ErrorAction SilentlyContinue ./build_pyinstaller
Remove-Item -Recurse -Force -ErrorAction SilentlyContinue ./dist_pyinstaller


# substitute X.Y.Z-devN with X.Y.Z.devN
$WHL_NAME = (Get-ChildItem -Path "dist" -Filter "*.whl").Name

## install alphamap
git clone alphamap@c65381a11a0f25fe04822846863520afc966b967
cd alphamap
pip install .
# Remove alphamap from requirements files
(Get-Content "requirements/requirements.txt") | Where-Object { $_ -notmatch "^alphamap==" } | Set-Content "requirements/requirements.txt"
(Get-Content "requirements/requirements_loose.txt") | Where-Object { $_ -notmatch "^alphamap$" } | Set-Content "requirements/requirements_loose.txt"

# re-build wheel
cd ..
rm -rf dist ./*.egg-info
python -m build


pip install "dist/$WHL_NAME[stable,gui-stable]"

# Creating the stand-alone pyinstaller folder
pip install pyinstaller
pyinstaller release/pyinstaller/alphaquant.spec --distpath dist_pyinstaller --workpath build_pyinstaller -y
