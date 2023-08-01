#!bash

# Initial cleanup
rm -rf dist
rm -rf build
FILE=alphaquant.pkg
if test -f "$FILE"; then
  rm alphaquant.pkg
fi
cd ../..
rm -rf dist
rm -rf build

# Creating a conda environment
conda create -n alphaquantinstaller python=3.8 -y
conda activate alphaquantinstaller

# Creating the wheel
python setup.py sdist bdist_wheel

# Setting up the local package
cd release/one_click_macos_gui
pip install "../../dist/alphaquant-0.0.6-py3-none-any.whl[stable]"

# Creating the stand-alone pyinstaller folder
pip install pyinstaller==4.10
pyinstaller ../pyinstaller/alphaquant.spec -y
conda deactivate

# If needed, include additional source such as e.g.:
# cp ../../alphaquant/data/*.fasta dist/alphaquant/data

# Wrapping the pyinstaller folder in a .pkg package
mkdir -p dist/alphaquant/Contents/Resources
cp ../logos/alpha_logo.icns dist/alphaquant/Contents/Resources
mv dist/alphaquant_gui dist/alphaquant/Contents/MacOS
cp Info.plist dist/alphaquant/Contents
cp alphaquant_terminal dist/alphaquant/Contents/MacOS
cp ../../LICENSE Resources/LICENSE
cp ../logos/alpha_logo.png Resources/alpha_logo.png
chmod 777 scripts/*

pkgbuild --root dist/alphaquant --identifier de.mpg.biochem.alphaquant.app --version 0.0.6 --install-location /Applications/alphaquant.app --scripts scripts alphaquant.pkg
productbuild --distribution distribution.xml --resources Resources --package-path alphaquant.pkg dist/alphaquant_gui_installer_macos.pkg
