#!bash

# Initial cleanup
rm -rf dist
rm -rf build
cd ../..
rm -rf dist
rm -rf build

# Creating a conda environment
conda create -n alphaquant_installer python=3.8 -y
conda activate alphaquant_installer

# Creating the wheel
python setup.py sdist bdist_wheel

# Setting up the local package
cd release/one_click_linux_gui
# Make sure you include the required extra packages and always use the stable or very-stable options!
pip install "../../dist/alphaquant-0.0.7-py3-none-any.whl[stable]"

# Creating the stand-alone pyinstaller folder
pip install pyinstaller==4.10
pyinstaller ../pyinstaller/alphaquant.spec -y
conda deactivate

# If needed, include additional source such as e.g.:
# cp ../../alphaquant/data/*.fasta dist/alphaquant/data
# WARNING: this probably does not work!!!!

# Wrapping the pyinstaller folder in a .deb package
mkdir -p dist/alphaquant_gui_installer_linux/usr/local/bin
mv dist/alphaquant dist/alphaquant_gui_installer_linux/usr/local/bin/alphaquant
mkdir dist/alphaquant_gui_installer_linux/DEBIAN
cp control dist/alphaquant_gui_installer_linux/DEBIAN
dpkg-deb --build --root-owner-group dist/alphaquant_gui_installer_linux/
