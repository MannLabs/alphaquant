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
cd release/one_click_windows_gui
# Make sure you include the required extra packages and always use the stable or very-stable options!
pip install "../../dist/alphaquant-0.0.8-py3-none-any.whl[stable]"

# Creating the stand-alone pyinstaller folder
pip install pyinstaller==4.10
pyinstaller ../pyinstaller/alphaquant.spec -y
conda deactivate

# If needed, include additional source such as e.g.:
# cp ../../alphaquant/data/*.fasta dist/alphaquant/data

# Wrapping the pyinstaller folder in a .exe package
"C:\Program Files (x86)\Inno Setup 6\ISCC.exe" alphaquant_innoinstaller.iss
# WARNING: this assumes a static location for innosetup
