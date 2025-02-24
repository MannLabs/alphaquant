# Build the install package for Windows.
# This script must be run from the root of the repository after running build_installer_windows.ps1

# A directory for handling the AlphaMap (no typo!) data, analogous to the Mac installer, is added via the .inno file

# Wrapping the pyinstaller folder in a .exe package
&  "C:\Program Files (x86)\Inno Setup 6\ISCC.exe" .\release\windows\alphaquant_innoinstaller.iss
