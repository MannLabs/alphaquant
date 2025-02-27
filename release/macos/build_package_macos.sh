#!/bin/bash
set -u -e

# Build the install package for MacOS.
# This script must be run from the root of the repository after running build_installer_macos.sh

PACKAGE_NAME=alphaquant
# BUILD_NAME is taken from environment variables, e.g. alphaquant-1.2.3-macos-darwin-arm64 or alphaquant-1.2.3-macos-darwin-x64
rm -rf ${BUILD_NAME}.pkg

# If needed, include additional source such as e.g.:
# cp ../../alphaquant/data/*.fasta dist/alphaquant/data

# Wrapping the pyinstaller folder in a .pkg package
CONTENTS_FOLDER=dist_pyinstaller/${PACKAGE_NAME}/Contents

mkdir -p ${CONTENTS_FOLDER}/Resources
cp release/logos/alpha_logo.icns ${CONTENTS_FOLDER}/Resources
mv dist_pyinstaller/alphaquant_gui ${CONTENTS_FOLDER}/MacOS
cp release/macos/Info.plist ${CONTENTS_FOLDER}
cp release/macos/alphaquant_terminal ${CONTENTS_FOLDER}/MacOS
cp ./LICENSE ${CONTENTS_FOLDER}/Resources/LICENSE
cp release/logos/alpha_logo.png ${CONTENTS_FOLDER}/Resources



# link _internal folder containing the python libraries to the Frameworks folder where they are expected
# to avoid e.g. "Failed to load Python shared library '/Applications/AlphaMap.app/Contents/Frameworks/libpython3.8.dylib'"
cd ${CONTENTS_FOLDER}
ln -s ./MacOS/_internal ./Frameworks
cd -

#make directory for AlphaMap. This is where AlphaMap stores downloaded data, such as fasta files
mkdir -p ${CONTENTS_FOLDER}/Frameworks/alphamap/data/

####
####Download all AlphaMap FASTA and CSV files from GitHub, which are needed for the further analyses. There is a lot of error checking to ensure that the files get actually added during the build
echo "Starting downloads of FASTA and CSV files..."
DOWNLOAD_LIST=$(curl -L -f https://api.github.com/repos/MannLabs/alphamap/contents/alphamap/data?ref=main)
if [ $? -ne 0 ]; then
    echo "Error: Failed to fetch file list from GitHub API"
    exit 1
fi

echo "$DOWNLOAD_LIST" | \
  grep "\"download_url\".*\.\(fasta\|csv\)\"" | \
  cut -d '"' -f 4 | \
  while read url; do
    if [ -z "$url" ]; then
        echo "Warning: Empty URL detected, skipping..."
        continue
    fi
    filename=$(basename $url)
    echo "Downloading $filename..."
    if ! curl -L -f "$url" -o "${CONTENTS_FOLDER}/Frameworks/alphamap/data/$filename"; then
        echo "Error: Failed to download $filename"
        exit 1
    fi
    echo "Successfully downloaded $filename"
  done

# Verify downloads
file_count=$(ls -1 ${CONTENTS_FOLDER}/Frameworks/alphamap/data/*.{fasta,csv} 2>/dev/null | wc -l)
echo "Downloaded $file_count files"
if [ $file_count -eq 0 ]; then
    echo "Error: No files were downloaded"
    exit 1
fi
####
###Download section complete


chmod 777 release/macos/scripts/*

pkgbuild --root dist_pyinstaller/${PACKAGE_NAME} --identifier de.mpg.biochem.${PACKAGE_NAME}.app --version 0.1.5 --install-location /Applications/${PACKAGE_NAME}.app --scripts release/macos/scripts ${PACKAGE_NAME}.pkg
productbuild --distribution release/macos/distribution.xml --resources release/macos/Resources --package-path ${PACKAGE_NAME}.pkg ${BUILD_NAME}.pkg
