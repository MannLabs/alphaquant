#!/bin/bash
set -e -u

INSTALLER_DATA_DIR=${INSTALLER_DATA_DIR:-build_pyinstaller_data}
ALPHAMAP_DATA_DIR="${INSTALLER_DATA_DIR}/alphamap/data"
ALPHAQUANT_RESOURCES_DIR="${INSTALLER_DATA_DIR}/alphaquant/resources"

curl_github() {
    local args=(
        -L
        -f
        --retry 3
        --retry-delay 2
        -H "Accept: application/vnd.github+json"
        -H "User-Agent: alphaquant-release-build"
    )

    local github_token="${GITHUB_TOKEN:-${GH_TOKEN:-}}"
    if [ -n "$github_token" ]; then
        args+=(-H "Authorization: Bearer ${github_token}")
    fi

    curl "${args[@]}" "$@"
}

download_datashare_zip() {
    local url=$1
    local label=$2
    local temp_zip

    temp_zip=$(mktemp)
    echo "Downloading ${label} from datashare..."
    if ! curl -L -f --retry 3 --retry-delay 2 "${url}" -o "${temp_zip}"; then
        echo "Error: Failed to download ${label} from datashare"
        rm -f "${temp_zip}"
        exit 1
    fi

    echo "Extracting ${label}..."
    unzip -o "${temp_zip}" -d "${ALPHAQUANT_RESOURCES_DIR}"
    rm -f "${temp_zip}"
}

rm -rf "${INSTALLER_DATA_DIR}"
mkdir -p "${ALPHAMAP_DATA_DIR}" "${ALPHAQUANT_RESOURCES_DIR}"

echo "Downloading AlphaMap FASTA and CSV files..."
if ! DOWNLOAD_LIST=$(curl_github https://api.github.com/repos/MannLabs/alphamap/contents/alphamap/data?ref=main); then
    echo "Error: Failed to fetch AlphaMap file list from GitHub API"
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

    filename=$(basename "$url")
    echo "Downloading ${filename}..."
    if ! curl_github "$url" -o "${ALPHAMAP_DATA_DIR}/${filename}"; then
        echo "Error: Failed to download ${filename}"
        exit 1
    fi
  done

alphamap_file_count=$(find "${ALPHAMAP_DATA_DIR}" -type f \( -name "*.fasta" -o -name "*.csv" \) | wc -l)
echo "Downloaded ${alphamap_file_count} AlphaMap files"
if [ "${alphamap_file_count}" -eq 0 ]; then
    echo "Error: No AlphaMap files were downloaded"
    exit 1
fi

download_datashare_zip "https://datashare.biochem.mpg.de/s/ezPzeqStEgDD8gg/download" "reference databases"
download_datashare_zip "https://datashare.biochem.mpg.de/s/stH9pmNe6O9CRHG/download" "phosphopred databases"

REFERENCE_DB_DIR="${ALPHAQUANT_RESOURCES_DIR}/reference_databases"
PHOSPHOPRED_DB_DIR="${ALPHAQUANT_RESOURCES_DIR}/phosphopred_databases"
HUMAN_PHOSPHO_FILE="${PHOSPHOPRED_DB_DIR}/human_uniprot_reviewed_phos_prob.tsv"

if [ ! -d "${REFERENCE_DB_DIR}" ]; then
    echo "Error: reference_databases directory not found at ${REFERENCE_DB_DIR}"
    exit 1
fi

if [ ! -d "${PHOSPHOPRED_DB_DIR}" ]; then
    echo "Error: phosphopred_databases directory not found at ${PHOSPHOPRED_DB_DIR}"
    exit 1
fi

if [ ! -f "${HUMAN_PHOSPHO_FILE}" ]; then
    echo "Error: Expected phosphopred database file not found at ${HUMAN_PHOSPHO_FILE}"
    exit 1
fi

echo "Installer data staged in ${INSTALLER_DATA_DIR}"
