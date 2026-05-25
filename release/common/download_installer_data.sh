#!/bin/bash
set -e -u

INSTALLER_DATA_DIR=${INSTALLER_DATA_DIR:-build_pyinstaller_data}
ALPHAMAP_DATA_DIR="${INSTALLER_DATA_DIR}/alphamap/data"
ALPHAQUANT_RESOURCES_DIR="${INSTALLER_DATA_DIR}/alphaquant/resources"

curl_with_retries() {
    local args=(
        -L
        -f
        --retry 3
        --retry-delay 2
        -H "User-Agent: alphaquant-release-build"
    )

    curl "${args[@]}" "$@"
}

list_alphamap_data_downloads() {
    python -c 'from urllib.parse import quote
from alphamap.organisms_data import all_organisms
base_url = "https://raw.githubusercontent.com/MannLabs/alphamap/main/alphamap/data/"
seen = set()
for organism in all_organisms.values():
    for key in ("fasta_name", "uniprot_name"):
        name = organism[key]
        if name not in seen:
            seen.add(name)
            print(f"{name}\t{base_url}{quote(name)}")'
}

download_datashare_zip() {
    local url=$1
    local label=$2
    local temp_zip

    temp_zip=$(mktemp)
    echo "Downloading ${label} from datashare..."
    if ! curl_with_retries "${url}" -o "${temp_zip}"; then
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
list_alphamap_data_downloads | \
  while IFS=$'\t' read -r filename url; do
    if [ -z "$filename" ] || [ -z "$url" ]; then
        echo "Warning: Empty URL detected, skipping..."
        continue
    fi

    echo "Downloading ${filename}..."
    if ! curl_with_retries "$url" -o "${ALPHAMAP_DATA_DIR}/${filename}"; then
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
