#!/bin/bash

# Define the directory where the script is located
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

download_file() {
    local url=$1
    local path=$2
    if [ ! -f "$path" ]; then
        echo "Downloading $(basename "$path")..."
        curl -o "$path" "$url"
        if [ $? -ne 0 ]; then
            echo "Failed to download $(basename "$path")."
            return 1
        fi

        if [[ "$url" == ftp* ]]; then
            # Check for FTP response code 450 (file unavailable or insufficient storage)
            curl -I "$url" 2>&1 | grep "450" && { echo "Server returned error 450 for $(basename "$path")."; return 1; }
        fi

        echo "$(basename "$path") has been downloaded to ${DIR}"
    else
        echo "$(basename "$path") already exists. Skipping download."
    fi
}

# URLs of the files to download
URL1="ftp://ftp.ncbi.nih.gov/pub/biosystems/biosystems.20170421/biosystems_pcsubstance.gz"
URL2="ftp://ftp.ncbi.nih.gov/pub/biosystems/biosystems.20170421/biosystems_taxonomy.gz"
URL3="https://files.opentreeoflife.org/ott/ott3.7/ott3.7.tgz"
URL4="https://files.opentreeoflife.org/synthesis/opentree14.9/opentree14.9tree.tgz"
URL5="https://redu.gnps2.org/dump"

# Download the files to the script's directory if they don't exist
download_file "$URL1" "${DIR}/biosystems_pcsubstance.gz"
download_file "$URL2" "${DIR}/biosystems_taxonomy.gz"
download_file "$URL3" "${DIR}/ott3.7.tgz"
download_file "$URL4" "${DIR}/opentree14.9tree.tgz"

# Download the redu.tsv file
if [ ! -f "${DIR}/redu.tsv" ]; then
    echo "Downloading redu.tsv..."
    wget -q -O "${DIR}/redu.tsv" "$URL5"
    if [ $? -eq 0 ]; then
        echo "redu.tsv has been downloaded to ${DIR}"
    else
        echo "Failed to download redu.tsv."
    fi
else
    echo "redu.tsv already exists. Skipping download."
fi

# Unzip the files in the same directory if they don't already exist
if [ -f "${DIR}/biosystems_pcsubstance.gz" ]; then
    echo "Unzipping biosystems_pcsubstance.gz..."
    gunzip -c "${DIR}/biosystems_pcsubstance.gz" > "${DIR}/biosystems_pcsubstance" || { echo "Failed to unzip biosystems_pcsubstance.gz"; exit 1; }
fi

if [ -f "${DIR}/biosystems_taxonomy.gz" ]; then
    echo "Unzipping biosystems_taxonomy.gz..."
    gunzip -c "${DIR}/biosystems_taxonomy.gz" > "${DIR}/biosystems_taxonomy" || { echo "Failed to unzip biosystems_taxonomy.gz"; exit 1; }
fi

if [ -f "${DIR}/ott3.7.tgz" ] && [ ! -f "${DIR}/taxonomy.tsv" ]; then
    echo "Listing contents of ott3.7.tgz..."
    tar -tzf "${DIR}/ott3.7.tgz"

    echo "Extracting taxonomy.tsv from ott3.7.tgz..."
    tar --extract --file="${DIR}/ott3.7.tgz" --directory="${DIR}" --strip-components=1 "ott3.7/taxonomy.tsv"
    if [ $? -eq 0 ]; then
        echo "taxonomy.tsv has been extracted."
    else
        echo "Failed to extract taxonomy.tsv."
    fi
fi

if [ -f "${DIR}/opentree14.9tree.tgz" ] && [ ! -f "${DIR}/labelled_supertree.tre" ]; then
    echo "Listing contents of opentree14.9tree.tgz..."
    tar -tzf "${DIR}/opentree14.9tree.tgz"

    echo "Extracting labelled_supertree.tre from opentree14.9tree.tgz..."
    tar --extract --file="${DIR}/opentree14.9tree.tgz" --directory="${DIR}" --strip-components=1 "opentree14.9_tree/labelled_supertree/labelled_supertree.tre"
    if [ $? -eq 0 ]; then
        echo "labelled_supertree.tre has been extracted."
        mv "${DIR}/labelled_supertree/labelled_supertree.tre" "${DIR}/" || { echo "Failed to move labelled_supertree.tre to ${DIR}"; exit 1; }
    else
        echo "Failed to extract labelled_supertree.tre."
    fi
fi

# Remove the original .gz and .tgz files if they have been successfully processed
if [ -f "${DIR}/biosystems_pcsubstance" ]; then
    rm "${DIR}/biosystems_pcsubstance.gz"
fi

if [ -f "${DIR}/biosystems_taxonomy" ]; then
    rm "${DIR}/biosystems_taxonomy.gz"
fi

if [ -f "${DIR}/taxonomy.tsv" ]; then
    rm "${DIR}/ott3.7.tgz"
fi

if [ -f "${DIR}/labelled_supertree.tre" ]; then
    rm "${DIR}/opentree14.9tree.tgz"
fi

echo "Compressed files have been removed."

echo "Files downloaded and uncompressed successfully."
