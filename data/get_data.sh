#!/bin/bash

# Define the directory where the script is located
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

download_file() {
    local url=$1
    local path=$2
    if [ ! -f "$path" ]; then
        echo "Downloading $(basename "$path")..."
        curl -o "$path" "$url" --verbose
        if [ $? -ne 0 ]; then
            echo "Failed to download $(basename "$path")."
            return 1
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

# Download the files to the script's directory if they don't exist
download_file "$URL1" "${DIR}/biosystems_pcsubstance.gz"
download_file "$URL2" "${DIR}/biosystems_taxonomy.gz"
download_file "$URL3" "${DIR}/ott3.7.tgz"

# Unzip the files in the same directory if they don't already exist
if [ ! -f "$DIR/biosystems_pcsubstance" ]; then
    echo "Unzipping biosystems_pcsubstance.gz..."
    gunzip -c "$DIR/biosystems_pcsubstance.gz" > "$DIR/biosystems_pcsubstance"
fi

if [ ! -f "$DIR/biosystems_taxonomy" ]; then
    echo "Unzipping biosystems_taxonomy.gz..."
    gunzip -c "$DIR/biosystems_taxonomy.gz" > "$DIR/biosystems_taxonomy"
fi

if [ ! -f "$DIR/taxonomy.tsv" ]; then
    echo "Listing contents of ott3.7.tgz..."
    tar -tzf "${DIR}/ott3.7.tgz"

    echo "Extracting taxonomy.tsv from ott3.7.tgz..."
    tar --extract --file="${DIR}/ott3.7.tgz" --directory="${DIR}" --strip-components=1 "ott3.7/taxonomy.tsv"
    if [ -f "${DIR}/taxonomy.tsv" ]; then
        echo "taxonomy.tsv has been extracted."
    else
        echo "Failed to extract taxonomy.tsv."
    fi
fi

# Remove the original .gz and .zip files if they have been successfully processed
if [ -f "$DIR/biosystems_pcsubstance" ]; then
    rm "$DIR/biosystems_pcsubstance.gz"
fi

if [ -f "$DIR/biosystems_taxonomy" ]; then
    rm "$DIR/biosystems_taxonomy.gz"
fi

if [ -f "$DIR/taxonomy.tsv" ]; then
    rm "${DIR}/ott3.7.tgz"
fi

echo "Compressed files have been removed."

echo "Files downloaded and uncompressed successfully."
