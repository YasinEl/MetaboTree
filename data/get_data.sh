#!/bin/bash

# Define the directory where the script is located
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

download_file() {
    local url=$1
    local path=$2
    echo "Downloading $(basename "$path")..."
    curl -o "$path" "$url"
    echo "$(basename "$path") has been downloaded to ${DIR}"
}

# URLs of the files to download
URL1="ftp://ftp.ncbi.nih.gov/pub/biosystems/biosystems.20170421/biosystems_pcsubstance.gz"
URL2="ftp://ftp.ncbi.nih.gov/pub/biosystems/biosystems.20170421/biosystems_taxonomy.gz"
URL3="https://files.opentreeoflife.org/ott/ott3.7/ott3.7.tgz"

# Download the files to the script's directory
download_file "$URL1" "${DIR}/biosystems_pcsubstance.gz"
download_file "$URL2" "${DIR}/biosystems_taxonomy.gz"
download_file "$URL3" "${DIR}/ott3.7.tgz"
download_file "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip" "${DIR}/taxdmp.zip"

# Unzip the files in the same directory
gunzip -c "$DIR/biosystems_pcsubstance.gz" > "$DIR/biosystems_pcsubstance"
gunzip -c "$DIR/biosystems_taxonomy.gz" > "$DIR/biosystems_taxonomy"

echo "Unzipping taxdmp.zip..."
unzip -j "${DIR}/taxdmp.zip" "names.dmp" -d "${DIR}"
echo "names.dmp has been extracted."

# List contents of ott3.7.tgz to verify the structure
echo "Listing contents of ott3.7.tgz..."
tar -tzf "${DIR}/ott3.7.tgz"

# Extract taxonomy.tsv
echo "Extracting taxonomy.tsv from ott3.7.tgz..."
tar --extract --file="${DIR}/ott3.7.tgz" --directory="${DIR}" --strip-components=1 "ott3.7/taxonomy.tsv"
if [ -f "${DIR}/taxonomy.tsv" ]; then
    echo "taxonomy.tsv has been extracted."
else
    echo "Failed to extract taxonomy.tsv."
fi

# Remove the original .gz and .zip files
rm "$DIR/biosystems_pcsubstance.gz"
rm "$DIR/biosystems_taxonomy.gz"
rm "${DIR}/taxdmp.zip"
rm "${DIR}/ott3.7.tgz"
echo "Compressed files have been removed."

echo "Files downloaded and uncompressed successfully."
