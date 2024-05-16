#!/bin/bash

# Check if the TSV file path is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 path_to_tsv_file"
    exit 1
fi

# Path to your TSV file
tsv_file=$1

# Extract the molecule name from the TSV file using awk
molecule_name=$(awk -F'\t' 'NR==1 {print $2}' "$tsv_file")

# Replace spaces with underscores in the molecule name to make it filename-safe
molecule_name_safe=$(echo "$molecule_name" | tr ' ' '_' | tr -d '[:punct:]')

# Rename the file
mv wikidata/sparql_query.csv "sparql_query_${molecule_name_safe}.csv"
