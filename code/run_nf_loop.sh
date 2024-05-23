#!/bin/bash

# Path to your TSV file
INPUT_TSV="data/plant_drugs_2.tsv"

# Read the TSV file and loop through each row
while IFS=$'\t' read -r usi pubchemid
do
    # Skip the header row
    if [ "$usi" != "usi" ]; then
        echo "Running Nextflow for USI: $usi, PubChemID: $pubchemid"
        nextflow run nf_workflow.nf --usi "$usi" --pubchemid "$pubchemid" -resume
    fi
done < "$INPUT_TSV"
