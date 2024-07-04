import json
import pandas as pd
import requests
import os
import argparse
import time

# Function to get PubChem ID using SMILES or compound name
def get_pubchem_id(smiles=None, name=None):
    if smiles:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/cids/JSON"
    elif name:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/cids/JSON"
    else:
        return None

    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        if 'IdentifierList' in data:
            return data['IdentifierList']['CID'][0]
    return None

# Function to process a single file and extract data
def process_file(file, existing_mibig_accessions):
    data = []
    with open(file) as f:
        content = json.load(f)
        mibig_accession = content['cluster'].get('mibig_accession')
        if mibig_accession in existing_mibig_accessions:
            print(f"Skipping already processed file with mibig_accession: {mibig_accession}")
            return data

        ncbi_tax_id = content['cluster'].get('ncbi_tax_id')
        organism_name = content['cluster'].get('organism_name')
        compounds = content['cluster'].get('compounds', [])
        for compound_index, compound in enumerate(compounds, start=1):
            compound_name = compound.get('compound', '')
            smiles = compound.get('chem_struct', '')
            pubchem_id = get_pubchem_id(smiles=smiles, name=compound_name)
            data.append({
                "mibig_accession": mibig_accession,
                "ncbi_tax_id": ncbi_tax_id,
                "organism_name": organism_name,
                "compound": compound_name,
                "smiles": smiles,
                "pubchem_id": pubchem_id
            })
            print(f"Processed compound {compound_index}/{len(compounds)} in file with mibig_accession: {mibig_accession}")
            # Ensure we do not exceed 5 requests per second
            time.sleep(0.2)
    return data

# Function to process files and extract data
def process_files(directory, output_file):
    # Initialize a list to store data
    data = []

    # List all JSON files in the directory
    files = [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith('.json')]
    total_files = len(files)

    # Check for existing summary file and read existing mibig_accessions
    if os.path.exists(output_file):
        existing_df = pd.read_csv(output_file)
        existing_mibig_accessions = set(existing_df['mibig_accession'])
    else:
        existing_df = pd.DataFrame()
        existing_mibig_accessions = set()

    # Extract relevant information from each file
    for file_index, file in enumerate(files, start=1):
        file_data = process_file(file, existing_mibig_accessions)
        if file_data:
            file_df = pd.DataFrame(file_data)
            if existing_df.empty:
                existing_df = file_df
            else:
                existing_df = pd.concat([existing_df, file_df], ignore_index=True)
            existing_df.to_csv(output_file, index=False)
            print(f"Summary written to {output_file} after processing file {file_index}/{total_files}")

    return existing_df

# Main function
def main():
    parser = argparse.ArgumentParser(description="Process MIBiG JSON files to extract compound information and PubChem IDs.")
    parser.add_argument('directory', type=str, help='Path to the directory containing JSON files')
    parser.add_argument('--output', type=str, default='mibig_summary.csv', help='Output CSV file name')

    args = parser.parse_args()
    directory = args.directory
    output_file = args.output

    print(f"Processing files in directory: {directory}")
    df = process_files(directory, output_file)

    print(f"Final summary table saved to {output_file}")

if __name__ == "__main__":
    main()
