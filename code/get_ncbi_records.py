import pandas as pd
import argparse


def main(ncbi_tax_path, ncbi_comp_path, cid_file):
    # Load data
    df_tax = pd.read_csv(ncbi_tax_path, sep='\t', names=['bsid', 'TaxID', 'score'], dtype={'bsid': int, 'TaxID': str, 'score': float})
    df_comp = pd.read_csv(ncbi_comp_path, sep='\t', names=['bsid', 'cid', 'score'], dtype={'bsid': 'Int64', 'cid': 'Int64', 'score': float})

    # Load PubChem IDs from the file
    with open(cid_file, 'r') as f:
        pubchem_ids = [int(line.strip()) for line in f]

    # Filter df_comp to include only rows where cid is in pubchem_ids
    filtered_comp = df_comp[df_comp['cid'].isin(pubchem_ids)]

    # Get unique bsids from filtered_comp
    unique_bsids = filtered_comp['bsid'].unique()

    # Filter df_tax to include only rows with bsids that are in unique_bsids
    filtered_tax = df_tax[df_tax['bsid'].isin(unique_bsids)]

    # Get unique TaxIDs from filtered_tax
    result_df = pd.DataFrame({'NCBI': filtered_tax['TaxID'].unique()})


    # Print and save results
    print(result_df)
    result_df.to_csv('ncbi_records.csv', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a CSV file of NCBI records based on given PubChem IDs.")
    parser.add_argument("ncbi_tax_path", help="Path to the TSV file containing NCBI taxonomy data.")
    parser.add_argument("ncbi_comp_path", help="Path to the TSV file containing NCBI compound data.")
    parser.add_argument("cid_file", help="Path to the TSV file containing a list of PubChem IDs.")
    args = parser.parse_args()

    main(args.ncbi_tax_path, args.ncbi_comp_path, args.cid_file)
