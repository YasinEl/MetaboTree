import pandas as pd
import argparse
from Bio import Entrez

def load_data(redu_path, sparql_path, ncbi_path):
    # Load the CSV files
    df_redu = pd.read_csv(redu_path)
    dt_sparql = pd.read_csv(sparql_path)
    df_ncbi = pd.read_csv(ncbi_path)
    
    # Split the 'NCBITaxonomy' column into 'NCBI' and 'tax_name'
    df_redu[['NCBI', 'tax_name']] = df_redu['NCBITaxonomy'].str.split('|', expand=True)
    
    # Convert NCBI column to integer
    df_redu['NCBI'] = df_redu['NCBI'].astype(int)

    # Collect existing NCBI IDs from df_redu
    existing_ncbi_ids = set(df_redu['NCBI'])

    # Prepare df_ncbi and dt_sparql by ensuring NCBI ID columns are integer
    df_ncbi['NCBI'] = df_ncbi['NCBI'].astype(int)
    dt_sparql['NCBI'] = dt_sparql['ncbiTaxonomyID'].astype(int)
    dt_sparql = dt_sparql[['NCBI']]
    dt_sparql = dt_sparql.drop_duplicates()

    # Filter new entries from df_ncbi that are not in existing_ncbi_ids
    new_ncbi_df = df_ncbi[~df_ncbi['ID'].isin(existing_ncbi_ids)]
    new_ncbi_df = new_ncbi_df.rename(columns={'ID': 'NCBI'})
    new_ncbi_df['tax_name'] = 'unknown'
    new_ncbi_df['DataSource'] = 'ncbi'

    # Filter new entries from dt_sparql that are not in existing_ncbi_ids
    new_sparql_df = dt_sparql[~dt_sparql['ncbiTaxonomyID'].isin(existing_ncbi_ids)]
    new_sparql_df = new_sparql_df.rename(columns={'ncbiTaxonomyID': 'NCBI', 'organismLabel': 'tax_name'})
    new_sparql_df['DataSource'] = 'sparql'

    # Combine all data avoiding any duplicates
    combined_new_entries = pd.concat([new_ncbi_df, new_sparql_df])
    combined_new_entries = combined_new_entries.drop_duplicates(subset=['NCBI'])

    df_redu = pd.concat([df_redu, combined_new_entries], ignore_index=True)

   
    return df_redu

def main(args):
    df_redu = load_data(args.redu_path, args.sparql_path, args.ncbi_path)
    df_redu.to_csv('all_organism_table.csv', index=False)
    print(f"Data saved to all_organism_table.csv")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Extend df_redu with additional NCBI IDs from CSV files.")
    parser.add_argument('redu_path', help='Path to the REDU CSV file')
    parser.add_argument('sparql_path', help='Path to the SPARQL CSV file')
    parser.add_argument('ncbi_path', help='Path to the NCBI CSV file')
    
    args = parser.parse_args()
    main(args)
