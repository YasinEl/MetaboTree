import pandas as pd
import argparse
from Bio import Entrez

def load_data(redu_path, sparql_path, ncbi_path):
    # Load the CSV files
    df_redu = pd.read_csv(redu_path)
    dt_sparql = pd.read_csv(sparql_path)
    df_ncbi = pd.read_csv(ncbi_path)
    
    # Convert NCBI column to integer
    df_redu['NCBI'] = df_redu['NCBI'].astype(int)

    # Collect existing NCBI IDs from df_redu
    existing_ncbi_ids = set(df_redu['NCBI'])

    # Prepare df_ncbi and dt_sparql by ensuring NCBI ID columns are integer
    df_ncbi['NCBI'] = df_ncbi['NCBI'].astype(int)
    df_ncbi = df_ncbi.drop_duplicates()
    df_ncbi['Database'] = 'NCBI_biosystems'
    dt_sparql['NCBI'] = dt_sparql['ncbiTaxonomyID'].astype(int)
    dt_sparql = dt_sparql[['NCBI']]
    dt_sparql = dt_sparql.drop_duplicates()
    dt_sparql['Database'] = 'Wikidata'


    df_databases = pd.concat([dt_sparql, df_ncbi])


    # Merge df_redu with df_databases to get the Database column
    df_merged = pd.merge(df_redu, df_databases, on='NCBI', how='left', suffixes=('', '_from_db'))

    # If there are duplicates in df_merged for the same NCBI, combine Database values
    df_merged['Database'] = df_merged.groupby('NCBI')['Database'].transform(lambda x: ','.join(x.dropna().unique()))

    # Drop duplicates and reset index
    df_merged = df_merged.drop_duplicates().reset_index(drop=True)

    # Find NCBIs present in df_databases but not in df_redu
    missing_ncbi = df_databases[~df_databases['NCBI'].isin(df_redu['NCBI'])]

    # Add missing NCBIs to df_merged
    df_final = pd.concat([df_merged, missing_ncbi], ignore_index=True, sort=False)

    # Fill NaN values in Database column with empty string
    df_final['Database'] = df_final['Database'].fillna('')
   
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
