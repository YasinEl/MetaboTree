import pandas as pd
import argparse
from Bio import Entrez

def load_data(redu_path, sparql_path, ncbi_path):
    df_redu = pd.read_csv(redu_path, sep='\t', low_memory=False)
    dt_sparql = pd.read_csv(sparql_path)
    df_ncbi = pd.read_csv(ncbi_path)

    # Convert NCBI column to integer and strip any whitespace
    df_redu['NCBI'] = df_redu['NCBI'].astype(str).str.strip().astype(int)

    # Print debug info
    print(f"df_redu NCBI unique count: {df_redu['NCBI'].nunique()}")
    print(f"df_redu NCBI unique values (sample): {df_redu['NCBI'].unique()[:10]}")

    # Collect existing NCBI IDs from df_redu
    existing_ncbi_ids = set(df_redu['NCBI'])
    print(f"Existing NCBI IDs in df_redu (sample): {list(existing_ncbi_ids)[:10]}")

    # Prepare df_ncbi and dt_sparql by ensuring NCBI ID columns are integer
    if not df_ncbi.empty:
        df_ncbi['NCBI'] = df_ncbi['NCBI'].astype(str).str.strip()
        df_ncbi['NCBI'] = df_ncbi['NCBI'].astype(int)
        df_ncbi = df_ncbi.drop_duplicates()
        df_ncbi['Database'] = 'NCBI_biosystems'

        # Print debug info
        print(f"df_ncbi NCBI unique count: {df_ncbi['NCBI'].nunique()}")
        print(f"df_ncbi NCBI unique values (sample): {df_ncbi['NCBI'].unique()[:10]}")
    else:
        print("df_ncbi is empty")

    if not dt_sparql.empty:
        dt_sparql['NCBI'] = dt_sparql['ncbiTaxonomyID'].astype(str).str.strip().astype(int)
        dt_sparql = dt_sparql[['NCBI']]
        dt_sparql = dt_sparql.drop_duplicates()
        dt_sparql['Database'] = 'Wikidata'

        # Print debug info
        print(f"dt_sparql NCBI unique count: {dt_sparql['NCBI'].nunique()}")
        print(f"dt_sparql NCBI unique values (sample): {dt_sparql['NCBI'].unique()[:10]}")
    else:
        print("dt_sparql is empty")

    # Concatenate dataframes
    df_databases = pd.concat([dt_sparql, df_ncbi])

    # Print debug info
    print(f"df_databases NCBI unique count: {df_databases['NCBI'].nunique()}")
    print(f"df_databases NCBI unique values (sample): {df_databases['NCBI'].unique()[:10]}")
    print(df_databases.head())

    # Merge df_redu with df_databases to get the Database column
    df_merged = pd.merge(df_redu, df_databases, on='NCBI', how='left', suffixes=('', '_from_db'))

    # Print debug info
    print(f"df_merged row count: {len(df_merged)}")
    print(df_merged.head())

    # Check for any unexpected issues in the NCBI columns
    print("df_redu NCBI values (sample):")
    print(df_redu['NCBI'].unique()[:10])

    print("df_databases NCBI values (sample):")
    print(df_databases['NCBI'].unique()[:10])

    print("df_merged NCBI values (sample):")
    print(df_merged['NCBI'].unique()[:10])

    # Check for rows where Database is not null after the merge
    non_null_merged = df_merged[df_merged['Database'].notna()]
    print(f"Non-null merged row count: {len(non_null_merged)}")
    print(non_null_merged.head())

    

    # If there are duplicates in df_merged for the same NCBI, combine Database values
    df_merged['Database'] = df_merged.groupby('NCBI')['Database'].transform(lambda x: ','.join(x.dropna().unique()))

    print(df_merged.columns)

    # Drop duplicates and reset index
    df_merged = df_merged.drop_duplicates().reset_index(drop=True)

    print(df_merged.columns)

    # Find NCBIs present in df_databases but not in df_redu
    missing_ncbi = df_databases[~df_databases['NCBI'].isin(df_redu['NCBI'])]

    # Add missing NCBIs to df_merged
    if not missing_ncbi.empty:
        df_merged = pd.concat([df_merged, missing_ncbi], ignore_index=True, sort=False)

    # Fill NaN values in Database column with empty string
    df_merged['Database'] = df_merged['Database'].fillna('')
   
    return df_merged

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
    df = main(args)
