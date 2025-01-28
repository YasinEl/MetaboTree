import pandas as pd
import argparse

def load_data(redu_path, sparql_path):

    print(sparql_path)
    df_redu = pd.read_csv(redu_path, sep='\t', low_memory=False)
    dt_sparql = pd.read_csv(sparql_path)

    # Convert NCBI column to integer and strip any whitespace
    df_redu['NCBI'] = df_redu['NCBI'].astype(str).str.strip().astype(int)

    # Collect existing NCBI IDs from df_redu
    existing_ncbi_ids = set(df_redu['NCBI'])
    print(f"Existing NCBI IDs in df_redu (sample): {list(existing_ncbi_ids)[:10]}")

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
        dt_sparql = pd.DataFrame({'NCBI': [], 'Database': []})

    # Concatenate dataframes in case we dont just have wikidata
    df_databases = dt_sparql

    # Merge df_redu with df_databases to get the Database column
    df_merged = pd.merge(df_redu, df_databases, on='NCBI', how='left', suffixes=('', '_from_db'))
    

    # If there are duplicates in df_merged for the same NCBI, combine Database values
    df_merged['Database'] = df_merged.groupby('NCBI')['Database'].transform(lambda x: ','.join(x.dropna().unique()))

    # Drop duplicates and reset index
    df_merged = df_merged.drop_duplicates().reset_index(drop=True)

    # Find NCBIs present in df_databases but not in df_redu
    missing_ncbi = df_databases[~df_databases['NCBI'].isin(df_redu['NCBI'])]

    # Add missing NCBIs to df_merged
    if not missing_ncbi.empty:
        df_merged = pd.concat([df_merged, missing_ncbi], ignore_index=True, sort=False)

    # Fill NaN values in Database column with empty string
    df_merged['Database'] = df_merged['Database'].fillna('')
   
    return df_merged

def main(args):
    df_redu = load_data(args.redu_path, args.sparql_path)
    df_redu.to_csv('all_organism_table.csv', index=False)
    print(f"Data saved to all_organism_table.csv")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Extend df_redu with additional NCBI IDs from CSV files.")
    parser.add_argument('redu_path', help='Path to the REDU CSV file')
    parser.add_argument('sparql_path', help='Path to the SPARQL CSV file')
    
    args = parser.parse_args()
    df = main(args)
