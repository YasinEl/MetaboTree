import pandas as pd
import argparse
import os

def load_data(redu_path):

    df_redu = pd.read_csv(redu_path, sep='\t', low_memory=False)

    # List all CSV files in the current directory starting with 'sparql_'
    sparql_path = [f for f in os.listdir('.') if f.startswith('sparql_') and f.endswith('.csv')]
    
    print(f"Found SPARQL files: {sparql_path}")

    # if sparql_path is a list of paths, load all of them and add column with filebasename after removing sparql_ prefix and .csv suffix
    if isinstance(sparql_path, list):
        df_list = []
        for path in sparql_path:
            df = pd.read_csv(path)
            basename = os.path.basename(path).replace('.csv', '')
            df['Source'] = basename
            df_list.append(df)
        dt_sparql = pd.concat(df_list, ignore_index=True)
    else:
        dt_sparql = pd.read_csv(sparql_path)

    dt_sparql['NCBI'] = dt_sparql['ncbi'].astype(str).str.strip().astype(int)
    dt_sparql['inchikey'] = dt_sparql['inchikey'].str.split('-').str[0]

    dt_sparql = dt_sparql[['Source', 'inchikey', 'NCBI']]

    dt_sparql = dt_sparql.drop_duplicates()

    # Replace spaces with underscores in the 'Source' column
    dt_sparql['Source'] = dt_sparql['Source'].str.replace(' ', '_')

    # Count inchikey by Source and NCBI
    dt_sparql = dt_sparql.groupby(['Source', 'NCBI']).size().reset_index(name='inchikey_count')
    
    # Use Source as columns and NCBI as index and inchikey_count as values
    dt_sparql = dt_sparql.pivot(index='NCBI', columns='Source', values='inchikey_count').fillna(0).astype(int)
    dt_sparql = dt_sparql.reset_index()

    # Convert NCBI column to integer and strip any whitespace
    df_redu['NCBI'] = df_redu['NCBI'].astype(str).str.strip().astype(int)

    # Merge df_redu with df_databases to get the Database column
    df_merged = pd.merge(df_redu, dt_sparql, on='NCBI', how='left', suffixes=('', '_from_db'))


    sparql_columns = [col for col in df_merged.columns if col.startswith('Sparql_')]
    df_merged[sparql_columns] = df_merged[sparql_columns].fillna(0).astype(int)
 
    return df_merged

def main(args):
    df_redu = load_data(args.redu_path)
    df_redu.to_csv('all_organism_table.csv', index=False)
    print(f"Data saved to all_organism_table.csv")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Extend df_redu with additional NCBI IDs from CSV files.")
    parser.add_argument('redu_path', help='Path to the REDU CSV file')
    
    args = parser.parse_args()
    df = main(args)
