import argparse
import time
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolHash
import psycopg2
#from tqdm import tqdm


def parse_arguments():
    parser = argparse.ArgumentParser(description="Find matching SMILES and export masst_results_table entries")
    parser.add_argument("smiles", help="SMILES string to search")
    parser.add_argument("--matching_peaks", type=int, help="Matching peaks")
    parser.add_argument("--masst_now_path", default = '', help="Intput novel MASST table")
    parser.add_argument("--output", default="results.tsv", help="Output TSV file")
    return parser.parse_args()


def check_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    return mol


def fetch_and_match_smiles(conn, target_smiles):
    # Check and convert the target SMILES to a molecule object
    target_mol = Chem.MolFromSmiles(target_smiles)
    if target_mol is None:
        raise ValueError(f"Invalid SMILES: {target_smiles}")
    #target_inchi_key = Chem.MolToInchiKey(target_mol)
    target_inchi_key = Chem.MolToInchiKey(target_mol).split('-')[0]

    with conn.cursor() as cur:
        print("Fetching SMILES and associated columns from gnps_library...")
        start_time = time.time()

        # Fetch all relevant columns in one step
        cur.execute("""
            SELECT inchikey_smiles, spectrum_id, collision_energy, adduct, ms_manufacturer, 
                   ms_mass_analyzer, gnps_library_membership, scan 
            FROM gnps_library
        """)
        
        rows = cur.fetchall()
        print(f"Fetched {len(rows)} rows in {time.time() - start_time:.2f} seconds.")
        
        # Convert the fetched rows to a DataFrame
        col_names = [desc[0] for desc in cur.description]
        df = pd.DataFrame(rows, columns=col_names)

        print("Matching SMILES (ignoring stereochemistry)...")
        start_matching_time = time.time()

        # Only keep rows with valid SMILES and perform the structure match (exact match based on InChIKey)
        df['inchi_key_first_block'] = df['inchikey_smiles'].apply(lambda inchi: inchi.split('-')[0] if inchi else None)
        df_matched = df[df['inchi_key_first_block'] == target_inchi_key]


        print(f"Found {len(df_matched)} matching SMILES in {time.time() - start_matching_time:.2f} seconds.")
        
        # If no matching SMILES were found, return an empty list
        if df_matched.empty:
            print("No matching SMILES found.")
            return []


        df_matched[['collision_energy', 'adduct', 'ms_manufacturer', 'ms_mass_analyzer', 'gnps_library_membership']] = df_matched[
            ['collision_energy', 'adduct', 'ms_manufacturer', 'ms_mass_analyzer', 'gnps_library_membership']
        ].fillna('unknown')

        # Group by the required columns and limit to at most 8 rows per group
        df_matched['row_num'] = df_matched.groupby(
            ['collision_energy', 'adduct', 'ms_manufacturer', 'ms_mass_analyzer', 'gnps_library_membership']
        ).cumcount() + 1

        
        # Keep only the first 8 rows per group
        df_limited = df_matched[df_matched['row_num'] <= 8].drop(columns=['row_num', 'inchi_key_first_block'])

        print(f"Limited to {len(df_limited)} rows after grouping.")

        print(df_limited['spectrum_id'].unique().tolist())

        
        # Return the unique spectrum_ids
        return df_limited['spectrum_id'].unique().tolist()



def fetch_masst_results(conn, spectrum_ids, matching_peaks):
    with conn.cursor() as cur:
        print("Fetching masst_results_table entries for matching spectrum_ids...")
        start_time = time.time()
        
        # Fetch rows by spectrum_id from masst_results_table
        query = """
        SELECT *
        FROM masst_results_table
        WHERE spectrum_id = ANY(%s)
        """
        cur.execute(query, (spectrum_ids,))
        rows = cur.fetchall()
        col_names = [desc[0] for desc in cur.description]
        print(f"Fetched {len(rows)} rows in {time.time() - start_time:.2f} seconds.")
        
        # Extract mri_ids to fetch corresponding mri values from unique_mri
        mri_ids = list({row[col_names.index("mri_id")] for row in rows if row[col_names.index("mri_id")] is not None})

        
         # Fetch mri values from unique_mri table in chunks
        chunk_size = 1000  # Adjust this chunk size as needed
        mri_map = {}
        for i in range(0, len(mri_ids), chunk_size):
            chunk = mri_ids[i:i + chunk_size]
            query_mri = """
            SELECT id, mri
            FROM unique_mri
            WHERE id = ANY(%s)
            """
            cur.execute(query_mri, (chunk,))
            mri_map.update({mri_id: mri for mri_id, mri in cur.fetchall()})

        # Add mri values to fetched rows in-memory
        rows = [
            dict(zip(col_names + ["mri"], row + (mri_map.get(row[col_names.index("mri_id")]),)))
            for row in rows
        ]
        
        # Filter in-memory by matching_peaks
        print(f"Filtering rows with matching_peaks >= {matching_peaks}...")
        start_filter_time = time.time()
        
        df = pd.DataFrame(rows, columns=col_names + ["mri"])

        # Convert 'matching_peaks' column to numeric, coercing errors to NaN
        df['matching_peaks'] = pd.to_numeric(df['matching_peaks'], errors='coerce')

        df = df.dropna(subset=['matching_peaks'])

        # Ensure matching_peaks is a float for comparison
        df['matching_peaks'] = df['matching_peaks'].astype(float)

        # Filter rows where 'matching_peaks' is greater than or equal to the threshold
        filtered_df = df[df['matching_peaks'] >= matching_peaks]

        # Convert the filtered DataFrame back to a list of rows (if needed)
        filtered_rows = filtered_df.to_dict(orient='records')
        print(f"Filtered down to {len(filtered_rows)} rows in {time.time() - start_filter_time:.2f} seconds.")
        
    return filtered_rows




def save_results_to_tsv(rows, output_file, masst_novel):
    print(f"Saving results to {output_file}...")
    start_time = time.time()
    df = pd.DataFrame(rows)

    df.rename(columns={'mri': 'USI'}, inplace=True)
    df.rename(columns={'cosine': 'Cosine'}, inplace=True) 
    df.rename(columns={'matching_peaks': 'Matching Peaks'}, inplace=True)

    # Add the novel MASST table
    if masst_novel != '':
        df_novel = pd.read_csv(masst_novel)
        if len(df_novel) > 0:
            df = pd.concat([df, df_novel], ignore_index=True)


    df.to_csv(output_file, index=False)
    print(f"Saved results in {time.time() - start_time:.2f} seconds.")


def main():
    args = parse_arguments()


    # Database connection
   # os.environ['PGPASSWORD'] = '9421'
    conn = psycopg2.connect(
        dbname="masst_records",
        user="yasel",
        password="9421",
        host="localhost",
        port=5432
    )

    try:

        # Step 1: Find matching SMILES (ignoring stereochemistry)
        matching_lib_ids = fetch_and_match_smiles(conn, args.smiles)

        if not matching_lib_ids:
            print("No matching SMILES found.")
            return

        # Step 3: Fetch masst_results_table entries for the matching lib_ids
        masst_results = fetch_masst_results(conn, matching_lib_ids, args.matching_peaks)

        # Step 4: Save results to TSV
        save_results_to_tsv(masst_results, args.output, args.masst_now_path)

    finally:
        conn.close()


if __name__ == "__main__":
    main()
