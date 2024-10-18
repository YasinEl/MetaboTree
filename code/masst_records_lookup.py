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
    parser.add_argument("--output", default="results.tsv", help="Output TSV file")
    return parser.parse_args()


def check_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    return mol


def fetch_distinct_smiles(conn):
    with conn.cursor() as cur:
        print("Fetching distinct SMILES and associated spectrum_ids from gnps_library...")
        start_time = time.time()
        cur.execute("""
            SELECT smiles, array_agg(spectrum_id) AS spectrum_ids 
            FROM (
                SELECT smiles, spectrum_id, 
                ROW_NUMBER() OVER (PARTITION BY collision_energy, Adduct, ms_manufacturer, ms_mass_analyzer, gnps_library_membership ORDER BY scan) AS row_num
                FROM gnps_library
            ) subquery
            WHERE row_num <= 8
            GROUP BY smiles
        """)
        rows = cur.fetchall()
        print(f"Fetched {len(rows)} distinct SMILES in {time.time() - start_time:.2f} seconds.")
    return rows






def match_smiles(target_mol, smiles_list):
    print("Checking for matching SMILES (ignoring stereochemistry) using InChIKey...")
    start_time = time.time()
    matching_spectrum_ids = []

    # Generate InChIKey for the target molecule
    target_inchi = Chem.MolToInchiKey(target_mol).split('-')[0]  # Take the first part of InChIKey for non-stereochemical match

    for smiles, spectrum_ids in smiles_list:  # Each row has SMILES and associated spectrum_ids
        if smiles is None:
            print("Skipping NoneType SMILES.")
            continue

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"Skipping invalid SMILES: {smiles}")
            continue

        # Generate InChIKey for the current molecule
        mol_inchi = Chem.MolToInchiKey(mol).split('-')[0]  # Take the first part for non-stereochemical match

        # Check if InChIKeys match (ignoring stereochemistry)
        if target_inchi == mol_inchi:
            matching_spectrum_ids.extend(spectrum_ids)  # Add all associated spectrum_ids for this SMILES

    print(f"Found {len(matching_spectrum_ids)} matching spectrum_ids in {time.time() - start_time:.2f} seconds.")
    return matching_spectrum_ids





def fetch_masst_results(conn, spectrum_ids, matching_peaks):
    with conn.cursor() as cur:
        print("Fetching masst_results_table entries for matching spectrum_ids...")
        start_time = time.time()

        print(len(spectrum_ids))

        print(spectrum_ids)
        
        # Fetch rows by spectrum_id
        query = """
        SELECT *
        FROM masst_results_table
        WHERE spectrum_id = ANY(%s)
        """
        cur.execute(query, (spectrum_ids,))
        rows = cur.fetchall()

        col_names = [desc[0] for desc in cur.description]
        print(f"Fetched {len(rows)} rows in {time.time() - start_time:.2f} seconds.")
        
        # Filter in-memory by matching_peaks
        print(f"Filtering rows with matching_peaks >= {matching_peaks}...")
        start_filter_time = time.time()
        
        df = pd.DataFrame(rows, columns=col_names)


        # Convert 'matching_peaks' column to numeric, coercing errors to NaN
        df['matching_peaks'] = pd.to_numeric(df['matching_peaks'], errors='coerce')

        df = df.dropna(subset=['matching_peaks'])

        # Ensure matching_peaks is an integer (or float) for comparison
        df['matching_peaks'] = df['matching_peaks'].astype(float)

        # Filter rows where 'matching_peaks' is greater than or equal to the threshold
        filtered_df = df[df['matching_peaks'] >= matching_peaks]


        # Convert the filtered DataFrame back to a list of rows (if needed)
        filtered_rows = filtered_df.to_dict(orient='records')
        print(f"Filtered down to {len(filtered_rows)} rows in {time.time() - start_filter_time:.2f} seconds.")
        
    return filtered_rows



def save_results_to_tsv(rows, output_file):
    print(f"Saving results to {output_file}...")
    start_time = time.time()
    df = pd.DataFrame(rows)

    df.rename(columns={'mri': 'USI'}, inplace=True)
    df.rename(columns={'cosine': 'Cosine'}, inplace=True) 
    df.rename(columns={'matching_peaks': 'Matching Peaks'}, inplace=True)


    df.to_csv(output_file, index=False)
    print(f"Saved results in {time.time() - start_time:.2f} seconds.")


def main():
    args = parse_arguments()

    # Check SMILES validity
    try:
        target_mol = check_smiles(args.smiles)
    except ValueError as e:
        print(e)
        return

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
        # Step 1: Fetch distinct SMILES from gnps_library
        smiles_list = fetch_distinct_smiles(conn)

        # Step 2: Find matching SMILES (ignoring stereochemistry)
        matching_lib_ids = match_smiles(target_mol, smiles_list)

        if not matching_lib_ids:
            print("No matching SMILES found.")
            return

        # Step 3: Fetch masst_results_table entries for the matching lib_ids
        masst_results = fetch_masst_results(conn, matching_lib_ids, args.matching_peaks)

        # Step 4: Save results to TSV
        save_results_to_tsv(masst_results, args.output)

    finally:
        conn.close()


if __name__ == "__main__":
    main()
