import argparse
import os
import time
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolHash
import psycopg2
from tqdm import tqdm
from formula_validation.Formula import Formula 
from rdkit.Chem import rdMolDescriptors
import re


def parse_arguments():
    parser = argparse.ArgumentParser(description="Find matching SMILES and export masst_results_table entries")
    parser.add_argument("--smiles", default='', help="SMILES string to search")
    parser.add_argument("--smiles_name", default='', help="SMILES name")
    parser.add_argument("--smiles_type", default='smiles', help="SMILES type")
    parser.add_argument("--match_type", default='exact', help="exact or substructure")
    parser.add_argument("--structure_file", default='', help="SMILES tsv file path")
    parser.add_argument("--matching_peaks", type=int, help="Matching peaks")
    parser.add_argument("--masst_now_path", default = '', help="Intput novel MASST table")
    parser.add_argument("--output", default="results.tsv", help="Output TSV file")
    return parser.parse_args()


def check_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    return mol


def fetch_and_match_smiles(conn, target_smiles, match_type='exact', smiles_name = 'only', smiles_type = 'smiles', formula_base = 'any', element_diff = 'any'):
    # Check and convert the target SMILES to a molecule object

    if smiles_type == 'smiles':
        target_mol = Chem.MolFromSmiles(target_smiles)
    elif smiles_type == 'smarts':
        target_mol = Chem.MolFromSmarts(target_smiles)

    if target_mol is None:
        raise ValueError(f"Invalid SMILES: {target_smiles}")
    target_inchi_key = Chem.MolToInchiKey(target_mol).split('-')[0]
    
    if formula_base != 'any':
        print(f"Formula base: {formula_base}")
        formula_base = Formula.formula_from_str(formula_base)

    if element_diff != 'any':
        print(f"Element diff: {element_diff}")
        element_diff = Formula.formula_from_str(element_diff)

    with conn.cursor() as cur:
        print("Fetching SMILES and associated columns from gnps_library...")
        start_time = time.time()
    
        # Fetch all relevant columns in one step
        cur.execute("""
            SELECT inchikey_smiles, smiles, spectrum_id, collision_energy, adduct, ms_manufacturer, 
                   ms_mass_analyzer, gnps_library_membership, scan 
            FROM gnps_library
        """)
        
        rows = cur.fetchall()
        print(f"Fetched {len(rows)} rows in {time.time() - start_time:.2f} seconds.")
        
        # Convert the fetched rows to a DataFrame
        col_names = [desc[0] for desc in cur.description]
        df = pd.DataFrame(rows, columns=col_names)
    
        df['inchi_key_first_block'] = df['inchikey_smiles'].apply(lambda inchi: inchi.split('-')[0] if inchi else None)

        if match_type == 'exact':
            print("Matching SMILES (ignoring stereochemistry)...")
            start_matching_time = time.time()
    
            # Only keep rows with valid SMILES and perform the structure match (exact match based on InChIKey)

            df_matched = df[df['inchi_key_first_block'] == target_inchi_key]
    
            print(f"Found {len(df_matched)} matching SMILES in {time.time() - start_matching_time:.2f} seconds.")
        else:
            print("Performing substructure search...")
            start_matching_time = time.time()
            # Get unique SMILES
            unique_smiles = df['smiles'].dropna().unique()
            print(f"Total unique SMILES to process: {len(unique_smiles)}")
            # Create a dict mapping SMILES to Mol
            smiles_to_mol = {}
            for smiles in unique_smiles:
                mol = Chem.MolFromSmiles(smiles)
                if mol is not None:
                    smiles_to_mol[smiles] = mol
                else:
                    print(f"Invalid SMILES encountered and skipped: {smiles}")

            # Perform substructure matching
            matching_smiles = []
            for smiles, mol in tqdm(smiles_to_mol.items(), desc="Substructure Matching", total=len(smiles_to_mol)):
                if mol.HasSubstructMatch(target_mol):
                    print('substructure matching successfull')


                    # Add formula difference matching
                    if formula_base != 'any':
                        formula_candidate = rdMolDescriptors.CalcMolFormula(mol)
                        formula_candidate = Formula.formula_from_str(formula_candidate)

                        try:
                            formula_diff_here = formula_candidate - formula_base
                        except:
                            continue

                        diff_comparisson = 'wrong'
                        try:
                            diff_comparisson = formula_diff_here - element_diff
                        except:
                            try:
                                diff_comparisson = element_diff - formula_diff_here
                            except:
                                diff_comparisson = 'wrong'



                        match_is = str(diff_comparisson) == '' or (set(re.sub(r'\d', '', str(diff_comparisson))) == {"H"})

                        if match_is == False:
                            continue
                    print(f"Match found: {smiles}")
                    matching_smiles.append(smiles)

            print(f"Found {len(matching_smiles)} matching SMILES in {time.time() - start_matching_time:.2f} seconds via substructure matching.")
            # Filter df to include only rows where 'inchikey_smiles' is in matching_smiles
            df_matched = df[df['smiles'].isin(matching_smiles)]
    
        # If no matching SMILES were found, return an empty list
        if df_matched.empty:
            print("No matching structures found.")
            return []
    
        df_matched[['collision_energy', 'adduct', 'ms_manufacturer', 'ms_mass_analyzer', 'gnps_library_membership']] = df_matched[
            ['collision_energy', 'adduct', 'ms_manufacturer', 'ms_mass_analyzer', 'gnps_library_membership']
        ].fillna('unknown')
    
        # Group by the required columns and limit to at most 8 rows per group
        df_matched['row_num'] = df_matched.groupby(
            ['collision_energy', 'adduct', 'ms_manufacturer', 'ms_mass_analyzer', 'gnps_library_membership', 'inchikey_smiles']
        ).cumcount() + 1
    
        # Keep only the first 8 rows per group
        df_limited = df_matched[df_matched['row_num'] <= 8]

        df_limited['smiles_name'] = smiles_name

    
        print(f"Limited to {len(df_limited)} rows after grouping.")
        print(df_limited['spectrum_id'].unique().tolist())
        
        # Return the unique spectrum_ids
        return df_limited[["spectrum_id", "smiles_name", "inchi_key_first_block"]]




def fetch_masst_results(conn, lib_table, matching_peaks):
    with conn.cursor() as cur:
        print("Fetching masst_results_table entries for matching spectrum_ids...")
        start_time = time.time()

        spectrum_ids = lib_table['spectrum_id'].unique().tolist()
        
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


        filtered_df = pd.merge(filtered_df, lib_table, on='spectrum_id', how='inner')

        # Convert the filtered DataFrame back to a list of rows (if needed)
        filtered_rows = filtered_df.to_dict(orient='records')
        print(f"Filtered down to {len(filtered_rows)} rows in {time.time() - start_filter_time:.2f} seconds.")
        
    return filtered_rows




def save_results_to_tsv(df, output_file, masst_novel):
    print(f"Saving results to {output_file}...")
    start_time = time.time()

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

        # check if the path args.structure_path points to an existing file
        print(args.structure_file)
        if os.path.isfile(args.structure_file):
            #read the tsv file
            df_smiles = pd.read_csv(args.structure_file, sep='\t')
        else:
            #Make dataframe with smiles value from smile sagument
            df_smiles = pd.DataFrame({'structure': [args.smiles],
                                      'match': [args.match_type],
                                      'type': [args.smiles_type],
                                      'name': [args.smiles_name]})
            
        all_results = []
        formula_base = 'any'
        element_diff = 'any'
        #take smiles, type and name in for loop per row
        for index, row in df_smiles.iterrows():
            print(index)
            smiles = row['structure']
            type = row['type']
            name = row['name']
            match = row['match']
            print(name)


            try:

                #Check if key "formula_base" is in the row
                if 'formula_base' in row.keys():
                    formula_base = row['formula_base']
                else:
                    formula_base = 'any'

                if 'element_diff' in row.keys():
                    element_diff = row['element_diff']
                else:
                    element_diff = 'any'
            
            except ValueError as e:
                print(f"Problem with formula_base {str(formula_base)} or element_diff {str(element_diff)}")
                print(e)


            try:
                # Step 1: Find matching SMILES (ignoring stereochemistry)
                matching_lib_ids = fetch_and_match_smiles(conn, smiles_name=name, target_smiles=smiles, match_type=match, smiles_type=type, formula_base=formula_base, element_diff=element_diff)

                if len(matching_lib_ids) == 0:
                    print("No matching SMILES found.")
                    continue  
                
                print(f"Found {len(matching_lib_ids)} matching SMILES.")
                # Step 2: Fetch masst_results_table entries for the matching lib_ids
                masst_results = fetch_masst_results(conn, matching_lib_ids, args.matching_peaks)

                print(f"Found {len(masst_results)} matching masst_results_table entries.")

                all_results.append(pd.DataFrame(masst_results))
            except ValueError as e:
                print(e)

        masst_results_final = pd.concat(all_results, ignore_index=True)

        #print columns
        print(masst_results_final.columns)

        # If we have multiple matches to the same spectrum keep the moelcule with the best match
        masst_results_final = masst_results_final.sort_values(by=['cosine', 'matching_peaks'], ascending=[False, False])
        masst_results_final = masst_results_final.drop_duplicates(subset=['mri_id', 'scan_id'])


        # Step 4: Save results to TSV
        save_results_to_tsv(masst_results_final, args.output, args.masst_now_path)

    finally:
        conn.close()


if __name__ == "__main__":
    main()
