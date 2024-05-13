import pandas as pd
import os
import argparse

def USI2MASST_matchCol(x):
    """Generate match columns based on the USI format."""
    elements = x.split(':')
    mzspec = elements[0]
    msvid = elements[1]
    filename = os.path.splitext(os.path.basename(elements[2]))[0]
    return ':'.join([mzspec, msvid, filename])

def write_empty_output(filename, columns):
    """Writes an empty DataFrame with specified columns to a file."""
    empty_df = pd.DataFrame(columns=columns)
    empty_df.to_csv(filename, sep='\t', index=False)

def main(input_redu_path, input_masst_path, input_sparql_path, input_ncbi_path):
    # Load data with low_memory=False to avoid Dtype warnings
    df_redu = pd.read_csv(input_redu_path, low_memory=False)
    print(f"Loaded REDU data with {df_redu.shape[0]} rows.")

    df_masst = pd.read_csv(input_masst_path, low_memory=False)
    print(f"Loaded MASST data with {df_masst.shape[0]} rows.")

    if not df_masst.empty:
        df_masst['match_col'] = df_masst['USI'].apply(USI2MASST_matchCol)
        print("Added 'match_col' based on 'USI'.")

        if 'Cosine' in df_masst.columns and 'Matching Peaks' in df_masst.columns:
            df_masst['max_cosine'] = df_masst.groupby('match_col')['Cosine'].transform('max')
            max_indices = df_masst.groupby('match_col')['Cosine'].idxmax()
            df_masst = df_masst.loc[max_indices]

            df_merged = pd.merge(df_masst, df_redu, on='match_col', how='inner')
            print(f"Merged MASST with REDU data with {df_merged.shape[0]} rows.")

            df_merged['n_detected'] = df_merged['Cosine'].notna().sum()
            df_merged['n_samples'] = len(df_merged)
            df_merged['p_detected'] = 100 * df_merged['n_detected'] / df_merged['n_samples']
            df_merged.to_csv('masst_by_ncbi_output.tsv', sep='\t', index=False)
            print("MASST results saved to 'masst_by_ncbi_output.tsv'.")
        else:
            print("Required columns 'Cosine' or 'Matching Peaks' are missing in MASST data.")
            write_empty_output('masst_by_ncbi_output.tsv', ['match_col', 'Cosine', 'Matching Peaks', 'n_detected', 'n_samples', 'p_detected', 'ID'])
    else:
        write_empty_output('masst_by_ncbi_output.tsv', ['match_col', 'Cosine', 'Matching Peaks', 'n_detected', 'n_samples', 'p_detected', 'ID'])

    df_sparql = pd.read_csv(input_sparql_path, low_memory=False)
    print(df_sparql.head)
    if len(df_sparql) > 0:
        df_sparql['ID'] = df_sparql['ncbiTaxonomyID'].astype(str)
        print(f"{len(df_sparql['ID'].unique())} NCBI ids reported in Wikidata results!")
        df_sparql = df_sparql[df_sparql['ID'].isin(df_redu['ID'].unique())]
        print(f"{len(df_sparql['ID'].unique())} NCBI ids present after ReDU matching!")
        df_sparql['present'] = True
        df_sparql.to_csv('sparql_by_ncbi_output.tsv', sep='\t', index=False)
        print("SPARQL results saved to 'sparql_by_ncbi_output.tsv'.")
    else:
        write_empty_output('sparql_by_ncbi_output.tsv', ['ncbiTaxonomyID', 'ID', 'present'])

    df_ncbi = pd.read_csv(input_ncbi_path, low_memory=False)
    if len(df_ncbi) > 0:
        df_ncbi['ID'] = df_ncbi['ID'].astype(str)
        df_ncbi = df_ncbi[df_ncbi['ID'].isin(df_redu['ID'].unique())]
        df_ncbi['present'] = True
        df_ncbi.to_csv('ncbi_by_ncbi_output.tsv', sep='\t', index=False)
        print("NCBI results saved to 'ncbi_by_ncbi_output.tsv'.")
    else:
        write_empty_output('ncbi_by_ncbi_output.tsv', ['ID', 'present'])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process NCBI data for matching and output results.")
    parser.add_argument("input_redu_path", help="Path to the REDU TSV file.")
    parser.add_argument("input_masst_path", help="Path to the MASST results TSV file.")
    parser.add_argument("input_sparql_path", help="Path to the SPARQL results TSV file.")
    parser.add_argument("input_ncbi_path", help="Path to the NCBI results TSV file.")
    args = parser.parse_args()

    main(args.input_redu_path, args.input_masst_path, args.input_sparql_path, args.input_ncbi_path)
