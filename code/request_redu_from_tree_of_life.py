import argparse
import pandas as pd
import numpy as np
import requests
import logging
import json
import time



def main(file_path, linage_path, ncbi_to_ToL_file_path):
    df = pd.read_csv(file_path, low_memory=False)
    df_linage = pd.read_csv(linage_path)
    df_ncbi_to_ToL = pd.read_csv(ncbi_to_ToL_file_path, low_memory=False)\

    print("Data frame loaded successfully.")

    # Filter to get unique NCBI IDs that need organism names fetched
    df['NCBI'] = df['NCBI'].astype(int)
    df_linage['NCBI'] = df_linage['NCBI'].astype(int)
    df_ncbi_to_ToL['NCBI'] = df_ncbi_to_ToL['NCBI'].astype(int)
    df_ncbi_to_ToL['uid'] = df_ncbi_to_ToL['uid'].astype(int)


    df = df.merge(df_ncbi_to_ToL, on='NCBI',  how='left')   
    df_linage = df_linage.merge(df_ncbi_to_ToL, on='NCBI', how='left', suffixes=('', '_dup'))

    # Fill missing values in original columns with values from the duplicated ones
    for col in df_linage.columns:
        if col.endswith('_dup'):
            original_col = col[:-4]
            df_linage[original_col] = df_linage[original_col].fillna(df_linage[col])

    # Drop the duplicated columns
    df_linage = df_linage.drop(columns=[col for col in df_linage.columns if col.endswith('_dup')])


    # Save the modified DataFrame to a CSV file
    df.to_csv('modified_tree_of_life.csv', index=False)
    df_linage.to_csv('linage_table.csv', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch and match names using the Open Tree of Life API.")
    parser.add_argument('redu_path', type=str, help="Path to the redu table file.")
    parser.add_argument('linage_path', type=str, help="Path to the redu table file.")
    parser.add_argument('ncbi_to_ToL_file_path', type=str, help="Path to the ncbi_to_ToL_file.")

    
    args = parser.parse_args()


    main(args.redu_path, args.linage_path, args.ncbi_to_ToL_file_path)
