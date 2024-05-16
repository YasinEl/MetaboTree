import argparse
import pandas as pd
import numpy as np
import requests
import logging
import json
import time



def main(file_path, ncbi_to_ToL_file_path):
    df = pd.read_csv(file_path, low_memory=False)
    df_ncbi_to_ToL = pd.read_csv(ncbi_to_ToL_file_path, low_memory=False)\

    print("Data frame loaded successfully.")

    # Filter to get unique NCBI IDs that need organism names fetched
    df['NCBI'] = df['NCBI'].astype(int)
    df_ncbi_to_ToL['NCBI'] = df_ncbi_to_ToL['NCBI'].astype(int)


    df = df.merge(df_ncbi_to_ToL, on='NCBI',  how='left')   


    # Save the modified DataFrame to a CSV file
    df.to_csv('modified_tree_of_life.csv', index=False)
    print("Modified DataFrame saved to modified_tree_of_life.csv")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch and match names using the Open Tree of Life API.")
    parser.add_argument('file_path', type=str, help="Path to the redu table file.")
    parser.add_argument('ncbi_to_ToL_file_path', type=str, help="Path to the ncbi_to_ToL_file.")

    
    args = parser.parse_args()


    main(args.file_path, args.ncbi_to_ToL_file_path)
