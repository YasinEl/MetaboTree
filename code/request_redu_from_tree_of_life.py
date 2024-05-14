import argparse
import pandas as pd
import numpy as np
import requests
import logging
from Bio import Entrez
import json
import time

# Configure logging
logging.basicConfig(level=logging.DEBUG, filename='treeoflife_detailed.log', filemode='w',
                    format='%(asctime)s - %(levelname)s - %(message)s')

def robust_request(url, data, headers, retries=5, backoff_factor=0.3):
    """Make a robust HTTP POST request with retries and exponential backoff."""
    logging.debug("Preparing to send request.")
    logging.debug(f"URL: {url}")
    logging.debug(f"Headers: {headers}")
    logging.debug(f"Payload: {data}")

    for attempt in range(retries):
        try:
            response = requests.post(url, data=data, headers=headers)
            logging.debug(f"Attempt {attempt + 1}: HTTP status code {response.status_code}")
            if response.status_code == 200:
                logging.info("Request successful.")
                return response.json()
            else:
                logging.error(f"Failed to fetch data: {response.status_code}, {response.text}")
                response.raise_for_status()  # To handle non-200 responses and raise an HTTPError
        except requests.exceptions.RequestException as e:
            logging.error(f"Request failed on attempt {attempt + 1}: {str(e)}")
            time.sleep(backoff_factor * (2 ** attempt))  # Exponential backoff

    logging.error("All retries failed.")
    return None


def match_names_to_ott(names):
    """Match names using the Open Tree of Life API."""
    url = "https://api.opentreeoflife.org/v3/tnrs/match_names"
    headers = {'Content-Type': 'application/json'}

    if isinstance(names, (pd.Series, pd.Index, np.ndarray)):
        names = names.tolist()
        logging.debug(f"Converted names to list: {names}")

    data = json.dumps({'names': names, 'do_approximate_matching': False})
    logging.debug(f"Sending data: {data}")

    response = robust_request(url, data, headers)
    if response:
        logging.info("Data retrieval successful.")
        return response
    else:
        logging.warning("Data retrieval returned None.")
        return None

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
