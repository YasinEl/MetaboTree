import argparse
import requests
from rdkit import Chem
import time
import sys
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry

def fetch_related_cids(cid):
    """Fetch all CIDs related by connectivity to the given CID."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/cids/JSON?cids_type=same_isotopes"
    print("fetching related CIDs")
    attempts = 0
    while attempts < 5:
        try:
            response = requests.get(url, timeout=8)
            print('got something')
            response.raise_for_status()  # Raises an HTTPError for bad responses
            data = response.json()
            print(data)
            return data['IdentifierList']['CID']
        except requests.exceptions.RequestException as e:
            attempts += 1
            print(f"Attempt {attempts} failed: {e}")
            if attempts == 2:
                print("Max attempts reached. Saving input CID.")
                return [cid]


def fetch_related_cids(cid):
    """Fetch all CIDs related by connectivity to the given CID."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/cids/JSON?cids_type=same_isotopes"
    print("fetching related CIDs")

    session = requests.Session()
    retries = Retry(total=5, backoff_factor=1, status_forcelist=[502, 503, 504])
    session.mount('http://', HTTPAdapter(max_retries=retries))
    session.mount('https://', HTTPAdapter(max_retries=retries))

    attempts = 0
    while attempts < 5:
        try:
            response = session.get(url, timeout=8)
            print('got something')
            response.raise_for_status()  # Raises an HTTPError for bad responses
            data = response.json()
            print(data)
            return data['IdentifierList']['CID']
        except requests.exceptions.RequestException as e:
            attempts += 1
            print(f"Attempt {attempts} failed: {e}")
            if attempts == 2:
                print("Max attempts reached. Saving input CID.")
                return [cid]
            

def fetch_main_name(cid):
    """Fetch the main name for a given CID."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/description/JSON"
    try:
        response = requests.get(url)
        response.raise_for_status()
        data = response.json()
        main_name = data['InformationList']['Information'][0]['Title']
        return main_name
    except requests.exceptions.RequestException as e:
        print(f"Error fetching main name for CID {cid}: {e}")
        return None


def fetch_main_name(cid):
    """Fetch the main name for a given CID."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/description/JSON"
    
    session = requests.Session()
    retries = Retry(total=5, backoff_factor=1, status_forcelist=[502, 503, 504])
    session.mount('http://', HTTPAdapter(max_retries=retries))
    session.mount('https://', HTTPAdapter(max_retries=retries))
    
    attempts = 0
    while attempts < 5:
        try:
            response = session.get(url, timeout=8)
            response.raise_for_status()
            data = response.json()
            main_name = data['InformationList']['Information'][0]['Title']
            return main_name
        except requests.exceptions.RequestException as e:
            attempts += 1
            print(f"Attempt {attempts} failed: {e}")
            if attempts == 2:
                print("Max attempts reached. Could not fetch main name.")
                return None
            
            
def save_cids_to_tsv(cids, output_file):
    """Save the list of CIDs to a TSV file."""
    with open(output_file, "w") as file:
        for cid in set(cids):
            file.write(f"{cid}\n")

def save_main_name_to_tsv(cid, output_file):
    """Save the main name of the given CID to a TSV file."""
    main_name = fetch_main_name(cid)
    with open(output_file, "w") as file:
        if main_name:
            file.write(f"{main_name}\n")

def main(cid):
    related_cids = fetch_related_cids(cid)
    print('fetch_related_cids -> done')
    save_cids_to_tsv(related_cids, "stereoisomers.tsv")
    print('save_cids_to_tsv -> done')
    save_main_name_to_tsv(cid, "main_name.tsv")
    print('save_main_name_to_tsv -> done')
    print(f"Saved {len(related_cids)} CIDs to stereoisomers.tsv")
    print(f"Saved the main name of CID {cid} to main_name.tsv")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch stereoisomer CIDs and their main name.")
    parser.add_argument('cid', type=int, help='PubChem Compound ID (CID) to process.')
    args = parser.parse_args()

    main(args.cid)
