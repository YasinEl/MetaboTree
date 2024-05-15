import argparse
import requests
from rdkit import Chem
import time
import sys

def fetch_related_cids(cid):
    """Fetch all CIDs related by connectivity to the given CID."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/cids/JSON?cids_type=same_isotopes"
    print("fetching related CIDs")
    attempts = 0
    while attempts < 2:
        try:
            response = requests.get(url, timeout=60)
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

def filter_cids_remove_isotopes(related_cids):
    """Filter out CIDs that contain isotopes in their structure."""
    non_isotope_cids = []
    i=1
    for cid in related_cids:
        print(cid)
        time.sleep(0.2)  # Delay to comply with the API rate limit of 5 requests per second
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/IsomericSMILES/JSON"
        try:
            response = requests.get(url)
            response.raise_for_status()
            data = response.json()
            smiles = data['PropertyTable']['Properties'][0]['IsomericSMILES']
            mol = Chem.MolFromSmiles(smiles)
            if not any(atom.GetIsotope() != 0 for atom in mol.GetAtoms()):
                non_isotope_cids.append(cid)
        except requests.exceptions.RequestException as e:
            print(f"Error processing CID {cid}: {e}")
            continue
        print(f"{i}/{len(related_cids)}")
        i=i+1
    return non_isotope_cids

def save_cids_to_tsv(cids, output_file):
    """Save the list of CIDs to a TSV file."""
    with open(output_file, "w") as file:
        for cid in set(cids):
            file.write(f"{cid}\n")

def main(cid):
    related_cids = fetch_related_cids(cid)
    #non_isotope_cids = filter_cids_remove_isotopes([cid] + related_cids)
    save_cids_to_tsv(related_cids, "stereoisomers.tsv")
    print(f"Saved {len(related_cids)} CIDs to stereoisomers.tsv")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch stereoisomer CIDs excluding isotopes.")
    parser.add_argument('cid', type=int, help='PubChem Compound ID (CID) to process.')
    args = parser.parse_args()

    main(args.cid)