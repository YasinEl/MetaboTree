import argparse
import requests
from rdkit import Chem
from rdkit.Chem import Draw

def get_smiles(pubchem_id):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{pubchem_id}/property/CanonicalSMILES/JSON"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        if 'PropertyTable' in data and 'Properties' in data['PropertyTable']:
            smiles = data['PropertyTable']['Properties'][0].get('CanonicalSMILES')
            return smiles
    return None

def main():
    parser = argparse.ArgumentParser(description='Fetch SMILES notation for a given PubChem CID and save it as a plotted molecule in SVG format.')
    parser.add_argument('cid', type=int, help='PubChem CID')
    parser.add_argument('output_file', type=str, help='Output SVG file')
    
    args = parser.parse_args()
    
    smiles = get_smiles(args.cid)
    if smiles:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            # Save the molecule as SVG
            svg = Draw.MolsToGridImage([mol], subImgSize=(300, 300), useSVG=True)
            with open(args.output_file, 'w') as f:
                f.write(svg)
            print(f'Molecule image for CID {args.cid} saved to {args.output_file}')
        else:
            print(f'Invalid SMILES notation: {smiles}')
    else:
        print(f'No SMILES notation found for CID {args.cid}')

if __name__ == '__main__':
    main()
