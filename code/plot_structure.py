import argparse
import requests
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from PIL import Image
import csv
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry

def get_smiles_and_name(pubchem_id):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{pubchem_id}/property/CanonicalSMILES,IUPACName/JSON"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        if 'PropertyTable' in data and 'Properties' in data['PropertyTable']:
            properties = data['PropertyTable']['Properties'][0]
            smiles = properties.get('CanonicalSMILES')
            name = properties.get('IUPACName')
            return smiles, name
    return None, None




def get_smiles_and_name(pubchem_id):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{pubchem_id}/property/CanonicalSMILES,IUPACName/JSON"
    
    session = requests.Session()
    retries = Retry(total=5, backoff_factor=1, status_forcelist=[502, 503, 504])
    session.mount('http://', HTTPAdapter(max_retries=retries))
    session.mount('https://', HTTPAdapter(max_retries=retries))
    
    attempts = 0
    while attempts < 5:
        try:
            response = session.get(url, timeout=8)
            response.raise_for_status()
            if response.status_code == 200:
                data = response.json()
                if 'PropertyTable' in data and 'Properties' in data['PropertyTable']:
                    properties = data['PropertyTable']['Properties'][0]
                    smiles = properties.get('CanonicalSMILES')
                    name = properties.get('IUPACName')
                    return smiles, name
            return None, None
        except requests.exceptions.RequestException as e:
            attempts += 1
            print(f"Attempt {attempts} failed: {e}")
            if attempts == 2:
                print("Max attempts reached. Could not fetch SMILES and name.")
                return None, None
            

def make_transparent(input_file, output_file):
    img = Image.open(input_file)
    img = img.convert("RGBA")
    
    datas = img.getdata()

    new_data = []
    for item in datas:
        # Change all white (also shades of whites)
        # pixels to transparent
        if item[0] > 200 and item[1] > 200 and item[2] > 200:
            new_data.append((255, 255, 255, 0))
        else:
            new_data.append(item)

    img.putdata(new_data)
    img.save(output_file, "PNG")

def main():
    parser = argparse.ArgumentParser(description='Fetch SMILES notation and name for a given PubChem CID, save it as a plotted molecule in PNG format, and store the name in a TSV file.')
    parser.add_argument('cid', type=int, help='PubChem CID')
    
    args = parser.parse_args()
    
    smiles, name = get_smiles_and_name(args.cid)
    if smiles:
        mol = Chem.MolFromSmiles(smiles)
        if mol:

            print('working on structure plot')
            # Draw the molecule with RDKit and save it as a high-resolution PNG
            d2d = Draw.MolDraw2DCairo(2500, 2500)  # Increase resolution as needed
            d2d.drawOptions().bgColor = (1, 1, 1, 0)  # Set the background color to white
            d2d.DrawMolecule(mol)
            d2d.FinishDrawing()
            print('FinishDrawing -> finshed')
            png = d2d.GetDrawingText()
            temp_output_file = 'temp_' + 'molecule.png'
            with open(temp_output_file, 'wb') as f:
                f.write(png)
            print('initial structure plot -> finshed')
            make_transparent(temp_output_file, 'molecule.png')
            print(f'Transparaent molecule image for CID {args.cid} saved to {'molecule.png'}')
            
            # Save the name to a TSV file
            with open('mol_name.tsv', 'a', newline='') as tsvfile:
                writer = csv.writer(tsvfile, delimiter='\t')
                writer.writerow([args.cid, name])
            print(f'Molecule name for CID {args.cid} saved to mol_name.tsv')
        else:
            print(f'Invalid SMILES notation: {smiles}')
    else:
        print(f'No SMILES notation or name found for CID {args.cid}')

if __name__ == '__main__':
    main()
