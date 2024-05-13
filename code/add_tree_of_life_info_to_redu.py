import argparse
import pandas as pd
import json

def load_matches(filename):
    """Load matches from a JSON file."""
    with open(filename, 'r') as file:
        data = json.load(file)
    return data

def update_redu_with_ott_ids(redu_filename, matches_dict, output_filename):
    """Update the REDU file with OTT IDs based on the matches and save to a new CSV file."""
    df_redu = pd.read_csv(redu_filename, low_memory=False)
    #df_redu[['NCBI_prefix', 'tax_name']] = df_redu['NCBITaxonomy'].str.split('|', expand=True)
    df_redu.drop('NCBI_prefix', axis=1, inplace=True)

    names = df_redu['tax_name'].dropna()
    matches = matches_dict["results"]

    ott_id_map = {}
    for item in matches:
        name = item['name']
        for match in item['matches']:
            if 'ott_id' in match['taxon']:
                ott_id_map[name] = int(match['taxon']['ott_id'])  

    df_redu['OTT_ID'] = df_redu['tax_name'].map(ott_id_map).astype(pd.Int64Dtype()) 
    df_redu['ott_id_use'] = df_redu.groupby(['tax_name', 'OTT_ID']).OTT_ID.transform(lambda x: 'ott' + ''.join(map(str, x.unique())))
    df_redu.to_csv(output_filename, index=False)
    print(f"Updated data saved to {output_filename}")

def main():
    parser = argparse.ArgumentParser(description="Update the REDU file with OTT IDs from JSON matches.")
    parser.add_argument('--json_path', type=str, required=True, help="Path to the JSON file containing matches.")
    parser.add_argument('--redu_path', type=str, required=True, help="Path to the REDU file to update.")
    args = parser.parse_args()

    matches_dict = load_matches(args.json_path)
    output_filename = 'redu_updated.csv'
    update_redu_with_ott_ids(args.redu_path, matches_dict, output_filename)

if __name__ == "__main__":
    main()
