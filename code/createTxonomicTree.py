import argparse
from ete3 import PhyloTree, NCBITaxa, Tree
import pandas as pd
import json
import random

# Set up argument parsing
def parse_args():
    parser = argparse.ArgumentParser(description="Convert ete3 tree to JSON format.")
    parser.add_argument("--input_csv", required=True, help="Input CSV file with lineage data.")
    parser.add_argument("--output_json", required=True, help="Output JSON file for the tree.")
    parser.add_argument("--output_nw", required=True, help="Output Newick format tree file.")
    return parser.parse_args()

# Function to convert ete3 tree to json
def get_json(node):
    if not hasattr(node, 'evoltype'):
        dup = random.sample(['N','Y'], 1)[0]
    elif node.evoltype == "S":
        dup = "N"
    elif node.evoltype == "D":
        dup = "Y"

    node.name = node.name.replace("'", '')

    json_data = { "NCBI": node.name,
                  "duplication": dup,
                  "type": "node" if node.children else "leaf",
                }
    if node.children:
        json_data["children"] = []
        for ch in node.children:
            json_data["children"].append(get_json(ch))
    return json_data

def main():
    args = parse_args()
    df = pd.read_csv(args.input_csv)
    id_list = df['id'].tolist()
    ncbi = NCBITaxa()
    tree = ncbi.get_topology(id_list)
    tree.write(format=1, outfile=args.output_nw)
    redu_tree = Tree(args.output_nw, format=1)
    redu_tree_json = get_json(redu_tree)
    with open(args.output_json, 'w') as outfile:
        json.dump(redu_tree_json, outfile, indent=5)

if __name__ == "__main__":
    main()