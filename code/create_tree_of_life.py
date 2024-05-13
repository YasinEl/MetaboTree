import argparse
import pandas as pd
from ete3 import Tree

def load_and_process_data(input_csv):
    """Load the CSV file and process the data to filter and create ott_id_use."""
    df = pd.read_csv(input_csv, dtype={'OTT_ID': 'Int64'})
    df = df.dropna(subset=['tax_name', 'OTT_ID'])
    return df

def load_and_filter_tree(tree_path, ott_ids):
    """Load the tree file and filter it based on OTT IDs."""
    tree = Tree(tree_path, format=1)  # format=1 assumes Newick format
    print("Available tip labels in the tree:", tree.get_leaf_names())  # Log tree tip labels for debugging
    # Find which OTT IDs are present in the tree
    valid_ids = set(ott_ids) & set(tree.get_leaf_names())
    print("OTT IDs to be kept:", valid_ids)  # Log valid OTT IDs
    tree.prune(list(valid_ids), preserve_branch_length=True)
    return tree

def save_tree(tree, output_path):
    """Save the filtered tree in Newick format."""
    tree.write(format=1, outfile=output_path)

def main(args):
    # Process data
    df_redu_taxas = load_and_process_data(args.input_csv)
    unique_ott_ids = df_redu_taxas['ott_id_use'].unique()

    # Load and filter the tree
    tree = load_and_filter_tree(args.tree_path, unique_ott_ids)

    # Output the tree and save it
    # print(tree)
    save_tree(tree, 'tree.nw')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process REDU data and filter phylogenetic tree based on OTT IDs.")
    parser.add_argument('input_csv', type=str, help="Path to the input CSV file containing REDU data.")
    parser.add_argument('tree_path', type=str, help="Path to the input tree file.")
    args = parser.parse_args()

    main(args)