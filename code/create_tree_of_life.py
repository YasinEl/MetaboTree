import argparse
import pandas as pd
from ete3 import Tree

def load_and_process_data(input_csv):
    """Load the CSV file and process the data to filter and create ott_id_use."""
    df = pd.read_csv(input_csv, dtype={'uid': 'Int64'}, low_memory=False)
    df = df.dropna(subset=['uid'])
    return df

def load_and_filter_tree(tree_path, ott_ids):
    """Load the tree file and filter it based on OTT IDs."""
    tree = Tree(tree_path, format=1)  # format=1 assumes Newick format
    # Convert ott_ids to strings for comparison with leaf names
    ott_ids_str = set(map(str, ott_ids))
    #tree_leaves = set(tree.get_leaf_names())
    all_nodes = set(node.name for node in tree.traverse())

    print('tree_leaves')
    print(list(all_nodes)[1:10])
    # Find which OTT IDs are present in the tree
    valid_ids = ott_ids_str & all_nodes
    invalid_ids = ott_ids_str - all_nodes
    if not valid_ids:
        with open('invalid_ott_ids.txt', 'w') as f:
            for ott_id in invalid_ids:
                f.write(f"{ott_id}\n")
        raise ValueError("No valid OTT IDs found in the tree. Check invalid_ott_ids.txt for the list of invalid OTT IDs.")
    tree.prune(list(valid_ids), preserve_branch_length=True)
    return tree

def save_tree(tree, output_path):
    """Save the filtered tree in Newick format."""
    tree.write(format=1, outfile=output_path)

def main(args):
    # Process data
    df_redu_taxas = load_and_process_data(args.input_csv)
    unique_ott_ids = df_redu_taxas['uid_leaf'].dropna().unique()
    
    print('unique_ott_ids')
    print(unique_ott_ids[1:10])

    # Load and filter the tree
    tree = load_and_filter_tree(args.tree_path, unique_ott_ids)

    return tree



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process REDU data and filter phylogenetic tree based on OTT IDs.")
    parser.add_argument('--input_csv', type=str, help="Path to the input CSV file containing REDU data.")
    parser.add_argument('--tree_path', type=str, help="Path to the input tree file.")
    parser.add_argument('--usi', type=str, help="Path to the input tree file.")
    parser.add_argument('--cid', type=str, help="Path to the input tree file.")
    args = parser.parse_args()

    tree = main(args)

    print(args.usi)
    tree.write(format=1, outfile='treeGraph_' + args.usi.split(':')[4] + '_' + args.cid + '.nw')
    
    #fwrite(dt_tree_ids_redu[!is.na(uid_leaf) & uid_leaf != '', c('FeatureID', 'tax_name', 'detected', 'Cosine', 'Matching Peaks', 'UBERONBodyPartName', 'NCBIDivision')], paste0(c('treeAnnotation_',  lib_id, '_', cid, '.tsv'), collapse = ''), sep = '\t')
