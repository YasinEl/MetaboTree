import argparse
import pandas as pd
from ete3 import Tree
import re

def load_and_process_data(input_csv):
    """
    Load the CSV file and process the data to extract OTT IDs.
    Rows with missing 'uid' are dropped.
    """
    df = pd.read_csv(input_csv, dtype={'uid': 'Int64'}, low_memory=False)
    df = df.dropna(subset=['uid'])
    return df

def load_and_filter_tree(tree_path, ott_ids):
    """
    Load the tree (Newick format) and prune it to only keep nodes
    whose names are in the provided ott_ids (as strings).
    """
    tree = Tree(tree_path, format=1)
    ott_ids_str = set(map(str, ott_ids))
    all_nodes = set(node.name for node in tree.traverse())
    print("tree_leaves (sample of all node names):", list(all_nodes)[:10])
    valid_ids = ott_ids_str & all_nodes
    invalid_ids = ott_ids_str - all_nodes
    print("Number of valid OTT IDs:", len(valid_ids))
    print("Number of invalid OTT IDs:", len(invalid_ids))
    if not valid_ids:
        with open('invalid_ott_ids.txt', 'w') as f:
            for ott_id in invalid_ids:
                f.write(f"{ott_id}\n")
        raise ValueError("No valid OTT IDs found in the tree. See invalid_ott_ids.txt for details.")
    tree.prune(list(valid_ids), preserve_branch_length=True)
    print("After pruning, tree has", len(list(tree.traverse())), "nodes and", len(tree.get_leaves()), "leaves.")
    return tree

def save_tree(tree, output_path):
    """Write the tree to a file in Newick format."""
    tree.write(format=1, outfile=output_path)

def remove_internal_ott_nodes_simple(tree, prefix="ott"):
    """
    Traverse the pruned tree in postorder.
    
    For every child of a node that is internal and whose name matches the pattern
    (i.e. "ott" followed only by digits), perform the following:
    
      - If the internal OTT node has exactly one child:
           * For each grandchild, adjust its branch length (add the OTT node's branch length)
           * Remove the OTT node by promoting its sole child
           * Also create a leaf copy of the OTT node (to preserve its label)
           * Replace the OTT node with both the promoted child and the leaf copy.
      - Otherwise, leave the node unchanged.
      
    We then rebuild the parent's children list using add_child() to ensure proper parent pointers.
    """
    pattern = re.compile(rf"^{prefix}\d+$")
    # Traverse in postorder so that children are processed before their parent.
    for node in tree.traverse("postorder"):
        new_children = []
        for child in node.children:
            if pattern.match(child.name) and not child.is_leaf():
                if len(child.children) == 1:
                    # Promote the sole child.
                    for grandchild in child.children:
                        grandchild.dist += child.dist
                        new_children.append(grandchild)
                    # Also add a leaf copy of the OTT node.
                    leaf_copy = child.copy()
                    leaf_copy.children = []
                    new_children.append(leaf_copy)
                    print(f"Processed internal OTT node '{child.name}' at parent '{node.name}' (one child): promoted its child and kept a leaf copy.")
                else:
                    new_children.append(child)
            else:
                new_children.append(child)
        # Rebuild the children list to update parent pointers.
        node.children = []
        for ch in new_children:
            node.add_child(ch)
    # Process the root separately.
    if pattern.match(tree.name) and not tree.is_leaf():
        new_children = []
        for child in tree.children:
            child.dist += tree.dist
            new_children.append(child)
        leaf_copy = tree.copy()
        leaf_copy.children = []
        new_root = Tree(name="artificial_root", format=1)
        new_root.add_child(leaf_copy)
        for child in new_children:
            new_root.add_child(child)
        tree = new_root
        print("Processed root internal OTT node; created artificial root.")
    return tree

def main(args):
    # Process CSV to obtain unique OTT IDs from column 'uid_leaf'
    df_redu_taxas = load_and_process_data(args.input_csv)
    unique_ott_ids = df_redu_taxas['uid_leaf'].dropna().unique()
    print("unique_ott_ids (sample):", unique_ott_ids[:10])
    
    # Load and filter the tree based on these OTT IDs
    tree = load_and_filter_tree(args.tree_path, unique_ott_ids)
    
    # Remove internal OTT nodes only if they have one child.
    tree = remove_internal_ott_nodes_simple(tree, prefix="ott")
    
    return tree

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Process REDU data and filter phylogenetic tree based on OTT IDs, ensuring that internal OTT nodes with one child become leaves."
    )
    parser.add_argument('--input_csv', type=str, help="Path to the input CSV file containing REDU data.")
    parser.add_argument('--tree_path', type=str, help="Path to the input tree file.")
    parser.add_argument('--usi', type=str, help="Identifier string for output filename (e.g. mzspec:... ).")
    parser.add_argument('--cid', type=str, help="Identifier string for output filename.")
    args = parser.parse_args()

    tree = main(args)

    print("Output identifier:", args.usi)
    output_filename = 'treeGraph_' + args.usi.split(':')[4] + '_' + args.cid + '.nw'
    save_tree(tree, output_filename)
    print("Tree written to", output_filename)
