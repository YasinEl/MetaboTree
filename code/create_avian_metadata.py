import argparse
import pandas as pd
from ete3 import Tree

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This script filters and updates a table based on a phylogenetic tree.")
    parser.add_argument("tree", help="Path to the phylogenetic tree file (in Newick format).")
    parser.add_argument("tsv", help="Path to the TSV file containing metadata.")
    args = parser.parse_args()

    # Load the metadata table
    df_metadata = pd.read_csv(args.tsv, sep="\t")

    # Ensure NCBI IDs are numeric
    df_metadata["NCBI"] = pd.to_numeric(df_metadata["NCBI"], errors="coerce")
    df_metadata.dropna(subset=["NCBI"], inplace=True)  # Drop rows with invalid NCBI IDs
    df_metadata["NCBI"] = df_metadata["NCBI"].astype(int)  # Convert to integers

    # Prefix "ncbi" to NCBI column values
    df_metadata["NCBI"] = "ncbi" + df_metadata["NCBI"].astype(str)

    # Load the phylogenetic tree
    tree = Tree(args.tree, format=1)

    # Add "ncbi" prefix to all tree node names (tips and internal) if the node name does not already start with ncbi
    for node in tree.traverse():
        if not node.name.startswith("ncbi"):
            node.name = "ncbi" + node.name

    # Get the list of tree tip labels
    tree_tips = list(tree.get_leaf_names())

    # Filter the table to keep only rows with NCBI values in the tree tips
    df_filtered = df_metadata[df_metadata["NCBI"].isin(tree_tips)].copy()

    # Find tips in the tree that are not in the metadata table
    tree_tips_set = set(tree_tips)
    table_tips_set = set(df_metadata["NCBI"])
    missing_tips = tree_tips_set - table_tips_set

    # Create new rows for missing tips with metabolomics_presence set to False
    new_rows = pd.DataFrame({
        "NCBI": list(missing_tips),
        "metabolomics_presence": [False] * len(missing_tips)
    })

    # Append the new rows to the filtered table
    df_updated = pd.concat([df_filtered, new_rows], ignore_index=True)

    # Rename columns: NCBI to FeatureID, FeatureID to otoid
    if "FeatureID" in df_updated.columns:
        df_updated.rename(columns={"FeatureID": "otoid", "NCBI": "FeatureID"}, inplace=True)
    else:
        df_updated.rename(columns={"NCBI": "FeatureID"}, inplace=True)
        df_updated["otoid"] = None  # Add an empty column for otoid if not in the original table

    # Reorder columns to ensure FeatureID is the first column
    column_order = ["FeatureID"] + [col for col in df_updated.columns if col != "FeatureID"]
    df_updated = df_updated[column_order]

    # Ensure FeatureID is character
    df_updated["FeatureID"] = df_updated["FeatureID"].astype(str)

    # Check for duplicates and print how many exist
    num_duplicates = df_updated.duplicated(subset=["FeatureID"]).sum()
    if num_duplicates > 0:
        print(f"Warning: {num_duplicates} duplicate FeatureIDs found in the metadata table.")

    # Check which nodes are in the tree but not the table
    missing_nodes = tree_tips_set - set(df_updated["FeatureID"])
    if missing_nodes:
        print(f"Warning: {len(missing_nodes)} nodes found in the tree but not in the table.")

    # Check which nodes are in the table but not the tree
    missing_nodes = set(df_updated["FeatureID"]) - tree_tips_set
    if missing_nodes:
        print(f"Warning: {len(missing_nodes)} nodes found in the table but not in the tree.")

    # Fill metabolomics_presence
    df_updated["metabolomics_presence"] = df_updated["metabolomics_presence"].fillna(True)

    # Save the updated table back to a TSV file
    output_file = "metadata_avian.tsv"
    df_updated.to_csv(output_file, sep="\t", index=False)

    # Save the updated tree back to a Newick file
    tree.write(format=1, outfile='tree_avian.nw')

    print(f"Updated table saved to {output_file}.")
    print(f"Updated tree saved to tree_avian.nw.")
